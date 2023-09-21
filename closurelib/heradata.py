# Bojan Nikolic <b.nikolic@mrao.cam.ac.uk> 2020
# Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021
"""
Handling of HERA datasets in uvh5 and similar format
"""

import numpy
import pandas
import h5py
import glob
import shutil
import os
import re

from astropy.time import Time
from pyuvdata import UVData

from . import clquants
from . import libtools as librarian


def jd_to_lst(jd, lat=-30.72138329631366, lon=21.428305555555557):
    """Convert julian day to apparent local sidereal time
    defaults are HERA coordinates

    Args:
        jd (float): julian day
        lat (float, optional): latitude. Defaults to -30.72138329631366.
        lon (float, optional): longitude. Defaults to 21.428305555555557.

    Returns:
        float: local sidereal time
    """

    t = Time(jd, format="jd", location=(lon, lat))

    return t.sidereal_time("apparent").to_value()


def LSTgrid(lstdata):
    """Generate an evenly spaced LST grid

    Args:
        lstdata (dict): {"start", "stop", "delta"} in hours

    Returns:
        array: LST grid
    """

    if lstdata["start"] < lstdata["stop"]:
        grid = numpy.arange(lstdata["start"], lstdata["stop"], lstdata["delta"])
    else:
        g1 = numpy.arange(lstdata["start"], 24, lstdata["delta"])
        g2 = numpy.arange(
            lstdata["delta"] - (24 - g1[-1]), lstdata["stop"], lstdata["delta"]
        )
        grid = numpy.hstack([g1, g2])
    return grid


def filesToJD(fnames, unique=True):
    """Extract julian days from a list of file names

    Args:
        fnames (list or str): file name or list of file names
        unique (bool): if True, return a unique and sorted list of julian days

    Returns:
        list: sorted list of julian days
    """

    if type(fnames) is str:
        jd = re.findall(r"\d+", fnames)[0]
    else:
        jd = [re.findall(r"\d+", fname)[0] for fname in fnames]

    if unique:
        return numpy.sort(numpy.unique(jd)).tolist()
    else:
        return jd
    

def filesToLST(fnames):
    """Extract local sideral times from a list of file names of LST-binned data

    Args:
        fnames (list or str): file name or list of file names

    Returns:
        list: sorted list of local sidereal times
    """

    if type(fnames) is str:
        lst = float(re.findall(r"\d+\.\d+", fnames)[0])
    else:
        lst = [float(re.findall(r"\d+\.\d+", fname)[0]) for fname in fnames]

    return numpy.array(lst) * 12 / numpy.pi


def frqFromFile(fpath):
    """Extract frequency array from a data file

    Args:
        fpath (str): path to data file

    Returns:
        array: frequencies
    """

    with h5py.File(fpath, "r") as fin:
        frq = fin["Header/frequency_array"][()]

    return frq


def getDataLST(root, ant1, ant2, lst, filename="*"):
    """Get LST-binned visibility data, given two antennas and an LST
    
    Args:
        root (str): path to LST-binned data
        ant1 (int): antenna number 1
        ant2 (int): antenna number 2
        lst (float): local sidereal time
        filename (str): filenames to search for
        
    Returns:
        numpy array: LST-binned visibilities of two antennas at the given LST.
    """
    
    # get file names
    fnames = numpy.array(glob.glob(root + "/" + filename))
    
    # find file containing LST
    lst_list = numpy.array(filesToLST(fnames))
    idx = numpy.where((lst_list - lst) <= 0)[0]
    fname = fnames[idx][numpy.argmax(lst_list[idx])]
        
    # read data 
    uv = UVData()
    uv.read(fname)
    uv_data = uv.get_data(ant1, ant2)
    
    
    lst_array = numpy.unique(uv.lst_array * 12 / numpy.pi)
    idx = numpy.argmin(numpy.abs(lst_array - lst))
    
    return uv_data[idx]


def scanTimeLib(jdrange, lstrange, outdir, namematch="zen.%.sum.uvh5", jdex=[]):
    """Scan librarian for data files of a given JD and LST range and extract metadata

    Args:
        jdrange (tuple): (min, max) julian day
        lstrange (tuple): (min, max) local sidereal time
        outdir (str): directory to write temporary data to
        namematch (str, optional): file name must match this string. Defaults to "zen.*.sum.uvh5"
        jdex (list, optional): julian days to exclude

    Raises:
        Exception: Failed to stage file from librarian

    Returns:
        DataFrame: metadata
    """

    fname = []
    jd = []
    lsts = []

    fnames = librarian.getFileNames(jdrange, lstrange, namematch)
    stagedir = os.path.join(outdir, "tmp")
    if not os.path.exists(stagedir):
        os.mkdir(stagedir)

    for i, fn in enumerate(fnames):
        with open("./printfile.out", "a") as f:
            f.write(f"processing file {i} of {len(fnames)}\n")

        success, fpath = librarian.stageFile(fn, stagedir)
        if not success:
            raise Exception(f"Failed to stage file {fn} from librarian!")
                
        if os.path.exists(fpath):
            UV = UVData()
            UV.read(fpath)
            ll = numpy.unique(UV.lst_array)
            fname += ([str(fpath)] * len(ll))
            jd += ([numpy.floor(UV.time_array).astype(int)[0]] * len(ll))
            lsts += list(ll * 12 / numpy.pi)

            os.remove(fpath)

    return pandas.DataFrame({"fname": fname, "jd": jd, "lsts": lsts})


def nearestTime(d, lstrange, delta=9.690134650097804):
    """Find **nearest** observation for a grid of LSTs

    Args:
        d (DataFrame): metadata as returned by scanTimeLib()
        lstrange (dict): [LST min, LST max] in hours

    Returns:
        Grid with filename containing closest integration; grid of integration lst, jd axis and lst grid axis
    """
    lg = numpy.arange(lstrange[0], lstrange[1], delta/3600)
    jg = d["jd"].unique()
    df = d.loc[(d["lsts"] > lstrange[0]-delta/3600) & (d["lsts"] < lstrange[1]+delta/3600)]
    dsize = [len(df.loc[df["jd"] == jd]) for jd in jg]

    jg = jg[numpy.where(dsize >= numpy.max(dsize)-1)]
    dt = numpy.array(d["fname"], dtype=numpy.string_).dtype
    res = numpy.zeros(shape=(len(jg), len(lg)), dtype=dt)
    dt = numpy.zeros(shape=(len(jg), len(lg)), dtype=numpy.double)
    
    for i, jd in enumerate(jg):
        dd = d[d.jd == jd].reset_index(drop=True)
        
        for k, lst in enumerate(lg):            
            mm = (dd.lsts-lst).abs().idxmin()
            
            if (dd.lsts-lst).abs().loc[mm] <= delta / 2:
                res[i,k] = dd.loc[mm].fname
                dt[i,k] = dd.loc[mm].lsts
            else:
                raise("LST discontinuity!") 
            
    return res, dt, jg, lg


def nearestTimeF(fnameout, *args):
    """Like nearestTime but save output to a npz file

    Args:
        fnameout (str): path to npz file
    """

    fg, lstg, jdaxis, lstaxis = nearestTime(*args)
    numpy.savez(fnameout, fnameg=fg, lstg=lstg, jdaxis=jdaxis, lstaxis=lstaxis)


def prepClosureDS(timesnpz, trlist, fnameout):
    """Prepare a HDF5 dataset for writing closure phase
    information. Separate prepare stage allows easy incremental
    filling in

    Args:
        timesnpz (dict):  information on LST/JD grid and nearest dataset for each point (see nearestTimeF)
        trlist (list): triads
        fnameout (str): path to HDF5 file
    """
    f = numpy.load(timesnpz)

    with h5py.File(fnameout, "w") as fout:
        gridshape = f["fnameg"].shape + (len(trlist), 1536, 2)
        gridshape_wf = f["fnameg"].shape + (1536, 2)
        fout.create_dataset("phase", shape=gridshape, dtype=numpy.float64)
        fout.create_dataset("flags", shape=gridshape, dtype=numpy.bool)
        fout.create_dataset("wf", shape=gridshape_wf, dtype=numpy.bool)
        fout.create_dataset("FRQ", shape=(1536,), dtype=numpy.float64)
        fout.create_dataset("triads", data=trlist)
        fout.create_dataset("JD", data=f["jdaxis"])
        fout.create_dataset("LST", data=f["lstaxis"])


def addClosureDS(timesnpz, outdir, fnameout, trlist, jdi, lsti, fldir=None, bisp=False):
    """Fill in closure phase information for JD index jdi and lst index lsti

    Args:
        timesnpz (dict): information on LST/JD grid and nearest dataset for each point (see nearestTimeF)
        outdir (str): directory to write temporary data to
        fnameout (str): path to HDF5 file
        trlist (list): triads
        jdi (int): JD index
        lsti (int): LST index
        fldir (str): Directory containing flagging files
        bisp (bool): If true, compute bispectrum instead of phase.
    """
    f = numpy.load(timesnpz)

    stagedir = os.path.join(outdir, "tmp")
    if not os.path.exists(stagedir):
        os.mkdir(stagedir)

    with h5py.File(fnameout, "a") as fout:
        fname = f["fnameg"][jdi, lsti].decode('utf-8')
        success, fpath = librarian.stageFile(fname, stagedir)
        
        if not success:
            raise Exception(f"Failed to stage file {fname} from librarian!")

        ff = h5py.File(fpath, "r")
        ffjd = numpy.unique(ff["Header/time_array"])
        fflst = jd_to_lst(ffjd)
        dlst = numpy.abs(fflst - fout["LST"][lsti])
        fflst[fflst > 12] = 24 - fflst[fflst > 12]
        ftime = dlst.argmin()
        jd = str(int(ffjd[0]))

        antfpath = os.path.join(fldir, jd, fname[:-4] + "antenna_flags.h5")
        wfpath = os.path.join(fldir, jd, fname[:-4] + "flag_waterfall.h5")
        
        if os.path.exists(antfpath):
            antf = h5py.File(antfpath)
        else:
            antf = None

        if os.path.exists(wfpath):
            wff = h5py.File(wfpath)
        else:
            wff = None

        if bisp:
            r = clquants.bispec(ff, ftime, trlist=trlist, antf=antf, wff=wff)
            fout["phase"][jdi, lsti] = r["bispec"]
        else:
            r = clquants.closurePh(ff, ftime, trlist=trlist, antf=antf, wff=wff)
            fout["phase"][jdi, lsti] = r["phase"]
        
        fout["flags"][jdi, lsti] = r["flags"]
        fout["wf"][jdi, lsti] = r["wf"]
        fout["FRQ"] = ff["Header/freq_array"]

        os.remove(fpath)


def addTrList(fnameout, trlist):
    """
    Insert triad information
    """
    with h5py.File(fnameout, "a") as fout:
        fout["triads"]=numpy.array(trlist)

def addFREQList(fnameout, frq):
    """
    Insert frequency information
    """
    with h5py.File(fnameout, "a") as fout:
        fout["FRQ"]=numpy.array(frq)
        
def prllSubset(fnamein, fnameout):
    """Extract only parallel polarisation products

    Args:
        fnamein (str): path to closure phase data file
        fnameout (str): path to wirte new data to
    """
    fin = h5py.File(fnamein, "r")

    with h5py.File(fnameout, "w") as fout:
        fout.create_dataset("phase", data=fin["phase"][..., 0:2])

        for ff in ["flags", "JD", "LST", "FRQ"]:
            fout.create_dataset(ff, data=fin[ff])


def visFormat(fnamein, fnameout):
    """Re-orient for best visualisation using napari

    Args:
        fnamein (str): path to closure phase data file
        fnameout (str): path to wirte new data to
    """
    fin = h5py.File(fnamein, "r")

    with h5py.File(fnameout, "w") as fout:
        z = fin["phase"].shape
        fout.create_dataset("phase", shape=(z[4], z[0], z[2], z[1], z[3]))

        for day in range(z[0]):
            d = numpy.moveaxis(fin["phase"][day], 3, 0)
            d = numpy.moveaxis(d, 2, 1)
            fout["phase"][:, day, :, :, :] = d + numpy.pi

        for ff in ["flags", "JD", "LST", "FRQ"]:
            fout.create_dataset(ff, data=fin[ff])
