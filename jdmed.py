""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Averaging over JD's

"""

import numpy as np
import argparse
import h5py

from multiprocessing import Pool

from closurelib import cptools as cp

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to closure data.", type=str)
parser.add_argument("-m", "--model", help="use model data.", action="store_true")
args = parser.parse_args()

if args.model:
    ifmodel = " model"
else:
    ifmodel = ""

for pol, polname in zip(range(2), ["XX", "YY"]):
    # load data
    with h5py.File(args.path, "r") as f:
        jd = f["JD"][()]
        lst = f["LST"][()]
        trlist = f["triads"][()]
        phase = f[f"phase{ifmodel}"][pol]
        nanflags = np.isnan(f["phase"][pol])

    Nlst = len(lst)

    # complex exponential of closure phase
    eicp = np.exp(1j * phase)
    del phase

    # apply flags and remove NaN-only slices
    eicp[nanflags] = np.nan
    eicp, jdidx = cp.remove_nan_slices(eicp, axis=0, return_index=True)
    eicp, tridx = cp.remove_nan_slices(eicp, axis=1, return_index=True)
    eicp = np.ma.array(eicp, mask=np.isnan(eicp))
    jd = jd[jdidx]
    trlist = trlist[tridx]
    
    # average complex closure phase in 2 and 4 JD bins respectively
    def geomed1(idx):
        return np.array([cp.geomed(eicp[i::2, :, idx], axis=(0,)) for i in range(2)])

    def geomed2(idx):
        return np.array([cp.geomed(eicp[i::4, :, idx], axis=(0,)) for i in range(4)])
    
    with Pool(processes=16) as pool:
        eicp1 = pool.map(geomed1, np.arange(Nlst))
    with Pool(processes=16) as pool:
        eicp2 = pool.map(geomed2, np.arange(Nlst))

    eicp1 = np.moveaxis(eicp1, 0, 2)
    eicp2 = np.moveaxis(eicp2, 0, 2)

    # write to file
    f = h5py.File(args.path, "a")

    for eicp, nbins in zip([eicp1, eicp2], [2, 4]):
        dname = f"eicp jdmed new ({nbins}) {polname}{ifmodel}"
        if dname in f.keys():
            del f[dname]
        f.create_dataset(dname, data=eicp)
    
    if f"triads {polname}" in f.keys():
        del f[f"triads {polname}"]
    f.create_dataset(f"triads {polname}", data=trlist)
    f.close()