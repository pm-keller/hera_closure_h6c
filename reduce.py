""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Initial Data Reduction

"""

import os
import h5py
import argparse
import numpy as np

from closurelib import cptools as cp


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inpath", help="Path to closure data.", type=str)
parser.add_argument(
    "-o", "--outpath", help="Path to write reduced closure data to.", type=str
)
parser.add_argument("--fmin", help="Minimum Frequency (MHz).", type=float)
parser.add_argument("--fmax", help="Maximum Frequency (MHz).", type=float)
args = parser.parse_args()

# load data
with h5py.File(args.inpath, "r") as f:
    jd = f["JD"][()]
    lst = f["LST"][()]
    trlist = f["triads"][()]
    frq = f["FRQ"][()]

# shape of reduced data set
Nlst = len(lst)

if os.path.exists(args.outpath):
    os.remove(args.outpath)
fw = h5py.File(args.outpath, "a")

for t in range(Nlst):
    with h5py.File(args.inpath, "r") as f:
        phase = f["phase"][:, t]
        flags = f["flags"][:, t]
        wf = f["wf"][:, t]

    # update flags
    wf = np.repeat(np.expand_dims(wf, axis=1), flags.shape[1], axis=1)
    flags = flags | wf | np.isnan(phase)
    phase[flags] = np.nan

    # select subband and rearrange axes in order (polarisation, JD, triad, LST, Frequency)
    phase, sb = cp.select(phase, frq, args.fmin, args.fmax, axis=2)
    phase = np.moveaxis(phase, (0, 1, 2, 3), (1, 2, 3, 0))

    # remove nan-only slices
    #phase, jd_idx = cp.remove_nan_slices(phase, axis=1, return_index=True)
    #phase, tr_idx = cp.remove_nan_slices(phase, axis=2, return_index=True)

    # add lst axis
    phase = phase[:, :, :, np.newaxis]
    

    if t == 0:
        #jd = jd[jd_idx]
        #trlist = trlist[tr_idx]
        shape = (2, len(jd), len(trlist), Nlst, len(sb))

        fw.create_dataset("JD", data=jd)
        fw.create_dataset("LST", data=lst)
        fw.create_dataset("FRQ", data=sb)
        fw.create_dataset("triads", data=trlist)
        fw.create_dataset("phase", data=phase, chunks=True, maxshape=shape)
    else:
        fw["phase"].resize((fw["phase"].shape[3] + 1), axis=3)
        fw["phase"][:, :, :, -1:] = phase

fw.close()
