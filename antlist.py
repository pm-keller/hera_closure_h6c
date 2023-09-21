""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Generate a list of antennas which can be used for the analysis of a field

"""

import sys
sys.path.append("/users/pkeller/code/H1C_IDR3.2/")

import numpy as np
import h5py
import os
import shutil

from pyuvdata import UVData

from closurelib import libtools as librarian

delta = 9.690134587635768 / 3600
lstdata_list = [{"start": 4.99, "stop": 5.01, "delta": delta, "name": "5h"}]

indir = "/lustre/aoc/projects/hera/pkeller/data/H6C/alignment/"
outdir = "/lustre/aoc/projects/hera/pkeller/data/array/"
stagedir = "/lustre/aoc/projects/hera/pkeller/data/H6C/tmp/"

for lstdata in lstdata_list:
    alignfname = os.path.join(indir, f"h6c_idr1_{lstdata['name']}.npz")
    fg = np.load(alignfname)["fnameg"].astype(str)

    print(fg[0, 0][:-2])
    
    success, fpath = librarian.stageFile(fg[0, 0], stagedir)
    print(fpath)
    UV = UVData()
    UV.read(fpath, file_type="uvh5")
    a = UV.ant_1_array
    os.remove(fpath)

    np.savetxt(os.path.join(outdir, f"antlist_h1c_{lstdata['name']}.dat"), a)
