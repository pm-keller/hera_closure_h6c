#!/lustre/aoc/projects/hera/pkeller/anaconda3/envs/closure_analysis/bin/python3
#SBATCH -p hera
#SBATCH -J mkcl
#SBATCH -o mkcl.out
#SBATCH -t 12:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Make closure phase or bispectrum data sets.

"""

import sys

sys.path.append("/users/pkeller/code/H6C/")

import os
import glob
import pandas
import numpy as np
from closurelib import heradata
from closurelib import triads

# time delta [s]
delta = 9.690134650097804 / 3600

# list of HERA fields
lstdata_list = [{"start": 4.995, "stop": 5.005, "delta": delta, "name": "5h"}]

# list of triad classes
triad_list = ["EQ14",]

# select triad class and HERA field
triad = triad_list[0]
lstdata = lstdata_list[0]

# generate a name
name = f"{triad}_{lstdata['name']}"

# specify i/o directories 
outdir = "/lustre/aoc/projects/hera/pkeller/data/H6C/"
fldir = "/lustre/aoc/projects/hera/h6c-analysis/IDR2/"
stagedir = os.path.join(outdir, "tmp")
alignfname = os.path.join(outdir, "alignment", f"h6c_idr1_{lstdata['name']}.npz")

# generate a triad list
antfile = f"/lustre/aoc/projects/hera/pkeller/data/array/antlist_h6c_{lstdata['name']}.dat"
trfile = f"/lustre/aoc/projects/hera/pkeller/data/array/trlist_{triad}.dat"
trlist = np.loadtxt(trfile).astype(int)
antlist = np.loadtxt(antfile)
trlist = triads.trlistSub(trlist, antlist)

# prepare closure phase data file
datafname = os.path.join(outdir, "sample", f"{name}.h5")

#os.remove(datafname)
if not os.path.isfile(datafname):
    heradata.prepClosureDS(alignfname, trlist, datafname)

# add aligned closure phases to data file
fg = np.load(alignfname)["fnameg"]
print("align, shape: ", fg.shape)

for j in range(fg.shape[0]):
    with open("/users/pkeller/code/H6C/printfile.out", "w") as f:
        print(f"processing day {j} of {fg.shape[0]}\n")
        f.write(f"processing day {j} of {fg.shape[0]}\n")

    for k in range(fg.shape[1]):
        print(f"{fg[j, k]}")
        heradata.addClosureDS(alignfname, outdir, datafname, trlist, j, k, fldir)

    for f in glob.iglob(r"/lustre/aoc/projects/hera/pkeller/data/H6C/tmp/*", recursive=True):
        if os.path.isfile(f):
            os.remove(f)
            
