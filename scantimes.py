#!/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python
#SBATCH -j scantimes
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH -o scantimes.out
#SBATCH -l walltime=48:00:00

""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Scan all visibility files and save metadata to file.

"""

import sys

sys.path.append("/users/pkeller/code/H6C/")

import os
import pandas
import numpy as np
from closurelib import heradata

jdrange = (2459860, 2459877)
jdlist = np.arange(*jdrange)
lstdata = {"start": 4.995, "stop": 5.005, "delta": 9.690134650097804 / 3600}
lstrange = (lstdata["start"], lstdata["stop"])
outdir = "/lustre/aoc/projects/hera/pkeller/data/H6C/scans/"

for i, jd in enumerate(jdlist):
    with open("./printfile.out", "w") as f:
        f.write(f"processing Julian Day {jd}, {i} of {len(jdlist)}\n")

    d = heradata.scanTimeLib((jd, jd + 1), lstrange, outdir)
    d.to_csv(os.path.join(outdir, f"{str(jd)}_scan.csv"))
