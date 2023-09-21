#!/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python
#SBATCH -j align
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH -o align.out
#SBATCH -l walltime=1:00:00

""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

LST alignment

"""

import sys
sys.path.append("/users/pkeller/code/H6C/")

import os
import pandas
from closurelib import heradata


indir = "/lustre/aoc/projects/hera/pkeller/data/H6C/scans/"
outdir = "/lustre/aoc/projects/hera/pkeller/data/H6C/alignment/"
lstdata = {"start": 4.995, "stop": 5.005, "delta": 9.690134650097804, "name": "5h"}

d = pandas.read_csv(os.path.join(indir, "h6c_idr1_scan.csv"))
alignfname = os.path.join(outdir, f"h6c_idr1_{lstdata['name']}.npz")
heradata.nearestTimeF(alignfname, d, [lstdata["start"], lstdata["stop"]], lstdata["delta"])