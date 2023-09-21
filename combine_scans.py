""" 

Pascal M. Keller <pmk46@mrao.cam.ac.uk> 2021/22
Cavendish Astrophysics, University of Cambridge, UK

Combine all scans to one file.

"""


import sys

sys.path.append("/users/pkeller/code/H6C/")

import os
import pandas
import numpy as np
from closurelib import heradata


# meta data
jdrange = (2459860, 2459877)
jdlist = np.arange(*jdrange)
outdir = "/lustre/aoc/projects/hera/pkeller/data/H6C/scans/"

# remove bad day
#jdlist = np.hstack([jdlist[:171], jdlist[172:]])

### COMBINE SCANS ####

masterf = os.path.join(outdir, "h6c_idr1_scan.csv")

if os.path.exists(masterf):
    os.remove(masterf)

fname = os.path.join(outdir, f"{str(jdlist[0])}_scan.csv")
data = pandas.read_csv(fname)
data = data.loc[:, ~data.columns.str.contains("^Unnamed")]
data["fname"] = [os.path.basename(path) for path in data["fname"].values]
data.to_csv(masterf)

for jd in jdlist[1:]:
    fname = os.path.join(outdir, f"{str(jd)}_scan.csv")
    df = pandas.read_csv(fname)
    masterdf = pandas.read_csv(masterf)
    data = pandas.concat([df, masterdf], ignore_index=True)
    data = data.loc[:, ~data.columns.str.contains("^Unnamed")]
    data["fname"] = [os.path.basename(path) for path in data["fname"].values]
    data.to_csv(masterf)