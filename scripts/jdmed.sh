#!/bin/bash
#SBATCH -p hera
#SBATCH -J jdmed
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk


pythonpath="/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3"

for trclass in EQ14
do 
    for field in 5h
    do  
        band=1
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        "${pythonpath}" jdmed.py -p "${outpath}"

        band=2
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        "${pythonpath}" jdmed.py -p "${outpath}"

        band=3
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        "${pythonpath}" jdmed.py -p "${outpath}"

        band=4
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        "${pythonpath}" jdmed.py -p "${outpath}"

        band=5
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        "${pythonpath}" jdmed.py -p "${outpath}"
    done
done
