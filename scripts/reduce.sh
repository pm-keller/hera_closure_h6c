#!/bin/bash
#SBATCH -p hera
#SBATCH -J reduce
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
        inpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}.h5"
        
        band=1
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        
        if [ -f "${inpath}" ]
        then
            "${pythonpath}" reduce.py -i "${inpath}" -o "${outpath}" --fmin 50.0946044896875 --fmax 62.3016357446875
        fi

        band=2
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        
        if [ -f "${inpath}" ]
        then
            "${pythonpath}" reduce.py -i "${inpath}" -o "${outpath}" --fmin 70.1141357396875 --fmax 87.4481201196875
        fi

        band=3
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        
        if [ -f "${inpath}" ]
        then
            "${pythonpath}" reduce.py -i "${inpath}" -o "${outpath}" --fmin 108.0780029271875 --fmax 124.6795654321875
        fi

        band=4
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        
        if [ -f "${inpath}" ]
        then
            "${pythonpath}" reduce.py -i "${inpath}" -o "${outpath}" --fmin 154.4647216771875 --fmax 175.0946044946875
        fi

        band=5
        outpath="/lustre/aoc/projects/hera/pkeller/data/H6C/sample/${trclass}_${field}_B${band}.h5"
        
        if [ -f "${inpath}" ]
        then
            "${pythonpath}" reduce.py -i "${inpath}" -o "${outpath}" --fmin 201.9500732396875 --fmax 223.0682373071875
        fi
    done
done