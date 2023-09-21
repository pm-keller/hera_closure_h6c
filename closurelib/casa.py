"""
Tools for working with CASA
"""

import os
import glob
import subprocess


def mk_casa_script(casa_str, path="casa_script.py"):
    """Generate a CASA script from a string

    Args:
        casa_str (list): list of strings to be combined in file
        path (str, optional): file to write casa_str to. Defaults to "casa_script.py".
    """
    f = open(path, "w")

    for s in casa_str:
        f.write(s)
        f.write("\n")

    f.close()


def del_casa_script(path="casa_script.py"):
    """Delete a CASA script

    Args:
        path (str, optional): path to CASA script. Defaults to "casa_script.py".
    """
    os.remove(path)


def run_casa(path="casa_script.py"):
    """Run a CASA script

    Args:
        path (str, optional): path to CASA script. Defaults to "casa_script.py".
    """
    subprocess.call(
        "casa --nologger --nologfile -c {:s}".format(os.path.abspath("casa_script.py")),
        shell=True,
        stdout=subprocess.PIPE,
    )

    for e in ["*.log", ".last"]:
        for filename in glob.glob(e):
            os.remove(filename)


def run(casa_str):
    """Run CASA on temporary script which is created and subsequently deleted

    Args:
        casa_str (list): list of strings constituting a CASA script
    """

    # generate CASA script
    mk_casa_script(casa_str)

    # run CASA on casa script
    run_casa()

    # delete CASA script
    del_casa_script()


def uvfitstoms(uvfits_file, mspath):
    """Convert UVFITS file to a CASA measurement set

    Args:
        uvfits_file (str): UVFITS file
        mspath (str): path to CASA measurement set
    """

    print("\nconverting UVFITS file {:s}".format(uvfits_file))
    print("to CASA measurement set {:s}\n".format(mspath))

    casa_str = ["importuvfits(fitsfile='{:s}', vis='{:s}')".format(uvfits_file, mspath)]
    run(casa_str)


def mkcl(clfile, clpath):
    """Generate a CASA component list

    Args:
        clfile (str): component list file
        clpath (str): path to CASA component list
    """

    print("\nconvert component list {:s}".format(clfile))
    print("to a CASA component list {:s}\n".format(clpath))

    casa_str = [
        "os.system('rm -rf ' + '{:s}')".format(clpath),
        "cl.done()",
        "execfile('{:s}')".format(clfile),
        "cl.rename('{:s}')".format(clpath),
        "cl.close()",
    ]

    run(casa_str)


def ft(mspath, clpath):
    """Run the CASA ft task (simulate visibilities)

    Args:
        mspath (str): path to CASA measurement set
        clpath (str): path to CASA component list
    """

    print("\nrunning CASA ft for component list {:s}".format(clpath))
    print("and writing the model into the MODEL_DATA column of")
    print("measurement set {:s}".format(mspath))

    casa_str = [
        "ft(vis='{:s}', complist='{:s}', usescratch=True) ".format(mspath, clpath),
    ]

    run(casa_str)
