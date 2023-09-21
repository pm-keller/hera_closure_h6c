# Bojan Nikolic <b.nikolic@mrao.cam.ac.uk> 2016,2017,2018 (as part of
# "heracasa" package), 2020 (as part of "closureps2")
"""
Closure quantities from visibility data
"""

import numpy


def rewrap(p):
    """
    Rewrap phase so that branch cut is along the negative real axis.
    """
    return numpy.arctan2(numpy.sin(p), numpy.cos(p))


def eitherWay(a1, a2, i, j):
    """
    Return rows where a1==i and a2==j OR a1==j and a2==i. Also return
    sign +1 if former or -1 if latter.
    """
    r1 = numpy.logical_and(a1 == i, a2 == j).nonzero()[0]
    if r1.shape[0]:
        return r1[0], +1.0
    else:
        r2 = numpy.logical_and(a1 == j, a2 == i).nonzero()[0]
        if r2.shape[0]:
            return r2[0], -1.0
        else:
            return 0, 0
        

def triadRows(a1, a2, tr):
    """
    Rows corresponding to single triad tr
    """
    i, j, k = tr
    p1, s1 = eitherWay(a1, a2, i, j)
    p2, s2 = eitherWay(a1, a2, j, k)
    p3, s3 = eitherWay(a1, a2, k, i)
    return ((p1, p2, p3), (s1, s2, s3))


def quadRows(a1, a2, qd):
    """
    Rows corresponding to a single quad qd
    """
    i, j, k, l = qd
    a12, s1 = eitherWay(a1, a2, i, j)
    a34, s2 = eitherWay(a1, a2, k, l)
    a13, s3 = eitherWay(a1, a2, i, k)
    a24, s4 = eitherWay(a1, a2, j, l)
    return (a12, a34, a13, a24)


def triads(a1, a2, alist):
    """
    List the rows corresponding to all triads in a list

    :param a1, a2: Arrays with antenna IDs for first and second antenna
    :param alist:  List of antenna IDs for which to generate the triads

    :returns: Tuple of (list of (tuple containing rows with data in the triad)),
              (list of tuples contianing antenna IDs in the triad),
              (list of signs to be used in computing a closure phase)

    """
    if len(alist) < 3:
        raise "Need at least three antennas to generate triads"
    nant = len(alist)
    rows = []
    tr = []
    signs = []
    for ni, i in enumerate(alist[:-2]):
        for nj, j in enumerate(alist[ni + 1 : -1]):
            for nk, k in enumerate(alist[ni + nj + 2 :]):
                ((p1, p2, p3), (s1, s2, s3)) = triadRows(a1, a2, (i, j, k))
                rows.append((p1, p2, p3))
                tr.append((i, j, k))
                signs.append((s1, s2, s3))
    return rows, tr, signs


def triadsList(a1, a2, trlist):
    """
    List the rows corresponding to specified triads

    :param a1, a2: Arrays with antenna IDs for first and second antenna

    :returns: Tuple of (list of (tuple containing rows with data in the triad)),
              (list of tuples contianing antenna IDs in the triad),
              (list of signs to be used in computing a closure phase)

    """
    rows = []
    trres = []
    signs = []
    for tr in trlist:
        ((p1, p2, p3), (s1, s2, s3)) = triadRows(a1, a2, tr)
        rows.append((p1, p2, p3))
        trres.append(tr)
        signs.append((s1, s2, s3))
    return rows, trres, signs


def closurePh(f, time, trlist=None, antf=None, wff=None):
    """The closure phase on a triad of antennas

    :param f: h5py data file

    :param time: time slot to read

    :param trlist: Explicit triads to compute the closure; if given
    only these triads will be returned

    :param antf: h5py antenna flags file

    :param wff: h5py waterfall flags file

    :returns: Dictionary: "phase": array containing phases; "tr": an
    array with the triad ids; "flags": array containing flags ,
    row-synchronous with phase; "JD": an array containing julian day of each integration

    """
    Ntimes = f["Header/Ntimes"][()]

    a1 = f["Header/ant_1_array"][time::Ntimes]
    a2 = f["Header/ant_2_array"][time::Ntimes]
    time_arr = f["Header/time_array"][time::Ntimes]
    ants = f["Header/antenna_numbers"]

    if trlist is not None:
        rows, tr, signs = triadsList(a1, a2, trlist)
    else:
        rows, tr, signs = triads(a1, a2, ants)

    ph = numpy.angle(
        f["Data/visdata"][time::Ntimes]["r"] + 1j * f["Data/visdata"][time::Ntimes]["i"]
    )
    fl = f["Data/flags"][time::Ntimes]

    jd = numpy.unique(f["Header/time_array"][time::Ntimes])
    if len(jd) > 1:
        raise RuntimeError("Inconsistent JD in integration")

    clp, flags = [], []
    for (p1, p2, p3), (s1, s2, s3) in zip(rows, signs):
        clp.append(rewrap(ph[p1] * s1 + ph[p2] * s2 + ph[p3] * s3))
        flags.append(fl[p1] | fl[p2] | fl[p3])

    flags = numpy.array(flags)

    trfl = []

    if antf is not None:
        afl = antf["Data/flag_array"][()]
        ants = antf["Header/ant_array"][()]
        times = antf["Header/time_array"][()]

        for (a1, a2, a3) in tr:
            i1 = numpy.where(ants == a1)[0]
            i2 = numpy.where(ants == a2)[0]
            i3 = numpy.where(ants == a3)[0]

            if len(i1) * len(i2) * len(i3) > 0:
                trfl.append(numpy.squeeze(afl[i1, :, time] | afl[i2, :, time] | afl[i3, :, time]))
            else:
                trfl.append(numpy.ones_like(afl[0, :, time]).astype(bool))

            print(numpy.shape(trfl[-1]), len(i1), len(i2), len(i3), len(i1) * len(i2) * len(i3))
    else:
        trfl = numpy.ones_like(flags[..., 0:2]).astype(bool)


    # get waterfall flags
    if wff is not None:
        wf = wff["Data/flag_array"][time]
    else:
        wf = numpy.ones((f["Header/Nfreqs"][()], f["Header/Npols"][()])).astype(bool)

    return {
        "phase": numpy.array(clp)[..., 0:2],
        "flags": flags[..., 0:2] | numpy.array(trfl),
        "wf": numpy.array(wf)[..., 0:2],
        "tr": numpy.array(tr),
        "JD": jd,
    }



def bispectrum(vis1, vis2, vis3, s1, s2, s3):
    """
    Return bispectrum of three visibilities
    """
    
    v1 = vis1.copy()
    v2 = vis2.copy()
    v3 = vis3.copy()
    
    if s1 == -1:
        v1 = v1.conjugate()
    if s2 == -1:
        v2 = v2.conjugate()
    if s3 == -1:
        v3 = v3.conjugate()
        
    return v1 * v2 * v3


def bispec(f, time, trlist=None, antf=None, wff=None):
    """The bispectrum on a triad of antennas

    :param f: h5py data file

    :param time: time slot to read

    :param trlist: Explicit triads to compute the closure; if given
    only these triads will be returned

    :param antf: h5py antenna flags file

    :param wff: h5py waterfall flags file

    :returns: Dictionary: "phase": array containing phases; "tr": an
    array with the triad ids; "flags": array containing flags ,
    row-synchronous with phase; "JD": an array containing julian day of each integration

    """
    Ntimes = f["Header/Ntimes"][()]

    a1 = f["Header/ant_1_array"][time::Ntimes]
    a2 = f["Header/ant_2_array"][time::Ntimes]

    ants = f["Header/antenna_numbers"]

    if trlist is not None:
        rows, tr, signs = triadsList(a1, a2, trlist)
    else:
        rows, tr, signs = triads(a1, a2, ants)

    vis = f["Data/visdata"][time::Ntimes]["r"] + 1j * f["Data/visdata"][time::Ntimes]["i"]
    fl = f["Data/flags"][time::Ntimes]

    # apply antenna flags
    if antf is not None:
        afl = antf["Data/flag_array"][()]
        ants = antf["Header/ant_array"][()]
        times = antf["Header/time_array"][()]
        ant_1_arr = f["Header/ant_1_array"][time::Ntimes]
        ant_2_arr = f["Header/ant_2_array"][time::Ntimes]
        time_arr = f["Header/time_array"][time::Ntimes]

        for i, (ant1, ant2, t) in enumerate(zip(ant_1_arr, ant_2_arr, time_arr)):
            i1 = numpy.where(ants == ant1)[0]
            i2 = numpy.where(ants == ant2)[0]
            it = numpy.where(numpy.isclose(times, t, atol=1e-7, rtol=0))[0]
            fl[i, :, 0:2] = fl[i, :, 0:2] | afl[i1, :, it] | afl[i2, :, it]

    # get waterfall flags
    if wff is not None:
        wf = wff["Data/flag_array"][time]
    else:
        wf = numpy.ones((f["Header/Nfreqs"][()], f["Header/Npols"][()])).astype(bool)

    jd = numpy.unique(f["Header/time_array"][time::Ntimes])
    if len(jd) > 1:
        raise RuntimeError("Inconsistent JD in integration")

    bisp, flags = [], []
    for (p1, p2, p3), (s1, s2, s3) in zip(rows, signs):
        bisp.append(bispectrum(vis[p1], vis[p2], vis[p3], s1, s2, s3))
        flags.append(fl[p1] | fl[p2] | fl[p3])

    return {
        "bispec": numpy.array(bisp)[..., 0:2],
        "flags": numpy.array(flags)[..., 0:2],
        "wf": numpy.array(wf)[..., 0:2],
        "tr": numpy.array(tr),
        "JD": jd,
    }