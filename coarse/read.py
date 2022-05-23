from numpy import pi, ones
from numpy import array, zeros, dot
from numpy.linalg import norm, inv

import os

from quantumEspresso.read import readRawExport, getPList, getCA5, getPA5
from quantumEspresso.read import sortByKpt, fileSort
from constants import *

"""Get coefficients C^n_{G+k} and PW momenta from pw_export output

to do:
* accept arbitrary k-point dimensions
* automatically detect <dir>/<prefix>.export
"""

def readExport(
        exroot2 = 'tmp/hBN.export',
        exroot3 = 'tmp/hBN.export',
        ):
    """Returns arrays containing PW coefficients for a pair of QE runs.

    * k_a2 much coarser than k3_a2, i.e., nk2 < nk3
    * k_a3 should contain k_a2
    * prob_a6[k2x, k2y, k3x, k3y, v, c] --> prob of |v,k2> --> |c,k3>
    * must divide by 2 areas corresponding to spacing of k2_a3 and k3_a3
    * assumes uniform k-point mesh. No reduction by symmetry!
    * assumes runs only differ in number of k-points

    Reads from pw_export.x output files wfc.{kpt}, grid.{kpt}
    * run in directory containing outdir
    * if k-points lie on mesh
      * p_a6[kx, ky, Kx, Ky, Kz, mu] --> p^mu (eV)
      * C_a6[kx, ky, n, Kx, Ky, Kz] --> C^n_{K+k}
    * else 
      * p_a5[k, Kx, Ky, Kz, mu] --> p^mu (eV)
      * C_a5[k, n, Kx, Ky, Kz] --> C^n_{K+k}
    * NOTE: K indices increase from [-Kmax, -Kmax+1, ..., Kmax]
            k indices increase from [0, 1, 2, ..., kdim]
            allows for things like q = (k2 - k1)%kdim; Qq = (k2 - k1)//kdim

    exroot{2,3}: directory containing pw_export outputs (str)
        * usually {outdir}/{prefix}.export
    """
    # Read coarse and dense runs.  Only k_a2, E_a2, and nk{x,y} differ
    k2_a2, E2_a2, b_a2, nk2, nk2x, nk2y, nb, ne, volume, area, encut \
            = readRawExport(exroot2)
    k3_a2, E3_a2, b_a2, nk3, nk3x, nk3y, nb, ne, volume, area, encut \
            = readRawExport(exroot3)

    # find max Kx, Ky, and Kz so that all levels have same dimension
    maxKx, maxKy, maxKz = 0, 0, 0
    minKx, minKy, minKz = 0, 0, 0

    # record all RLVs listed in all coarse and dense grid files
    p2_list, minKx, minKy, minKz, maxKx, maxKy, maxKz = getPList(
                exroot2, k2_a2, b_a2, minKx, minKy, minKz, maxKx, maxKy, maxKz)
    p3_list, minKx, minKy, minKz, maxKx, maxKy, maxKz = getPList(
                exroot3, k3_a2, b_a2, minKx, minKy, minKz, maxKx, maxKy, maxKz)

    # array ranges
    nKx = maxKx - minKx + 1 
    nKy = maxKy - minKy + 1 
    nKz = maxKz - minKz + 1 

    # transfer momenta to an array
    p2_a5 = getPA5(p2_list, nk2, nKx, nKy, nKz)
    p3_a5 = getPA5(p3_list, nk3, nKx, nKy, nKz)

    # record all coefficients listed in all wfc files
    C2_a5 = getCA5(p2_list, exroot2, nk2, nb, nKx, nKy, nKz)
    C3_a5 = getCA5(p3_list, exroot3, nk3, nb, nKx, nKy, nKz)

    # treat autogenerated k-point meshes differently from k-point lines
    if nk2x > 0:
        C2_a6, p2_a6, E2_a3, k2_a3 = sortByKpt(C2_a5, p2_a5, k2_a2, E2_a2, nb,
                nk2x, nk2y, nKx, nKy, nKz)
        C3_a6, p3_a6, E3_a3, k3_a3 = sortByKpt(C3_a5, p3_a5, k3_a2, E3_a2, nb,
                nk3x, nk3y, nKx, nKy, nKz)
        return C2_a6, p2_a6, E2_a3, k2_a3, \
                C3_a6, p3_a6, E3_a3, k3_a3, \
                b_a2, ne, volume, area, encut

    else:
        return C2_a5, p2_a5, E2_a2, k2_a2, \
                C3_a5, p3_a5, E3_a2, k3_a2, \
                b_a2, ne, volume, area, encut


def readExport_kpts(
        infile = 'input',
        exroot = 'tmp/hBN.export',
        ):
    """Returns k-points used for QE run in fractional coordinates.

    Can read location of exroot from input file so that it can be run in same.
    diractory as newprob.txt file.
    """
    print(f'Reading from {infile}:')
    try:
        with open(f'{infile}') as f:

            for line in f:
                if 'exroot' in line:
                    exroot = line.split()[-1]
                    print(f'  exroot = {exroot}')

    except FileNotFoundError:
        print(f'  Could not find exroot location from {infile}.')
        pass

    print(f'  Looking for index.xml in {exroot}.')
    with open(f'{exroot}/index.xml') as f:
        for line in f:

            if '<Kpoints' in line:
                kline_list = line.split('"')
                nk = int(kline_list[1])
                nkx = int(kline_list[5])
                nky = int(kline_list[7])
                nb = int(f.readline().split('"')[1])

            if '<Cell' in line:
                alat = float(f.readline().split('"')[1])  # bohr

            if '<a3' in line:
                b_a2 = zeros((3, 3))
                for n in range(3):
                    b_a2[n] = array(f.readline().split('"')[1].split())\
                            .astype(float)  # bohr^{-1}
                invb_a2 = inv(b_a2)  # bohr

            if '<k' in line:
                blat = 2*pi/alat  # bohr^{-1}
                k_a2 = zeros((nk, 3))
                for k in range(nk):
                    k_ar = array([float(val) for val in
                            f.readline().split()])*blat  # bohr^{-1}
                    k_a2[k] = dot(k_ar, invb_a2)  # recip lat (fractional)

    # treat autogenerated k-point meshes differently from k-point lines
    if nkx > 0:

        # sort k-points in order of increasing crystal momenta (ky runs fastest)
        k_l2 = k_a2.tolist()
        for k, k_list in enumerate(k_l2):
            k_list.append(k)
        k_l2.sort()
        sort_list = [k_list[-1] for k_list in k_l2]
        sk_a2 = k_a2[sort_list]

        # sort k-points in order of increasing crystal momenta (ky runs fastest)
        return sk_a2.reshape((nkx, nky, 3))

    else:
        return k_a2

#--------------------------------- SCRATCH ------------------------------------

#                volume = float(line.split()[2].split('=')[1].strip('"')) \
#                        * bohrtoInvEV**3  # ev^{-3}
#            try:
#            except ValueError:
#                print("file = '{}', band = {}, nK = {}, K = ({}, {}, {})"\
#                        .format(wfc, n+1, len(pk_dict), Kx, Ky, Kz))
