from numpy import pi, ones
from numpy import array, zeros, dot
from numpy.linalg import norm, inv

import os

from constants import *

"""Get coefficients C^n_{G+k} and PW momenta from pw_export output

to do:
* accept arbitrary k-point dimensions
* automatically detect <dir>/<prefix>.export
"""

def readExport(
        exroot = 'tmp/hBN.export',
        ):
    """Returns arrays containing momenta and PW coefficients.

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

    exroot: directory containing pw_export outputs (str)
        * usually {outdir}/{prefix}.export
    """
    with open(f'{exroot}/index.xml') as f:
        for line in f:

            if '<Kpoints' in line and 'nk' not in locals():
                kline_list = line.split('"')
                nk = int(kline_list[1])
                nkx = int(kline_list[5])
                nky = int(kline_list[7])
                nb = int(f.readline().split('"')[1])

            if '<Cell' in line and 'volume' not in locals():
                a_eq, v_eq = f.readline().split()[1:3]
                alat = float(a_eq.split('"')[1])*bohrtoInvEV  # eV^{-1} 
                volume = float(v_eq.split('"')[1])*bohrtoInvEV**3  # eV^{-3}

                a_a2 = zeros((3, 3))
                for n in range(3):
                    a_a2[n] = array(f.readline().split('"')[1].split())\
                            .astype(float)  # bohr
                a_a2 *= bohrtoInvEV  # eV^{-1}
                area = volume / norm(a_a2[2])  # eV^{-2}

                b_a2 = zeros((3, 3))
                for n in range(3):
                    b_a2[n] = array(f.readline().split('"')[1].split())\
                            .astype(float)  # bohr^{-1}
                b_a2 *= invBohrtoEV  # eV

            if '<k' in line and 'k_a2' not in locals():
                blat = 2*pi/alat
                k_a2 = zeros((nk, 3))
                for k in range(nk):
                    k_a2[k] = array([float(val) for val in
                            f.readline().split()])*blat  # eV

            if '<Cutoff' in line and 'encut' not in locals():
                encut = float(line.split('"')[1]) * rydtoEV  # eV

            if '<Charge' in line and 'ne' not in locals():
                ne = int(round(float(line.split('"')[1])))

            if '<Eigenvalues' in line and 'E_a2' not in locals():
                E_a2 = zeros((nk, nb))
                for k in range(nk):
                    f.readline()
                    for b in range(nb):
                        E_a2[k, b] = float(f.readline())*rydtoEV  # eV
                    f.readline()

    # find max Kx, Ky, and Kz so that all levels have same dimension
    maxKx, maxKy, maxKz = 0, 0, 0
    minKx, minKy, minKz = 0, 0, 0

    # record all RLVs listed in all grid files
    grid_list = ['{}/{}'.format(exroot, d) for d in os.listdir(exroot)
            if 'grid.' in d]
    grid_list.sort(key = fileSort)
    p_list = []  # p_dict[k] --> pk_dict

    for k_ar, grid in zip(k_a2, grid_list):
        pk_dict = {}  # pk_dict[(Kx, Ky, Kz)] --> p_ar

        with open(grid) as f:
            for _ in range(5):
                f.readline()
            nK = int(f.readline().split()[1].split('=')[1].strip('"'))
            for _ in range(nK + 3):
                f.readline()

            # store RLVs
            K_list = []
            for n in range(nK):
                K_ar = [int(val) for val in f.readline().split()]
                Kx, Ky, Kz = K_ar

                if Kx > maxKx:
                    maxKx = Kx
                elif Kx < minKx:
                    minKx = Kx

                if Ky > maxKy:
                    maxKy = Ky
                elif Ky < minKy:
                    minKy = Ky

                if Kz > maxKz:
                    maxKz = Kz
                elif Kz < minKz:
                    minKz = Kz

                p_ar = dot(K_ar, b_a2) + k_ar  # eV
                pk_dict[(Kx, Ky, Kz)] = p_ar

        p_list.append(pk_dict)

    # array ranges
    nKx = maxKx - minKx + 1 
    nKy = maxKy - minKy + 1 
    nKz = maxKz - minKz + 1 

    # transfer momenta to an array
    p_a5 = ones((nk, nKx, nKy, nKz, 3))
    for k, pk_dict in enumerate(p_list):
        for (Kx, Ky, Kz), p_ar in pk_dict.items():
            p_a5[k, Kx, Ky, Kz] = p_ar

    # record all coefficients listed all wfc files
    wfc_list = ['{}/{}'.format(exroot, d) for d in os.listdir(exroot)
            if d[:4] == 'wfc.']
    wfc_list.sort(key = fileSort)
    C_a5 = zeros((nk, nb, nKx, nKy, nKz), dtype=complex)

    for k, wfc, pk_dict in zip(range(nk), wfc_list, p_list):
        with open(wfc) as f:
            for line in f:
                if 'igwx=' in line:
    
                    # store planewaves coefficients
                    for n in range(nb):
                        f.readline()
                        for Kx, Ky, Kz in pk_dict:
                            reC, ImC = [float(val) for val in
                                    f.readline().split(',')]
                            C_a5[k, n, Kx, Ky, Kz] = reC + 1j*ImC
                        f.readline()

    # treat autogenerated k-point meshes differently from k-point lines
    if nkx > 0:

        # sort k-points in order of increasing crystal momenta (ky runs fastest)
        k_l2 = k_a2.tolist()
        for k, k_list in enumerate(k_l2):
            k_list.append(k)
        k_l2.sort()
        sort_list = [k_list[-1] for k_list in k_l2]

        # sort other values to match k-point ordering
        sk_a2 = k_a2[sort_list]
        sE_a2 = E_a2[sort_list]
        sp_a5 = p_a5[sort_list]
        sC_a5 = C_a5[sort_list]

        k_a3 = sk_a2.reshape((nkx, nky, 3))
        E_a3 = sE_a2.reshape((nkx, nky, nb))
        p_a6 = sp_a5.reshape((nkx, nky, nKx, nKy, nKz, 3))
        C_a6 = sC_a5.reshape((nkx, nky, nb, nKx, nKy, nKz))

        return C_a6, p_a6, E_a3, k_a3, b_a2, ne, volume, area, encut

    else:
        return C_a5, p_a5, E_a2, k_a2, b_a2, ne, volume, area, encut


def fileSort(s):
    """Key to sort grid.{k} and wfc.{k} in order of k-points listed in k_a2."""
    _, n = s.split('/')[-1].split('.')
    return int(n)


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
