from numpy import array, zeros, cross, dot, pi, arange
from numpy import newaxis as na
from numpy.linalg import norm, inv
from os import makedirs

from constants import hbar, c

"""Prepares various types of quantumEspresso calculations."""

invAtoEV = hbar*c*1e10  # inverse angstrom to eV

def prep_kcon(qmax=2000, qmin=300, infile='qscf.in'):
    """Prepares input files for k-point mesh convergence.

    mesh increases by integer value of
    q{max,min}: kpt spacing in eV for {coars,fin}est mesh (pos float)
    infile: qe input file with automatic k-point mesh (str)
    """
    # get ratios of RLV magnitudes
    with open(infile) as f:
        for line in f:
            if 'CELL_PARAMETERS' in line:
                a_l2 = []
                for n in range(3):
                    a_l2.append([float(x) for x in f.readline().split()])
                b_a2 = get_recip(array(a_l2))  # RLV (A^{-1})
                b_ar = norm(b_a2, axis=1)*invAtoEV  # RLV magnitudes (eV)
                max_idx = b_ar.argmax()
                bmax = b_ar[max_idx]  # largest RLV magnitude (eV)
                rat_ar = b_ar/bmax  # ratios against largest RLV
                print(rat_ar)
                break

    # get smallest and largets kpt mesh dims along largest RLV
    mindim = round(bmax/qmax)
    maxdim = round(bmax/qmin)

    # prepare kpts for every mesh from smallest to largest
    dim_ar = arange(mindim, maxdim+1)
    k_a2 = (dim_ar[:,na]*rat_ar[na,:]).round().astype(int)

    return k_a2

    # find line that defines k-point mesh
    line_list = read_input(infile)
    for n, line in enumerate(line_list):
        if 'K_POINTS' in line:
            kline_idx = n+1

    # write new input file for each mesh
    for k_ar in k_a2:
        k_line = '{} {} {} 0 0 0\n'.format(*k_ar)
        outdir = '{:0>2d}x{:0>2d}x{:0>2d}'.format(*k_ar)
        makedirs(outdir, exist_ok=True)
        outpath = f'{outdir}/{infile}'
        with open(outpath, 'w') as f:
            for n, line in enumerate(line_list):
                if n==kline_idx:
                    f.write(k_line)
                else:
                    f.write(line)


def mod_kpts(k_list = [3, 3, 1], infile='qscf.in'):
    """Prepares input file with new k-point mesh dimensions.

    nkx: number of k-points along first RLV (pos int)
    infile: qe input file with automatic k-point mesh (str)
    """
    line_list = read_input(infile)

    newline_list = []
#     for line in line_list:
    with open(infile) as f:
        for line in f:

            if 'K_POINTS (automatic)' in line:
                newkpts_list = []
    
    return b_a2

def get_recip(a_a2):
    """returns reciprical cell as a rank 2 array."""
    b_a2 = zeros([3, 3])
    
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        cross_product = cross(a_a2[j], a_a2[k])
        volume = dot(a_a2[i], cross_product)
        b_a2[i] = 2 * pi * cross_product / volume

    return b_a2

def read_input(infile='qscf.in'):
    """returns input files as list of lines."""
    with open(infile) as f:
        line_list = f.readlines()
    return line_list
        
#----------------------------------SCRATCH ------------------------------------

#     maxdim_fl = bmax/qmax  # largest k-point mesh dim (RLV frac)
#     maxdim = round(maxdim_fl)  # make largest mesh dim an integer
#     rat = maxdim/maxdim_fl
#     k_ar = (b_ar/qmax*rat).round().astype(int)
