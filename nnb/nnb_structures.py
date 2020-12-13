from ase import Atoms
from ase.build import make_supercell
from ase.io import write
import numpy as np
from labutil.src.objects import *

nnb_lattice_parameter = 4.7

BH4_bond_length = 1.17
NH2_bond_length = 0.98
structures_folder_path = 'cif_files/'

def get_tetrahedron_coords(bond_len):
    """
    returns vertices of tetrahedron centered at (0, 0, 0) with given 'bond length'
    returns list of lists [[x1, y1, z1], [x2, y2, z2], ...]
    """
    unit_tetrahedron = [[np.sqrt(8/9), 0, -1/3],
     [-np.sqrt(2/9), np.sqrt(2/3), -1/3],
     [-np.sqrt(2/9), -np.sqrt(2/3), -1/3],
     [0, 0, 1]]
    return [[bond_len*coord for coord in vertex] for vertex in unit_tetrahedron]


def make_na2_nh2_bh4_unitcell(write_file=False):

    name = 'Na2NH2BH4'
    BH4_tetrahedron = get_tetrahedron_coords(BH4_bond_length/nnb_lattice_parameter)
    NH2_tetrahedron = get_tetrahedron_coords(NH2_bond_length/nnb_lattice_parameter)

    Na1_pos = (0.5, 0, 0.5)
    Na2_pos = (0, 0.5, 0.5)

    N_pos = (0.5, 0.5, 0.5)
    NH1_pos = tuple(N_pos[i] + NH2_tetrahedron[0][i] for i in range(3))
    NH2_pos = tuple(N_pos[i] + NH2_tetrahedron[1][i] for i in range(3))

    B_pos = (0, 0, 0)
    BH1_pos = tuple(B_pos[i] + BH4_tetrahedron[0][i] for i in range(3))
    BH2_pos = tuple(B_pos[i] + BH4_tetrahedron[1][i] for i in range(3))
    BH3_pos = tuple(B_pos[i] + BH4_tetrahedron[2][i] for i in range(3))
    BH4_pos = tuple(B_pos[i] + BH4_tetrahedron[3][i] for i in range(3))


    fractional_positions = [Na1_pos, Na2_pos,
                            N_pos, NH1_pos, NH2_pos,
                            B_pos, BH1_pos, BH2_pos, BH3_pos, BH4_pos]

    absolute_positions = [tuple(nnb_lattice_parameter*coord for coord in fractional_position)
                          for fractional_position in fractional_positions]

    nnb_unitcell = Atoms(name,
                         positions=absolute_positions,
                         cell=3*[nnb_lattice_parameter]+3*[90],
                         pbc=[1, 1, 1])

    if write_file: write(structures_folder_path + 'nnb_unitcell.cif', nnb_unitcell)

    return nnb_unitcell


def make_na2_nh2_bh4_nxnxn(n, write_file=False):
    nnb_unitcell = make_na2_nh2_bh4_unitcell()
    nnb_nxnxn = make_supercell(nnb_unitcell, np.identity(3)*n)
    if write_file: write(structures_folder_path + 'nnb_'+str(n)+'x'+str(n)+'x'+str(n)+'.cif', nnb_nxnxn)
    return nnb_nxnxn


def make_initial_state():
    return Struc(ase2struc(make_na2_nh2_bh4_nxnxn(2)))


def make_final_state():
    final_state = make_na2_nh2_bh4_nxnxn(2)
    final_state.positions[0] = [nnb_lattice_parameter*coord for coord in [0.5, 0.5, 0]]
    return Struc(ase2struc(final_state))


if __name__ == '__main__':
    make_na2_nh2_bh4_nxnxn(2, True)