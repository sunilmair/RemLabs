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


def make_na2_nh2_bh4_unitcell_ase(write_file=False):

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


def make_na2_nh2_bh4_nxnxn_ase(n, write_file=False):
    nnb_unitcell = make_na2_nh2_bh4_unitcell_ase()
    nnb_nxnxn = make_supercell(nnb_unitcell, np.identity(3)*n)
    if write_file: write(structures_folder_path + 'nnb_'+str(n)+'x'+str(n)+'x'+str(n)+'.cif', nnb_nxnxn)
    return nnb_nxnxn

def make_initial_unitcell_state():
    return Struc(ase2struc(make_na2_nh2_bh4_unitcell_ase()))


def make_initial_2x2x2_neb_state():
    return Struc(ase2struc(make_na2_nh2_bh4_nxnxn_ase(2)))


def make_final_2x2x2_neb_state():
    final_state = make_na2_nh2_bh4_nxnxn_ase(2)
    final_state.positions[0] = [nnb_lattice_parameter*coord for coord in [0.5, 0.5, 0]]
    return Struc(ase2struc(final_state))


def make_nnb_relaxed_structure_ase(write_file=False):
    name = 'Na2NH2BH4'
    Na1_pos = (2.287693187, 0.118479557, 2.223873283)
    Na2_pos = (0.126266155, 2.302452064, 2.610427732)

    N_pos = (2.177422414, 2.340089543, 1.609427022)
    NH1_pos = (2.791658178, 1.894258226, 0.909980651)
    NH2_pos = (1.586837022, 2.855421231, 0.937307902)

    B_pos = (-0.002127017, 0.050662897, 0.244878644)
    BH1_pos = (1.108110027, 0.077054367,  -0.236510738)
    BH2_pos = (-0.664951296, 0.990285950, -0.139625958)
    BH3_pos = (-0.609179425, -0.964593720, -0.007842135)
    BH4_pos = (0.096958556, 0.132031215, 1.452070786)

    positions = [Na1_pos, Na2_pos,
                 N_pos, NH1_pos, NH2_pos,
                 B_pos, BH1_pos, BH2_pos, BH3_pos, BH4_pos]

    cell = [[4.325199939, 0.064686593, -0.127337173],
            [0.071390375, 4.429954029, 0.028242632],
            [-0.140528465, 0.019613305, 4.154985085]]

    nnb_relaxed_unitcell = Atoms(name,
                                 positions=positions,
                                 cell=cell,
                                 pbc=[1,1,1])

    if write_file: write(structures_folder_path + 'nnb_relaxed_unitcell.cif', nnb_relaxed_unitcell)

    return nnb_relaxed_unitcell

def make_nnb_relaxed_structure_2_ase(write_file=False):
    name = 'Na2NH2BH4'
    Na1_pos = (2.126845646, 0.027255500, 2.021437822)
    Na2_pos = (-0.030716908, 2.199824696, 2.075527861)

    N_pos = (2.110284350, 2.235160446, 1.278166977)
    NH1_pos = (2.657228411, 1.691358084, 0.592610032)
    NH2_pos = (1.588207029, 2.792272903, 0.583835092)

    B_pos = (0.000000000, 0.000000000, 0.000000000)
    BH1_pos = (1.161388732, 0.031772517, -0.325051421)
    BH2_pos = (-0.584702037, 0.976547059, -0.421391252)
    BH3_pos = (-0.576965713, -0.999650076, -0.373722808)
    BH4_pos = (-0.057832071, 0.020226134, 1.214530881)

    positions = [Na1_pos, Na2_pos,
                 N_pos, NH1_pos, NH2_pos,
                 B_pos, BH1_pos, BH2_pos, BH3_pos, BH4_pos]

    cell = [[4.381420580, 0.055575440, -0.085954569],
            [0.066295472, 4.400713833 , 0.022225508],
            [-0.127729288, -0.001064441, 4.128830213]]

    nnb_relaxed_unitcell = Atoms(name,
                                 positions=positions,
                                 cell=cell,
                                 pbc=[1,1,1])

    if write_file: write(structures_folder_path + 'nnb_relaxed_unitcell.cif', nnb_relaxed_unitcell)

    return nnb_relaxed_unitcell


if __name__ == '__main__':
    make_nnb_relaxed_structure_2_ase(True)
