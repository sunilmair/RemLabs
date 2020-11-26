from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase import Atoms
from ase.build import make_supercell
from ase.io import write
import numpy as np

# Si_alat = 5.427
Si_alat = 6.4012
structures_folder_path = 'cif_files/'

def make_Si_unitcell(write_file=False):
    """
    creates an Si unit cell with no Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    if write_file : write(structures_folder_path + 'unitcell_Si.cif', unitcell)
    structure = Struc(ase2struc(unitcell))
    structure.content['species']['Li'] = {'mass' : 6.939, 'kind' : 2}
    return structure


def make_Li_filled_unitcell(write_file=False):
    """
    creates a conventional Si unitcell of Si with all tetrahedral sites occupied by Li
    """
    Li_fractional_sites = [(0.5, 0.5, 0.5),
                           (0.25, 0.75, 0.25), (0.75, 0.25, 0.25), (0.25, 0.25, 0.75), (0.75, 0.75, 0.75),
                           (0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5)]

    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3*[Si_alat]+3*[90])
    unitcell.extend(Atoms('Li'+str(len(Li_fractional_sites)),
                          positions=[tuple([frac_coord*Si_alat for frac_coord in Li_frac_site])
                                     for Li_frac_site in Li_fractional_sites]))
    if write_file : write('Li_filled_unitcell.cif', unitcell)
    return Struc(ase2struc(unitcell))


def make_unitcell_central_Li(write_file=False):
    """
    creates an Si unit cell with the central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    unitcell.extend(Atoms('Li', positions=[tuple(3 * [0.5 * Si_alat])]))
    if write_file : write('unitcell_central_Li.cif', unitcell)
    return Struc(ase2struc(unitcell))


def make_unitcell_neighbor_central_Li(write_file=False):
    """
    creates an Si unit cell with a neighbor to the central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    unitcell.extend(Atoms('Li', positions=[tuple(3 * [0.75 * Si_alat])]))
    if write_file : write('unitcell_neighbor_central_Li.cif', unitcell)
    return Struc(ase2struc(unitcell))


def make_2x2x2_supercell_1x1x1_central_Li(write_file=False):
    """
    creates a 2x2x2 Si supercell with the 1, 1, 1 cell's central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    supercell = make_supercell(unitcell, np.identity(3) * 2)
    supercell.extend(Atoms('Li', positions=[tuple(3 * [0.5 * Si_alat])]))
    if write_file : write('2x2x2_supercell_1x1x1_central_Li.cif', supercell)
    return Struc(ase2struc(supercell))


def make_2x2x2_supercell_1x1x1_neighbor_central_Li(write_file=False):
    """
    creates a 2x2x2 Si supercell with the 1, 1, 1 cell's neighbor to the central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    supercell = make_supercell(unitcell, np.identity(3) * 2)
    supercell.extend(Atoms('Li', positions=[tuple(3 * [0.75 * Si_alat])]))
    if write_file : write('2x2x2_supercell_1x1x1_neighbor_central_Li.cif', supercell)
    return Struc(ase2struc(supercell))


def make_3x3x3_supercell_central_Li(write_file=False):
    """
    creates a 3x3x3 Si supercell with a single central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    supercell = make_supercell(unitcell, np.identity(3)*3)
    supercell.extend(Atoms('Li', positions=[tuple(3*[1.5*Si_alat])]))
    if write_file : write('3x3x3_supercell_central_Li.cif', supercell)
    return Struc(ase2struc(supercell))


def make_3x3x3_supercell_neighbor_central_Li(write_file=False):
    """
    creates a 3x3x3 Si supercell with a single site that neighbors the central tetrahedral site occupied by Li
    """
    unitcell = crystal('Si', [(0, 0, 0)], spacegroup=227, cellpar=3 * [Si_alat] + 3 * [90])
    supercell = make_supercell(unitcell, np.identity(3) * 3)
    supercell.extend(Atoms('Li', positions=[tuple(3 * [1.75 * Si_alat])]))
    if write_file : write('3x3x3_supercell_neighbor_central_Li.cif', supercell)
    return Struc(ase2struc(supercell))


if __name__ == "__main__":
    make_Si_unitcell(True)
    #make_Li_filled_unitcell(True)
    #make_unitcell_central_Li(True)
    #make_unitcell_neighbor_central_Li(True)
    #make_2x2x2_supercell_1x1x1_central_Li(True)
    #make_2x2x2_supercell_1x1x1_neighbor_central_Li(True)
    #make_3x3x3_supercell_central_Li(True)
    #make_3x3x3_supercell_neighbor_central_Li(True)
