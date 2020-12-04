from LiSi_structures import *
from LiSi_LAMMPS_templates import *
import time

project_full_path = '/home/modeler/RemLabs/Project'
def Si_lattice_constant_calculation():
    """
    Runs an optimization calculation on an Si unitcell to find the equilibrium lattice parameter using the 2NN MEAM
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, 'Si_lattice_constant', timestamp)

    runpath = Dir(path=path)
    struc = make_Si_unitcell()
    inparam = {'OUTFILE': path}

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=relaxation_calculation_template, inparam=inparam)
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


def Si_unitcell_central_Li_relaxation():
    """
    Runs an optimization calculation on an Si unitcell with an Li in the central tetrahedral site
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, 'Si_unitcell_central_Li_relaxation', timestamp)

    runpath = Dir(path=path)
    struc = make_unitcell_central_Li()
    inparam = {'OUTFILE': path}

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=relaxation_calculation_template, inparam=inparam)
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


def Si_3x3x3_supercell_central_Li_relaxation():
    """
    Runs an optimization calculation on an Si unitcell with an Li in the central tetrahedral site
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, 'Si_3x3x3_supercell_central_Li_relaxation', timestamp)

    runpath = Dir(path=path)
    struc = make_3x3x3_supercell_central_Li()
    inparam = {'OUTFILE': path}

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=relaxation_calculation_template, inparam=inparam)
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice

if __name__ == "__main__":
    #Si_lattice_constant_calculation()
    #Si_unitcell_central_Li_relaxation()
    Si_3x3x3_supercell_central_Li_relaxation()
