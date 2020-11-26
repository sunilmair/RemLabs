from LiS_structures import *
from LiS_LAMMPS_templates import *
import time

def Si_lattice_constant_calculation():
    """
    Runs a optimization calculation on an Si unitcell to find the equilibrium lattice parameter using the 2NN MEAM
    """
    runpath = Dir(path=os.path.join('Si_lattice_constant', time.strftime(%Y%m%d-%H%M%S)))
    struc = make_Si_unitcell()
    output_file = lammps_run(struc=struc, runpath=runpath, intemplate=Si_lattice_constant_template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


if __name__ == "__main__":
    Si_lattice_constant_calculation()