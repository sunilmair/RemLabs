from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import matplotlib.pyplot as plt


def make_struc(size):
    """
    Creates the crystal structure using ASE.
    :param size: supercell multiplier
    :return: structure object converted from ase
    """
    alat = 4.10
    unitcell = crystal('Al', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    multiplier = numpy.identity(3) * size
    supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(supercell))
    return structure


def compute_dynamics(size, timestep, nsteps, temperature):
    """
    Make an input template and select potential and structure, and input parameters.
    Return a pair of output file and RDF file written to the runpath directory.
    """
    intemplate = """
    # ---------- Initialize simulation ---------------------
    units metal
    atom_style atomic
    dimension  3
    boundary   p p p
    read_data $DATAINPUT

    pair_style eam/alloy
    pair_coeff * * $POTENTIAL  Al

    velocity  all create $TEMPERATURE 87287 dist gaussian

    # ---------- Describe computed properties------------------
    compute msdall all msd
    thermo_style custom step pe ke etotal temp press density c_msdall[4]
    thermo $TOUTPUT

    # ---------- Specify ensemble  ---------------------
    fix  1 all nve
    #fix  1 all nvt temp $TEMPERATURE $TEMPERATURE $TDAMP

    # --------- Compute RDF ---------------
    compute rdfall all rdf 100 1 1
    fix 2 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfall[*] file $RDFFILE mode vector

    # --------- Run -------------
    timestep $TIMESTEP
    run $NSTEPS
    """

    potential = ClassicalPotential(ptype='eam', element='Al', name='Al_zhou.eam.alloy')
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab4/Problem1", "size_" + str(size)))
    struc = make_struc(size=size)
    inparam = {
        'TEMPERATURE': temperature,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,                 # how often to write thermo output
        'TDAMP': 50 * timestep,       # thermostat damping time scale
        'RDFFRAME': int(nsteps / 4),   # frames for radial distribution function
    }
    outfile = lammps_run(struc=struc, runpath=runpath, potential=potential,
                                  intemplate=intemplate, inparam=inparam)
    output = parse_lammps_thermo(outfile=outfile)
    rdffile = get_rdf(runpath=runpath)
    rdfs = parse_lammps_rdf(rdffile=rdffile)
    return output, rdfs


def md_run():
    output, rdfs = compute_dynamics(size=3, timestep=0.001, nsteps=1000, temperature=300)
    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    ## ------- plot output properties
    #plt.plot(simtime, temp)
    #plt.show()
    plt.plot(simtime, press)
    plt.show()

    # ----- plot radial distribution functions
    for rdf in rdfs:
        plt.plot(rdf[0], rdf[1])
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    md_run()
