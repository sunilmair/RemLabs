from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import matplotlib.pyplot as plt
import numpy as np


def make_struc(size):
    """
    Creates the crystal structure using ASE.
    :param size: supercell multiplier
    :return: structure object converted from ase
    """
    alat = 4.09
    unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
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

    pair_style eam
    pair_coeff * * $POTENTIAL

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

    potential = ClassicalPotential(ptype='eam', element='Ag', name='Ag_u3.eam')
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "RemLabs/Lab4/Problem1A", "timestep_" + str(timestep)))
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


def md_run(size=3, timestep=0.001):
    savepath = '/home/modeler/RemLabs/Lab4/Problem1/timestep_' + str(timestep) + '/'
    output, rdfs = compute_dynamics(size, timestep, nsteps=1000, temperature=300)
    print(output)
    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    ## ------- plot output properties
    fig1, ax1 = plt.subplots()
    ax1.plot(simtime, temp)


    plt.plot(simtime, temp)
    #plt.show()
    plt.savefig(savepath + 'temp.png')
    plt.plot(simtime, press)
    #plt.show()
    plt.savefig(savepath + 'press.png')

    # ----- plot radial distribution functions
    for rdf in rdfs:
        plt.plot(rdf[0], rdf[1])
    #plt.show()
    plt.savefig(savepath + 'rdf.png')


def md_analyze_timestep(total_time, timestep_smallest, timestep_largest, num_runs, size=3, temperature=300):
    timesteps = np.logspace(np.log10(timestep_smallest), np.log10(timestep_largest), num_runs)
    data = []
    mean_energies = []
    fig1, ax1 = plt.subplots()
    for timestep in timesteps:
        savepath = '/home/modeler/RemLabs/Lab4/Problem1A/timestep_' + str(timestep) + '/'
        nsteps = int(np.ceil(total_time/timestep))

        output, rdfs = compute_dynamics(size, timestep, nsteps, temperature)
        output = output.astype(np.float)
        [simtime, pe, ke, energy, temp, pres, dens, msd] = output

        mean_energies.append(np.mean(energy))

        #single run plots here
        fig2, ax2 = plt.subplots()
        ax2.plot(simtime*timestep, energy)
        fig2.savefig(savepath + 'energy.png')

        data.append(output)
        ax1.plot(simtime*timestep, energy, label=str(timestep)[:4])

    fig1.legend()

    fig3, ax3 = plt.subplots()
    ax3.plot(timesteps, mean_energies)

    fig1.savefig('/home/modeler/RemLabs/Lab4/Problem1A/energy_vs_time.png')
    fig3.savefig('/home/modeler/RemLabs/Lab4/Problem1A/mean_energies.png')

if __name__ == '__main__':
    # put here the function that you actually want to run
    md_analyze_timestep(10, 0.001, 0.02, 3)
