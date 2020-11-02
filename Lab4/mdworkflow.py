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
    #fix  1 all nve
    fix  1 all nvt temp $TEMPERATURE $TEMPERATURE $TDAMP

    # --------- Compute RDF ---------------
    compute rdfall all rdf 100 1 1
    fix 2 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfall[*] file $RDFFILE mode vector

    # --------- Run -------------
    timestep $TIMESTEP
    run $NSTEPS
    """

    potential = ClassicalPotential(ptype='eam', element='Ag', name='Ag_u3.eam')
    #runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "RemLabs/Lab4/Problem1A", "timestep_" + str(timestep)))
    #runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "RemLabs/Lab4/Problem1C", "size_" + str(size)))
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "RemLabs/Lab4/Problem2A_size_" + str(size) , "temp_" + str(temperature)))
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


def md_analyze_timestep(total_time, timestep_smallest, timestep_largest, num_runs, size=3, temperature=300):
    timesteps = np.logspace(np.log10(timestep_smallest), np.log10(timestep_largest), num_runs)
    data = []
    mean_energies = []
    fig1, ax1 = plt.subplots(1, 2, figsize=(18, 6))
    for timestep in timesteps:
        savepath = '/home/modeler/RemLabs/Lab4/Problem1A/timestep_' + str(timestep) + '/'
        nsteps = int(np.ceil(total_time/timestep))

        output, rdfs = compute_dynamics(size, timestep, nsteps, temperature)
        output = output.astype(np.float)
        [simtime, pe, ke, energy, temp, pres, dens, msd] = output

        mean_energies.append(np.mean(energy[1:]))

        #single run plots here
        fig2, ax2 = plt.subplots(1, 2, figsize=(18, 6))
        ax2[0].plot(simtime*timestep, temp)
        ax2[0].set_title('Temp vs Time')
        ax2[0].set_xlabel('Time (ps)')
        ax2[0].set_ylabel('Temp (K)')

        ax2[1].plot(simtime*timestep, pe, color='tab:red')
        ax2[1].set_title('Energy vs Time')
        ax2[1].set_xlabel('Time (ps)')
        ax2[1].set_ylabel('PE (units)', color='tab:red') #change energy units
        ax2[1].tick_params(axis='y', labelcolor='tab:red')

        ax3 = ax2[1].twinx()
        ax3.plot(simtime*timestep, ke, color='tab:blue')
        ax3.set_ylabel('KE (units)', color='tab:blue') #change energy units
        ax3.tick_params(axis='y', labelcolor='tab:blue')

        fig2.tight_layout()
        fig2.savefig(savepath + 'temp.png')

        data.append(output)

        ax1[0].plot(simtime*timestep, energy, label=str(timestep)[:7])

    ax1[0].legend()
    ax1[0].set_title('Total Energy vs Time')
    ax1[0].set_xlabel('Time (ps)')
    ax1[0].set_ylabel('Total Energy (units)') #change energy units

    ax1[1].plot(timesteps, mean_energies, marker='o')
    ax1[1].set_title('Mean Energy vs Timestep')
    ax1[1].set_xlabel('Timestep (ps)')
    ax1[1].set_ylabel('Mean Energy (units)') #change energy units

    fig1.savefig('/home/modeler/RemLabs/Lab4/Problem1A/energies.png')


def md_analyze_supercell_size(sizes, timestep=0.001, nsteps=10000, temperature=300):
    fig1, ax1 = plt.subplots(1, 2, figsize=(12, 6))
    temp_stds = []
    for size in sizes:
        savepath = '/home/modeler/RemLabs/Lab4/Problem1C/size_' + str(size) + '/'

        output, rdfs = compute_dynamics(size, timestep, nsteps, temperature)
        output = output.astype(np.float)

        [simtime, pe, ke, energy, temp, pres, dens, msd] = output
        temp_stds.append(np.std(temp[1:]))

        ax1[0].plot(simtime*timestep, temp, label=str(size))

    ax1[0].legend()
    ax1[0].set_title('Temp vs Time')
    ax1[0].set_xlabel('Time (ps)')
    ax1[0].set_ylabel('Temp (K)')

    ax1[1].plot(sizes, temp_stds, marker='o')
    ax1[1].set_title('Temp Std Dev vs Supercell Size')
    ax1[1].set_xlabel('Supercell Size')
    ax1[1].set_ylabel('Temp Std Dev (K)')

    fig1.savefig('/home/modeler/RemLabs/Lab4/Problem1C/temps_sizes.png')


def md_analyze_melt_nvt(min_temp, max_temp, temp_step, size, timestep, nsteps):
    temperatures = np.arange(min_temp, max_temp, temp_step)
    means = []
    fig1, ax1 = plt.subplots(1, 2, figsize=(12, 6))
    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6))
    for temperature in temperatures:
        savepath = '/home/modeler/RemLabs/Lab4/Problem2A_size_' + str(size) + '/temp_' + str(temperature) + '/'

        output, rdfs = compute_dynamics(size, timestep, nsteps, temperature)
        output = output.astype(np.float)

        [simtime, pe, ke, energy, temp, pres, dens, msd] = output
        mean = []
        for col in output:
            mean.append(np.mean(col[1:]))
        means.append(mean)

        ax1[0].plot(simtime*timestep, msd, label=str(temperature))
        ax1[1].plot(simtime*timestep, energy, label=str(temperature))

    ax1[0].legend()
    ax1[0].set_title('MSD  vs Time')
    ax1[0].set_xlabel('Time (ps)')
    ax1[0].set_ylabel('MSD (units)') # change units

    ax1[1].legend()
    ax1[1].set_title('Energy vs Time')
    ax1[1].set_xlabel('Time (ps)')
    ax1[1].set_ylabel('Energy (units)') # change units

    fig1.savefig('/home/modeler/RemLabs/Lab4/Problem2A_size_' + str(size) + '/time.png')

    ax2.plot(temperatures, [entry[7] for entry in means], color='tab:red')
    ax2.set_title('MSD and Energy vs Temp')
    ax2.set_xlabel('Temp (K)')
    ax2.set_ylabel('MSD (units)', color='tab:red') # change  units
    ax2.tick_params(axis='y', labelcolor='tab:red')

    ax3 = ax2.twinx()
    ax3.plot(temperatures, [entry[3] for entry in means], color='tab:blue')
    ax3.set_ylabel('Energy (units)', color='tab:blue') # change units
    ax3.tick_params(axis='y', labelcolor='tab:blue')

    fig2.tight_layout()
    fig2.savefig('/home/modeler/RemLabs/Lab4/Problem2A_size_' + str(size) + '/temp.png')


if __name__ == '__main__':
    # put here the function that you actually want to run
    #md_analyze_timestep(10, 0.001, 0.02, 8) # 1A and 1B
    #md_analyze_supercell_size([3, 4, 5]) # 1C
    md_analyze_melt_nvt(200, 4000, 200, 7, 0.001, 10000)
