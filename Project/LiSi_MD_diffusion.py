from LiSi_structures import *
from LiSi_LAMMPS_templates import *
import time
import matplotlib.pyplot as plt

project_full_path = '/home/modeler/RemLabs/Project'

def Si_n3_supercell_run_MD(n, T, timestep, nsteps, filepath):
    """
    Runs an npt ensemble to track the MSD of a Li in a 3x3x3 Si supercell to calculate the diffusion coefficient
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, filepath,
                        'n'+str(n)+'_T'+str(T)+'_timestep'+'{:2e}'.format(timestep)+'_steps'+str(nsteps), timestamp)

    runpath = Dir(path=path)
    struc = make_n3_supercell_1x1x1_central_Li(n)
    inparam = {
        'OUTFILE': path,
        'TEMPERATURE': T,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,  # how often to write thermo output
        'TDAMP': 50 * timestep,  # thermostat damping time scale
    }

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=MD_npt,
                             inparam=inparam)
    output = parse_lammps_thermo(outfile=output_file)
    output = output.astype(np.float)

    [simtime, pe, ke, energy, temp, press, dens, msd] = output

    return output


def evaluate_timestep():
    """
    Use mean energy (after equilibration) as a convergence metric for timestep size
    """
    filepath = 'eval_timestep'
    timestep_list = np.logspace(np.log10(0.0005), np.log10(0.005), 20)
    total_time = 15
    equilibration_time = 2
    n = 3
    T = 1800

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 6))
    mean_energy_list = []

    for timestep in timestep_list:
        output = Si_n3_supercell_run_MD(n, T, timestep, int(np.ceil(total_time/timestep)), filepath)
        [simtime, pe, ke, energy, temp, press, dens, msd] = output
        ax_left.plot([timestep*simtimestep for simtimestep in simtime], energy, label='{:.3}'.format(timestep))
        mean_energy_list.append(np.mean(
            [energy[i] for i in range(len(simtime)) if simtime[i] > equilibration_time]))

    ax_right.plot(timestep_list, mean_energy_list, marker='o')
    ax_right_min, ax_right_max = ax_right.get_ylim()
    ax_right_twin = ax_right.twinx()
    ax_right_twin.set_ylim(ax_right_min - mean_energy_list[0], ax_right_max - mean_energy_list[0])

    ax_left.set_xlabel('Time (ps)')
    ax_left.set_ylabel('Total Energy') # units?
    ax_left.set_title('Total Energy vs Time')
    ax_left.legend(loc='center right')

    ax_right.set_xlabel('Timestep (ps)')
    ax_right.set_ylabel('Mean Energy') # units?
    ax_right_twin.set_ylabel('Energy Convergence')
    ax_right.set_title('Mean Energy vs Timestep')

    fig.savefig(filepath+'/timestep_convergence-' + time.strftime('%Y%m%d-%H%M%S'))


def Si_n3_supercell_run_equil_MD(n, T, timestep, equilnsteps, production_time, filepath, intemplate):
    """
    Runs an npt ensemble to track the MSD of a Li in a 3x3x3 Si supercell to calculate the diffusion coefficient
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, filepath,
                        'n'+str(n)+'_T'+str(T)+'_timestep'+'{:2e}'.format(timestep)+
                        '_eqsteps'+str(equilnsteps)+'_time'+str(production_time), timestamp)

    runpath = Dir(path=path)
    struc = make_n3_supercell_1x1x1_central_Li(n)
    nsteps = int(np.ceil(production_time/timestep))
    inparam = {
        'OUTFILE': path,
        'TEMPERATURE': T,
        'EQUILNSTEPS': equilnsteps,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,  # how often to write thermo output
        'TDAMP': 50 * timestep,  # thermostat damping time scale
    }

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=intemplate,
                             inparam=inparam)
    output = parse_lammps_thermo(outfile=output_file)
    output = output[2:]
    output = [element for element in output]
    output = np.array(output)
    output = output.astype(np.float)
    outrows = np.transpose(np.array(output))

    [simtime, pe, ke, energy, temp, press, dens, msdli, msdsi] = outrows

    return outrows


def test_equil_run(n, T, timestep, equilnsteps, production_time, filepath, intemplate):
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 6))
    output = Si_n3_supercell_run_equil_MD(n, T, timestep, equilnsteps, production_time, filepath, intemplate)
    [simtime, pe, ke, energy, temp, press, dens, msdli, msdsi] = output
    ax_left.plot([timestep*simtimestep for simtimestep in simtime], energy)
    ax_right.plot([timestep*simtimestep for simtimestep in simtime], msdli, label='Li')
    ax_right.plot([timestep*simtimestep for simtimestep in simtime], msdsi, label='Si')
    ax_right.legend()

    ax_left.set_xlabel('Time(ps)')
    ax_left.set_ylabel('Energy')
    ax_left.set_title('Energy vs Time')

    ax_right.set_xlabel('Time(ps)')
    ax_right.set_ylabel('MSD (A)')
    ax_right.set_title('MSD vs Time')

    fig.savefig(filepath+'/test-'+ time.strftime('%Y%m%d-%H%M%S'))


def Si_n3_supercell_run_equil_MD_rseed(n, T, timestep, equilnsteps, production_time, filepath, intemplate, rseed):
    """
    Runs an npt ensemble to track the MSD of a Li in a 3x3x3 Si supercell to calculate the diffusion coefficient
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, filepath,
                        'n'+str(n)+'_T'+str(T)+'_timestep'+'{:2e}'.format(timestep)+
                        '_eqsteps'+str(equilnsteps)+'_time'+str(production_time)+'_rseed'+str(rseed), timestamp)

    runpath = Dir(path=path)
    struc = make_n3_supercell_1x1x1_central_Li(n)
    nsteps = int(np.ceil(production_time/timestep))
    inparam = {
        'OUTFILE': path,
        'TEMPERATURE': T,
        'EQUILNSTEPS': equilnsteps,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,  # how often to write thermo output
        'TDAMP': 50 * timestep,  # thermostat damping time scale
        'RSEED': rseed,
    }

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=intemplate,
                             inparam=inparam)
    output = parse_lammps_thermo(outfile=output_file)
    output = output[2:]
    output = [element for element in output]
    output = np.array(output)
    output = output.astype(np.float)
    outrows = np.transpose(np.array(output))

    return outrows


def get_MSD(n, T, production_time, num_runs, filepath):
    timestep = 0.003
    equilnsteps = 3200

    msdli_list = []

    for i in range(num_runs):
        print(i)

        output = Si_n3_supercell_run_equil_MD_rseed(n, T, timestep, equilnsteps, production_time, filepath, MD_equilibrate_npt_track_MSD_rseed, int(3320 + i))
        [simtime, pe, ke, energy, temp, press, dens, msdli, msdsi] = output
        msdli_list.append(msdli)

    return simtime, msdli_list

def calc(n, Tstart, Tstop, numT, production_time, num_runs, filepath):
    T_list = np.linspace(Tstart, Tstop, numT)
    msdli_list_list = []
    D_list = []
    simtime_list = []

    for T in T_list:
        print(T)

        simtime, msdli_list = get_MSD(n, T, production_time, num_runs, filepath+'/T'+str(T))
        msdli_list_list.append(msdli_list)
        simtime_list.append(simtime)

    fig, axs = plt.subplots(numT, 1, figsize=(18, 6))
    for i in range(numT):

        for j in range(len(msdli_list_list[i])):
            axs[i].plot(simtime_list[i], msdli_list_list[i][j], linewidth=0.75, alpha=0.5)

        avg_msdli = []
        for j in range(len(simtime_list[i])):
            avg_msdli.append(np.mean([list[j] for list in msdli_list_list[i]]))
        axs[i].plot(simtime_list[i], avg_msdli, linewidth=2)

        np_simtime = np.asarray(simtime_list[i])
        np_avg_msdli = np.asarray(avg_msdli)
        num = np.sum([np_simtime[i]*np_avg_msdli[i] for i in range(len(simtime[i]))])
        den = np.sum([element**2 for element in np_simtime])
        slope = num/den

        axs[i].plot(simtime_list[i], [slope*simtime_element for simtime_element in simtime_list[i]], linewidth=2)
        #get D, abandon histogram idea
    fig.savefig(filepath+'/T_series')




if __name__ == "__main__":
    #evaluate_timestep()

    #test_equil_run(3, 1600, 0.003, 3200, 600, 'test_equil_nvt_nn_day2', MD_equilibrate_nvt_track_MSD)
    #test_equil_run(3, 1600, 0.003, 3200, 600, 'test_equil_npt_nn_day2', MD_equilibrate_npt_track_MSD)
    #test_equil_run(3, 1600, 0.003, 3200, 600, 'test_equil_npt_then_nvt_nn_day2', MD_equilibrate_npt_track_MSD_nvt)

    calc(3, 1600, 1800, 2, 15, 2, 'arrhenius')
