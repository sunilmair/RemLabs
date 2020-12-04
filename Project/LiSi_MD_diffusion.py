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
    filepath = 'eval_timestep_recheckdelete'
    timestep_list = np.logspace(np.log10(0.005), np.log10(0.05), 20)
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


def Si_n3_supercell_run_equil_MD(n, T, timestep, equilnsteps, nsteps, filepath):
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
        'EQUILNSTEPS': equilnsteps,
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


def test_equil_run(n, T, timestep, equilnsteps, nsteps, filepath):
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 6))
    output = Si_n3_supercell_run_equil_MD(n, T, timestep, equilnsteps, nsteps, filepath)
    output = output.astype(np.float)
    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    ax_left.plot([timestep*simtimestep for simtimestep in simtime], energy)
    ax_right.plot([timestep*simtimestep for simtimestep in simtime], msd)
    fig.savefig(filepath+'/test-'+ time.strftime('%Y%m%d-%H%M%S'))


if __name__ == "__main__":
    evaluate_timestep()
    #test_equil_run(3, 1800, 0.003, 3200, 5000, 'test_equil')
