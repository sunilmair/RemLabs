from LiSi_structures import *
from LiSi_LAMMPS_templates import *
import time
import matplotlib.pyplot as plt

project_full_path = '/home/modeler/RemLabs/Project'

def Si_3x3x3_supercell_run_MD(temperature, timestep, nsteps):
    """
    Runs an npt ensemble to track the MSD of a Li in a 3x3x3 Si supercell to calculate the diffusion coefficient
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, 'Si_3x3x3_supercell_MD', str(temperature), timestamp)

    runpath = Dir(path=path)
    struc = make_3x3x3_supercell_central_Li()
    inparam = {
        'OUTFILE': path,
        'TEMPERATURE': temperature,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,  # how often to write thermo output
        'TDAMP': 50 * timestep,  # thermostat damping time scale
    }

    output_file = lammps_run(struc=struc, runpath=runpath, potential=False, intemplate=MD_npt_track_MSD,
                             inparam=inparam)
    output = parse_lammps_thermo(outfile=output_file)

    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    print(output[0])
    print(output[-1])
    print(output[4])

    return output


def Si_3x3x3_supercell_Li_MSD_vs_time(Tstart, Tstop, numT):
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 6))
    for T in np.linspace(Tstart, Tstop, numT):
        output = Si_3x3x3_supercell_run_MD(T, 0.005, 1000)
        ax_left.plot(output[0], output[4])
        ax_right.plot(output[0], output[-1], label=str(T))
    ax_left.set_ylabel('Temperature')
    ax_right.set_ylabel('Li MSD')
    plt.legend()
    fig.savefig('Si_3x3x3_supercell_MD/Temp_and_MSD-' + time.strftime('%Y%m%d-%H%M%S'))


if __name__ == "__main__":
    Si_3x3x3_supercell_Li_MSD_vs_time(1400, 1800, 3)
