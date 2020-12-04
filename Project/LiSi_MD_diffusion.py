from LiSi_structures import *
from LiSi_LAMMPS_templates import *
import time

project_full_path = '/home/modeler/RemLabs/Project'

def Li_diffusion_constant(temperature, timestep, nsteps):
    """
    Runs an npt ensemble to track the MSD of a Li in a 3x3x3 Si supercell to calculate the diffusion coefficient
    """
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    path = os.path.join(project_full_path, 'Li_diffusion_constant_T=' + str(temperature), timestamp)

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

    # [simtime, pe, ke, energy, temp, press, dens, msd] = output

    return output



if __name__ == "__main__":
    Li_diffusion_constant(600, 0.005, 5000)
    Li_diffusion_constant(700, 0.005, 5000)
    Li_diffusion_constant(800, 0.005, 5000)
    Li_diffusion_constant(900, 0.005, 5000)
