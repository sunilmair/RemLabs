from labutil.src.plugins.pwscf import *
from ase.io import write
from ase import Atoms
import matplotlib.pyplot as plt



def make_struc(alat, displacement=0):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    lattice = alat * numpy.identity(3)
    symbols = ['Pb', 'Ti', 'O', 'O', 'O']
    sc_pos = [[0,0,0], [0.5,0.5,0.5 + displacement], [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0]]
    perov = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(perov))
    return structure



def compute_energy(alat, nk, ecut, displ=0):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {'Pb': PseudoPotential(ptype='uspp', element='Pb', functional='LDA', name='Pb.pz-d-van.UPF'),
                  'Ti': PseudoPotential(ptype='uspp', element='Ti', functional='LDA', name='Ti.pz-sp-van_ak.UPF'),
                  'O': PseudoPotential(ptype='uspp', element='O', functional='LDA', name='O.pz-rrkjus.UPF')}
    struc = make_struc(alat=alat, displacement=displ)
    # fix the Pb and Ti atoms in place during relaxation
    constraint = Constraint(atoms={'0': [0,0,0], '1': [0,0,0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=True)
    dirname = 'PbTiO3_a_{}_ecut_{}_nk_{}_displ_{}'.format(alat, ecut, nk, displ)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem2/a/lat1", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none'
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.7,
            'conv_thr': 1e-7,
        },
        'IONS': {
            'ion_dynamics': 'bfgs'
        },
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, constraint=constraint, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def compute_Ti_disp_energy(alat, nk, ecut, displ=0):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    #print(displ)
    pseudopots = {'Pb': PseudoPotential(ptype='uspp', element='Pb', functional='LDA', name='Pb.pz-d-van.UPF'),
                  'Ti': PseudoPotential(ptype='uspp', element='Ti', functional='LDA', name='Ti.pz-sp-van_ak.UPF'),
                  'O': PseudoPotential(ptype='uspp', element='O', functional='LDA', name='O.pz-rrkjus.UPF')}
    struc = make_struc(alat=alat, displacement=displ)
    # fix the Pb and Ti atoms in place during relaxation
    constraint = Constraint(atoms={'0': [0,0,0], '1': [0,0,0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=True)
    dirname = 'PbTiO3_a_{}_ecut_{}_nk_{}_displ_{}'.format(alat, ecut, nk, displ)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem2/b1", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'relax',
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none'
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.7,
            'conv_thr': 1e-7,
        },
        'IONS': {
            'ion_dynamics': 'bfgs'
        },
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, constraint=constraint, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def compute_Ti_relax(alat, nk, ecut, displ=0):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {'Pb': PseudoPotential(ptype='uspp', element='Pb', functional='LDA', name='Pb.pz-d-van.UPF'),
                  'Ti': PseudoPotential(ptype='uspp', element='Ti', functional='LDA', name='Ti.pz-sp-van_ak.UPF'),
                  'O': PseudoPotential(ptype='uspp', element='O', functional='LDA', name='O.pz-rrkjus.UPF')}
    struc = make_struc(alat=alat, displacement=displ)
    # fix the Pb atoms in place during relaxation
    constraint = Constraint(atoms={'0': [0,0,0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=True)
    dirname = 'PbTiO3_a_{}_ecut_{}_nk_{}_displ_{}'.format(alat, ecut, nk, displ)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem2/c", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'relax',
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none'
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.7,
            'conv_thr': 1e-7,
        },
        'IONS': {
            'ion_dynamics': 'bfgs'
        },
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, constraint=constraint, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def lattice_scan():
    nk = 4
    ecut = 30
    alat_list = numpy.linspace(3.8, 4.0, 11)
    print(alat_list)
    energy_list = []
    for alat in alat_list:
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        energy_list.append(output['energy'])
        print(output)
    print(alat_list)
    print(energy_list)
    plt.plot(alat_list, energy_list)
    plt.show()


def problem2_ecut_conv():
    nk = 4
    ecut_list = numpy.arange(45.0, 100.0, 10)
    alat = 3.97
    output = [compute_energy(alat=alat, ecut=ecut, nk=nk)['energy'] for ecut in ecut_list]

    print(ecut_list)
    print(output)


def problem2a_lattice_scan():
    nk = 4
    ecut = 90
    alat_list = numpy.linspace(3.8, 4.0, 11)

    output = [compute_energy(alat=alat, ecut=ecut, nk=nk)['energy'] for alat in alat_list]

    print(alat_list)
    print(output)


def problem2b_Ti_disp_scan():
    nk = 4
    ecut = 70
    alat = 3.88
    disp_list = numpy.linspace(0, 0.05, 11)

    output = [compute_Ti_disp_energy(alat=alat, ecut=ecut, nk=nk, displ=disp)['energy'] for disp in disp_list]

    print(disp_list)
    print(output)


def problem2c_Ti_relax():
    nk = 4
    ecut = 90  # change based on conv
    alat = 3.88  # change based on eq
    disp = 0.0275 # change based on eq

    output = compute_Ti_relax(alat=alat, ecut=ecut, nk=nk, displ=disp)
    print(output)

if __name__ == '__main__':
    # put here the function that you actually want to run
    #problem2_ecut_conv()
    #problem2a_lattice_scan()
    #problem2b_Ti_disp_scan()
    problem2c_Ti_relax()
