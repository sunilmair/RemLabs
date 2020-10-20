from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
from ase.build import *
import numpy
import matplotlib.pyplot as plt



def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Fe', 'hcp', a=alat)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Fe.pbe-nd-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Fe': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname),
                  'Co': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname)
                  }
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
            'nspin': 2,
            'starting_magnetization(1)': 0.7,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def relax_energy(alat, nk, ecut, forc_conv_thr, press_conv_thr):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Fe.pbe-nd-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Fe': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname),
                  'Co': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname)
                  }
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[2*nk, 2*nk, nk], option='automatic', offset=False)
    dirname = 'Fe_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'vc-relax',
            'forc_conv_thr': forc_conv_thr,
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
            'nspin': 2,
            'starting_magnetization(1)': 0.7,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {
            'ion_dynamics': 'bfgs',
        },
        'CELL': {
            'cell_dofree': 'all',
            'cell_dynamics': 'bfgs',
            'press': 0.0,
            'press_conv_thr': press_conv_thr,
        },
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def lattice_scan():
    nk = 3
    ecut = 30
    alat = 3.0
    output = compute_energy(alat=alat, ecut=ecut, nk=nk)
    print(output)


def relax_scan():
    nk = 3
    ecut = 30
    alat = 3.0

    forc_conv_thr = 0.001
    press_conv_thr = 0.5

    output = relax_energy(alat=alat, ecut=ecut, nk=nk, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr)
    print(output)


def relax_scan_nk_conv():
    nk_list = [3, 4, 5, 6]
    ecut = 30
    alat = 2.5

    forc_conv_thr = 0.001
    press_conv_thr = 0.5

    output = [relax_energy(alat=alat, ecut=ecut, nk=nk, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr) for nk in nk_list]
    print(output)


if __name__ == '__main__':
    # put here the function that you actually want to run
    #lattice_scan()
    relax_scan_nk_conv()
