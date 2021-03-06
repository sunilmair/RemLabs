from labutil.src.plugins.pwscf import *
from ase import Atoms
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
    #fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure


def make_bcc_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Fe', 'bcc', a=alat)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    #fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure


def make_anti_bcc_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = Atoms('Fe2', [(0, 0, 0), (alat/2, alat/2, alat/2)], cell=[alat, alat, alat, 90, 90, 90], pbc=True)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    #fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure

def make_hcp_struc(alat, clat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Fe', 'hcp', a=alat, c=clat)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    #fecell.set_atomic_numbers([26, 27])
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


def compute_hcp_energy(alat, clat, nk, ecut):
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
    struc = make_hcp_struc(alat=alat, clat=clat)
    kpts = Kpoints(gridsize=[2*nk, 2*nk, nk], option='automatic', offset=False)
    dirname = 'Fe_hcp_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1B/hcpv2", dirname))
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


def compute_bcc_energy(alat, nk, ecut):
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
    struc = make_bcc_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_bcc_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1B/bcc", dirname))
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


def compute_bcc_single_atom_mag_energy(alat, nk, ecut, mag_state):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    if (mag_state == 'non'):
        mag_start = 0
    elif (mag_state == 'ferro'):
        mag_start = 1
    else:
        mag_start = 0.7

    potname = 'Fe.pbe-nd-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Fe': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname),
                  'Co': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname)
                  }
    struc = make_bcc_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_'+mag_state+'_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1C/"+"mag_state", dirname))
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
            'starting_magnetization(1)': mag_start,
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

def compute_bcc_anti_mag_energy(alat, nk, ecut):
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
    struc = make_anti_bcc_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_anti_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1C/anti", dirname))
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
            'starting_magnetization(1)': 1.0,
            'starting_magnetization(2)': -1.0,
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
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1a_hcp_conv_v1", dirname))
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


def relax_bcc_energy(alat, nk, ecut, forc_conv_thr, press_conv_thr):
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
    struc = make_bcc_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1_bcc", dirname))
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


def relax_scan_bcc_nk_conv():
    nk_list = [14, 16, 18]
    ecut = 30
    alat = 2.86

    forc_conv_thr = 0.001
    press_conv_thr = 0.5

    output = [relax_bcc_energy(alat=alat, ecut=ecut, nk=nk, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr) for nk in nk_list]
    print(output)


def volume_scan_hcp():
    nk = 7
    ecut = 30

    c_over_a = 4.393156594/2.542281076
    equil_a = 2.542281076

    max_a_multiplier = 1.05
    min_a_multiplier = 0.85

    alat_list = numpy.linspace(min_a_multiplier*equil_a, max_a_multiplier*equil_a, 6)

    output = [compute_hcp_energy(alat, c_over_a*alat, nk, ecut) for alat in alat_list]
    print(output)


def volume_scan_bcc():
    nk = 14
    ecut = 30

    equil_a = 2.845479126

    max_a_multiplier = 1.05
    min_a_multiplier = 0.9

    alat_list = numpy.linspace(min_a_multiplier*equil_a, max_a_multiplier*equil_a, 5)

    output = [compute_bcc_energy(alat, nk, ecut) for alat in alat_list]
    print(output)


def mag_evaluate():
    nk = 14
    ecut = 30

    equil_a = 2.845479126

    output = [compute_bcc_single_atom_mag_energy(equil_a, nk, ecut, 'default'),
              compute_bcc_single_atom_mag_energy(equil_a, nk, ecut, 'non'),
              compute_bcc_single_atom_mag_energy(equil_a, nk, ecut, 'ferro'),
              compute_bcc_anti_mag_energy(equil_a, nk, ecut)]

    print(output)

if __name__ == '__main__':
    # put here the function that you actually want to run
    #relax_scan()
    #relax_scan_nk_conv()
    volume_scan_hcp()
