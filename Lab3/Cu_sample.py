from labutil.src.plugins.pwscf import *
from ase import Atoms
from ase.spacegroup import crystal
from ase.io import write
from ase.build import *
import numpy
import matplotlib.pyplot as plt



def make_fcc_struc(element, alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    unit_cell = bulk(element, 'fcc', a=alat)
    # check how your cell looks like
    #write('s.cif', unit_cell)
    print(unit_cell.cell)
    structure = Struc(ase2struc(unit_cell))
    return structure

def make_CuAu_struc(alat, clat):
    unit_cell = Atoms('CuAu', [(0,0,0), (alat/2, alat/2, clat/2)], cell=[alat, alat, clat, 90, 90, 90], pbc=True)
    structure = Struc(ase2struc(unit_cell))
    return structure


def relax_Cu_energy(alat, nk, ecut, forc_conv_thr, press_conv_thr):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    print(nk)
    pseudopots = {'Cu': PseudoPotential(ptype='uspp', element='Cu', functional='LDA', name='Cu.pz-d-rrkjus.UPF')}
    struc = make_fcc_struc(element='Cu', alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Cu_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem3a/Cu2", dirname))
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


def relax_Au_energy(alat, nk, ecut, forc_conv_thr, press_conv_thr):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    print(nk)
    pseudopots = {'Au': PseudoPotential(ptype='uspp', element='Au', functional='LDA', name='Au.pz-d-rrkjus.UPF')}
    struc = make_fcc_struc(element='Au', alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Au_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem3a/Au", dirname))
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


def relax_CuAu_energy(alat, clat, nk, ecut, forc_conv_thr, press_conv_thr):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {'Cu': PseudoPotential(ptype='uspp', element='Cu', functional='LDA', name='Cu.pz-d-rrkjus.UPF'),
                 'Au': PseudoPotential(ptype='uspp', element='Au', functional='LDA', name='Au.pz-d-rrkjus.UPF')}
    struc = make_CuAu_struc(alat=alat, clat=clat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'CuAu_a_{}_c_{}_ecut_{}_nk_{}'.format(alat, clat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem3c/CuAu_relaxv1", dirname))
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


def compute_CuAu_energy(alat, clat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    print(nk)
    pseudopots = {'Cu': PseudoPotential(ptype='uspp', element='Cu', functional='LDA', name='Cu.pz-d-rrkjus.UPF'),
                 'Au': PseudoPotential(ptype='uspp', element='Au', functional='LDA', name='Au.pz-d-rrkjus.UPF')}
    struc = make_CuAu_struc(alat=alat, clat=clat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'CuAu_a_{}_c_{}_ecut_{}_nk_{}'.format(alat, clat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem3c/CuAu_convv1", dirname))
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
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def problem3a_Cu_relax():
    alat = numpy.sqrt(2)*2.561
    nk_list = numpy.arange(4, 16, 2)
    ecut = 40

    forc_conv_thr = 0.001
    press_conv_thr = 0.1

    output = [relax_Cu_energy(alat=alat, nk=nk, ecut=ecut, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr)['energy']
              for nk in nk_list]

    print(output)


def problem3a_Au_relax():
    alat = numpy.sqrt(2)*2.95
    nk_list = numpy.arange(4, 16, 2)
    ecut = 40

    forc_conv_thr = 0.001
    press_conv_thr = 0.1

    output = [relax_Au_energy(alat=alat, nk=nk, ecut=ecut, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr)['energy']
              for nk in nk_list]

    print(output)


def problem3c_CuAu_relax():
    average_parameter = (3.55078382 + 4.047458838)/2
    alat = average_parameter/numpy.sqrt(2)
    clat = average_parameter
    nk = 12
    ecut = 40

    forc_conv_thr = 0.001
    press_conv_thr = 0.1

    output = relax_CuAu_energy(alat=alat, clat=clat, nk=nk, ecut=ecut, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr)
    print(output)


def problem3c_CuAu_nk_conv():
    alat = 2.809846895
    clat = 3.527086544
    nk_list = numpy.arange(4, 16, 2)
    ecut = 40

    output = [compute_CuAu_energy(alat=alat, clat=clat, nk=nk, ecut=ecut)['energy'] for nk in nk_list]
    print(output)


if __name__ == '__main__':
    #problem3a_Cu_relax()
    #problem3a_Au_relax()
    #problem3c_CuAu_relax()
    problem3c_CuAu_nk_conv()
