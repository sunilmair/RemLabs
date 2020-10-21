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
    nk_list = numpy.arange(4, 18, 2)
    ecut = 40

    forc_conv_thr = 0.001
    press_conv_thr = 0.1

    output = [relax_Au_energy(alat=alat, nk=nk, ecut=ecut, forc_conv_thr=forc_conv_thr, press_conv_thr=press_conv_thr)['energy']
              for nk in nk_list]

    print(output)


if __name__ == '__main__':
    #problem3a_Cu_relax()
    problem3a_Au_relax()
