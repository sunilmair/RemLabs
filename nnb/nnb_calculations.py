from nnb_structures import *
from labutil.src.plugins.pwscf import *

pseudopots = {'Na': PseudoPotential(element='Na', name='na_pbe_v1.5.uspp.F.UPF'),
              'N': PseudoPotential(element='N', name='N.pbe-n-radius_5.UPF'),
              'H': PseudoPotential(element='H', name='H.pbe-rrkjus_psl.1.0.0.UPF'),
              'B': PseudoPotential(element='B', name='b_pbe_v1.4.uspp.F.UPF')}


def scf_calculation(struc, nk, ecut, dirname):
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], 'RemLabs', 'nnb', dirname))
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
        },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
    })
    output_file = run_qe_pwscf(runpath=runpath, struc=struc, pseudopots=pseudopots,
                               params=input_params, kpoints=kpts)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def relax_calculation(struc, nk, ecut, forc_conv_thr, press_conv_thr, dirname):
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], 'RemLabs', 'nnb', dirname))
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
            'degauss': 0.02,
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

def test_scf_calculation():
    #struc = make_initial_unitcell_state()
    struc = make_initial_2x2x2_neb_state()
    nk = 2
    ecut = 10
    dirname = 'test_scf_2'
    scf_calculation(struc, nk, ecut, dirname)


def test_relax_calculation():
    struc = make_initial_state()
    nk = 2
    ecut = 10
    forc_conv_thr = 0.001
    press_conv_thr = 0.1
    dirname = 'test_relax'
    relax_calculation(struc, nk, ecut, forc_conv_thr, press_conv_thr, dirname)


if __name__ == '__main__':
    test_scf_calculation()
