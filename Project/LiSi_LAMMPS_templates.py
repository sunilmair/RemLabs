
relaxation_calculation_template = """
# ---------- 1. Initialize Simulation ----------
units metal
atom_style atomic
dimension 3
boundary p p p
read_data $DATAINPUT

# ---------- 2. Specify Interatomic Potential ----------
pair_style meam/c
pair_coeff * * library.meam Li Si LiSi.meam Li Si

# ---------- 3. Optimization of Atomic Positions Allowing Unit Cell Changes ----------
write_dump all atom $OUTFILE/initial_dump.atom

thermo_style custom step pe lx ly lz press pxx pyy pzz

fix 1 all box/relax iso 0.0 vmax 0.001

min_style cg
minimize 1e-10 1e-10 1000 10000

write_dump all atom $OUTFILE/final_dump.atom

# ---------- 4. Output Variables ----------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
"""

MD_npt = """
# ---------- 1. Initialize simulation ----------
units metal
atom_style atomic
dimension 3
boundary p p p
read_data $DATAINPUT

pair_style meam/c
pair_coeff * * library.meam Li Si LiSi.meam Li Si

group Li type 1

velocity all create $TEMPERATURE 3320 dist gaussian

# ---------- 2. Describe computed properties----------
compute msdli Li msd
thermo_style custom step pe ke etotal temp press density c_msdli[4]
thermo $TOUTPUT

# ---------- 3. Specify ensemble  ----------
fix 1 all npt temp $TEMPERATURE $TEMPERATURE $TDAMP tchain 2 iso 1.0 1.0 1.0 pchain 2

# ---------- 4. Run -------------
timestep $TIMESTEP
run $NSTEPS
"""

MD_equilibrate_npt_track_MSD = """
# ---------- 1. Initialize simulation ----------
units metal
atom_style atomic
dimension 3
boundary p p p
read_data $DATAINPUT

pair_style meam/c
pair_coeff * * library.meam Li Si LiSi.meam Li Si

group Li type 1

velocity all create $TEMPERATURE 3320 dist gaussian

# ---------- 2. Specify equilibration ensemble  ----------
fix 1 all npt temp $TEMPERATURE $TEMPERATURE $TDAMP tchain 2 iso 1.0 1.0 1.0 pchain 2

# ---------- 3. Run equilibration -------------
timestep $TIMESTEP
run $EQUILNSTEPS

# ---------- 4. Describe computed properties----------
compute msdli Li msd
thermo_style custom step pe ke etotal temp press density c_msdli[4]
thermo $TOUTPUT

# ---------- 3. Run production -------------
timestep $TIMESTEP
run $NSTEPS
"""
