
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
