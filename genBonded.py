from mbuild.compound import Compound
from mbuild.formats.gromacs import save_gromacs

system = Compound.load('temp.mol2')
system = system.to_trajectory()

system.top.find_forcefield_terms()
save_gromacs(system)
