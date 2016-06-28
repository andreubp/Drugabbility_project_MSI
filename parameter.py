from htmd import *
from htmd.parameterize import Configuration, Parameterisation
glutamate = Molecule('./choline.mol2')
config = Configuration()
config.FileName = "./choline.mol2"
config.JobName = "cho-2"
config.NetCharge = "1"
param = Parameterisation(config=config)
paramfiles = param.getParameters()
shutil.copyfile(paramfiles['RTF'], "choline.rtf")
shutil.copyfile(paramfiles['PRM'], "choline.prm")
shutil.copyfile(paramfiles['PDB'], "choline.pdb")
