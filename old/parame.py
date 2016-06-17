from htmd import *
from htmd.parameterize import Configuration, Parameterisation
glutamate = Molecule('./ethanol.mol2')
config = Configuration()
config.FileName = "./ethanol.mol2"
config.JobName = "glu-1"
config.NetCharge = "1"
param = Parameterisation(config=config)
paramfiles = param.getParameters()
shutil.copyfile(paramfiles['RTF'], "glutamate.rtf")
shutil.copyfile(paramfiles['PRM'], "glutamate.prm")
shutil.copyfile(paramfiles['PDB'], "glutamate-final.pdb")
