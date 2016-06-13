
# coding: utf-8

from htmd import *

'''
prot = Molecule('ethtryp/trypsin.pdb')
prot.center()
lig = Molecule('ethtryp/ethanol.pdb')
poses, scores = dock(prot, lig)


# ### Visualize the docked poses

mol = Molecule()
mol.append(prot)
for i, p in enumerate(poses):
    mol.append(p)
mol.view(sel='protein', style='NewCartoon', hold=True)
mol.view(sel='resname MOL', style='Licorice', color=1)


# ## Build systems from docked poses


molbuilt = []
for i, p in enumerate(poses):
    prot = Molecule('ethtryp/trypsin.pdb')
    prot.filter('chain A and (protein or water or resname CA)')
    prot.set('segid', 'P', sel='protein and noh')
    prot.set('segid', 'W', sel='water')
    prot.set('segid', 'CA', sel='resname CA')
    prot.center()
    from htmd.molecule.util import maxDistance
    D = maxDistance(prot, 'all')
    
    ligand = p
    ligand.set('segid','L')
    ligand.set('resname','MOL')
    
    mol = Molecule(name='combo')
    mol.append(prot)
    mol.append(ligand) #ligand does not have to rotate nor move as it is already docked.
    
    D = D + 15
    smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
    topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf', './ethtryp/ethanol.rtf']
    params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', './ethtryp/ethanol.prm']

    molbuilt.append(charmm.build(smol, topo=topos, param=params, outdir='./docked/build/{}/'.format(i+1), saltconc=0.15))
    if i==1: # For time purposes lets only build the two first
        break





# ## Equilibrate the build systems

from htmd.protocols.equilibration_v1 import Equilibration
from natsort import natsorted
md = Equilibration()
md.numsteps = 1000
md.temperature = 298

builds = natsorted(glob('./docked/build/*/'))
for i, b in enumerate(builds):
    md.write(b, 'docked/equil/{}/'.format(i+1))

mdx = AcemdLocal()
mdx.submit(glob('./docked/equil/*/'))
mdx.wait()


# ## Create the production folder

from htmd.protocols.production_v1 import Production
from natsort import natsorted

adapt = Production()
adapt.acemd.bincoordinates = 'output.coor'
adapt.acemd.extendedsystem  = 'output.xsc'
adapt.acemd.binvelocities=None
adapt.acemd.binindex=None
adapt.acemd.run='50ns'
adapt.temperature = 300
equils = natsorted(glob('docked/equil/*/'))
print(equils)
for i, b in enumerate(equils):
    md.write(b, 'docked/generators/{}/'.format(i+1))

mdx = AcemdLocal()
mdx.submit(glob('./docked/generators/*/'))
mdx.wait()
'''

adapt = AdaptiveRun()
print ('AdaptiveRun has already started')
adapt.project='projecttrial3'
adapt.nmin = 2
adapt.nmax = 4
adapt.nepochs = 30
adapt.ticadim = 0
adapt.metricsel1 = 'name CA'
print('give the input directory')
adapt.generatorspath = './docked/generators/2/'
adapt.inputpath='./input6' #this name has to change with each run
adapt.resultspath='./adaprunresults'
print('start the proccess')
adapt.app = AcemdLocal()
adapt.run()

