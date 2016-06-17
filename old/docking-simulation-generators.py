
# coding: utf-8

# In[ ]:

from htmd import *
'''
# # Using docking to generate starting poses for simulations
# 
# Download the files for this tutorial from this [link](http://pub.htmd.org/nc983hu3brda/ethtryp.tar.gz)

# ## Dock the protein with the ligand

# In[ ]:

prot = Molecule('ethtryp/trypsin.pdb')
prot.center()
lig = Molecule('ethtryp/ethanol.pdb')
print(lig,prot)
poses, scores = dock(prot, lig)


# ### Visualize the docked poses

# In[ ]:

mol = Molecule()
mol.append(prot)
for i, p in enumerate(poses):
    mol.append(p)
#mol.view(sel='protein', style='NewCartoon', hold=True)
#mol.view(sel='resname MOL', style='Licorice', color=1)


# ## Build systems from docked poses

# In[ ]:
moltbuilt=[]
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
    mol.append(ligand)
    
    D = D + 15
    smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
    topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf', './ethtryp/ethanol.rtf']
    params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', './ethtryp/ethanol.prm']

    moltbuilt.append(charmm.build(smol, topo=topos, param=params, outdir='./docked/build/{}/'.format(i+1), saltconc=0.15))
    if i==4: # For time purposes lets only build the two first
        break


# ## Equilibrate the build systems

# In[ ]:

from htmd.protocols.equilibration_v1 import Equilibration
from natsort import natsorted
md = Equilibration()
md.numsteps = 1000
md.temperature = 298
builds=natsorted(glob('docked/build/*/'))
for i,b in enumerate(builds):
    md.write(b,'docked/equil/{}/'.format(i+1))

mdx = AcemdLocal()
mdx.submit(glob('./docked/equil/*/'))
mdx.wait()
'''
#PRODUCTION OF GENERATORS. IT HAS WORKED!!!!!!!

#all inside loop:
from htmd import *
from natsort import natsorted
from htmd.protocols.production_v1 import Production
equils=natsorted(glob('docked/equil/*/'))
for i,b in enumerate(equils):
	md= Production()
	md.acemd.bincoordinates = 'output.coor'
	md.acemd.extendedsystem  = 'output.xsc'
	md.acemd.binvelocities=None
	md.acemd.binindex=None
	md.acemd.run='50ns'
	md.temperature = 300
	equils=natsorted(glob('docked/equil/*/'))
	md.write('./docked/equil/{}/'.format(i+1), 'docked/generators/{}/'.format(i+1))

mdx = AcemdLocal()
mdx.submit(glob('./docked/generators/*/'))
mdx.wait()




#ADAPTIVE RUN

md = AdaptiveRun()
md.nmin=4
md.nmax=8
md.nepochs = 12
md.app = AcemdLocal()
md.generatorspath='generators/*/'
md.dryrun = True  # creates everything but does not submit anything
md.metricsel1 = 'name CA'
md.metricsel2 = '(resname BEN) and ((name C7) or (name C6))' 
md.metrictype = 'contacts'
md.ticadim = 3 #metricticadmin?
md.updateperiod = 14400 #14400 # execute every 4 hours
md.run()

