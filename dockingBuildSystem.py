
# coding: utf-8

# In[4]:

from htmd import *
os.chdir('./') # Don't use this command


# # Using docking to generate starting poses for simulations
# 
# Download the files for this tutorial from this [link](http://pub.htmd.org/nc983hu3brda/bentryp.tar.gz)

# ## Dock the protein with the ligand

# In[5]:

prot = Molecule('bentryp/trypsin.pdb')
prot.center()
lig = Molecule('bentryp/benzamidine.pdb')
poses, scores = dock(prot, lig)


# ### Visualize the docked poses

# In[6]:
print('Now it will fail')
mol = Molecule()
mol.append(prot)
for i, p in enumerate(poses):
    mol.append(p)
mol.view(sel='protein', style='NewCartoon', hold=True)
mol.view(sel='resname MOL', style='Licorice', color=1)


# ## Build systems from docked poses

# In[7]:

molbuilt = []
for i, p in enumerate(poses):
    prot = Molecule('bentryp/trypsin.pdb')
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
    topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf', './bentryp/benzamidine.rtf']
    params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', './bentryp/benzamidine.prm']

    molbuilt.append(charmm.build(smol, topo=topos, param=params, outdir='./docked/build/{}/'.format(i+1), saltconc=0.15))
    if i==1: # For time purposes lets only build the two first
        break


# In[11]:

len(molbuilt)


# In[9]:

from ipywidgets.widgets import Box
w = []
for i, m in enumerate(molbuilt):
    #m.view(sel='protein', style='NewCartoon', hold=True)
    #m.view(sel='water', style='Lines', hold=True)
    h = m.view(sel='resname MOL', style='Licorice', color=0)
    w.append(h)
b = Box(children=(w[0],w[1]))
b


# ## Equilibrate the build systems

# In[21]:

from htmd.protocols.equilibration_v1 import Equilibration
from natsort import natsorted
md = Equilibration()
md.numsteps = 1000
md.temperature = 298

builds = natsorted(glob('./docked/build/*/'))
for i, b in enumerate(builds):
    md.write(b, 'docked/equil/{}/'.format(i+1))


# In[22]:

mdx = AcemdLocal()
mdx.submit(glob('./docked/equil/*/'))


# In[ ]:

mdx.wait()


# ## Create the production folder

# In[ ]:

from htmd.protocols.production_v1 import Production
md = Production()
md.acemd.run = '50ns'
md.temperature = 298

equils = natsort(glob('docked/equil/*/'))
for i, b in enumerate(equils):
    md.write(b, 'docked/generators/{}/'.format(i+1))


# In[ ]:

mdx = AcemdLocal()
mdx.submit(glob('./docked/generators/*/'))


# In[ ]:

mdx.wait()

