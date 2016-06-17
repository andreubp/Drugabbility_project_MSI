from htmd import *
from htmd.molecule.util import uniformRandomRotation


glob('parametering/*')

print ("Adding Molecule...")
prot = Molecule('parametering/trypsin.pdb')
prot.filter('chain A and (protein or water or resname CA)')
prot.set('segid', 'P', sel='protein')
prot.set('segid', 'W', sel='water')
prot.set('segid', 'CA', sel='resname CA')

prot.center()

print ("Adding Ligand...")

ligand = Molecule('parametering/benzamidine.pdb')
ligand.center()

ligand.rotateBy(uniformRandomRotation())
from htmd.molecule.util import maxDistance
D = maxDistance(prot, 'all')
D += 10

ligand.moveBy([D, 0, 0])
ligand.rotateBy(uniformRandomRotation())

ligand.set('segid','L')
ligand.set('resname','MOL')

mol = Molecule(name='combo')
mol.append(prot)
mol.append(ligand)
mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
mol.reps.add(sel='resname MOL', style='Licorice')

print("solvating...")
D += 5
smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
smol.reps.add(sel='water', style='Lines')

charmm.listFiles()

topos  = ['parametering/benzamidine.rtf','top/top_all36_prot.rtf','top/top_all36_lipid.rtf','top/top_water_ions.rtf']
params = ['parametering/benzamidine.prm','par/par_all36_prot.prm','par/par_all36_lipid.prm','par/par_water_ions.prm']

#topos  = ['top/top_all22star_prot.rtf', 'parametering/benzamidine.rtf']
#params = ['par/par_all22star_prot.prm', 'parametering/benzamidine.prm']

molbuilt = charmm.build(smol, topo=topos, param=params, outdir='/tmp/build')
