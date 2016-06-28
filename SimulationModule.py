from htmd import *
from htmd.molecule.util import maxDistance
from htmd.protocols.equilibration_v1 import Equilibration
from htmd.protocols.production_v1 import Production
from natsort import natsorted
import argparse

parser = argparse.ArgumentParser(description="Druggability Project")
parser.add_argument('-l', '--ligand',
dest='ligand',
action='store',
default=None,
required=False,
help='Ligand path')

parser.add_argument('-p', '--prot',
dest='prot',
action='store',
default=None,
required=False,
help='Protein path')

parser.add_argument('-rtf', '--rtf',
dest='rtf',
action='store',
default=None,
required=False,
help='rtf path')

parser.add_argument('-prm', '--prm',
dest='params',
action='store',
default=None,
required=False,
help='Params path')

args = parser.parse_args()

def simulate(pdbpath,ligandpath,path_ligand_rtf,path_ligand_prm,nbuilds=4,run_time=50,minsim=6,maxsim=8,numbep=12,dimtica=3,sleeping=14400):
    prot = Molecule(pdbpath) 
    prot.filter('protein or water or resname CA')
    prot.set('segid', 'P', sel='protein and noh')
    prot.set('segid', 'W', sel='water')
    prot.set('segid', 'CA', sel='resname CA')
    D = maxDistance(prot, 'all')
    D = D + 15
    prot.center()
    lig = Molecule(ligandpath)
    poses, scores = dock(prot, lig)
    print('\nDocking finished.')
    building(prot,poses,D,path_ligand_rtf,path_ligand_prm,nbuilds)
    print('\nAll systems build.')
    Equilibrate()
    print('All systems equilibrated.Entering production, this could take days of running...')
    Produce(run_time)
    print('Finished producing. Starting the adaptive run, this could take days of running...')
    adaptive(minsim,maxsim,numbep,dimtica,sleeping)


def building(prot,poses,D,path_ligand_rtf,path_ligand_prm,nbuilds=4):
    moltbuilt=[]
    for i, p in enumerate(poses):
        ligand = p
        ligand.set('segid','L')
        ligand.set('resname','MOL')
        mol = Molecule(name='combo')
        mol.append(prot)
        mol.append(ligand)

        smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
        topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf',path_ligand_rtf] #'./ethtryp/ethanol.rtf'
        params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', path_ligand_prm] #'./ethtryp/ethanol.prm'

        moltbuilt.append(charmm.build(smol, topo=topos, param=params, outdir='./docked/build/{}/'.format(i+1), saltconc=0.15))
        if i==nbuilds:
            break


def Equilibrate():
    md = Equilibration()
    md.numsteps = 1000
    md.temperature = 298
    builds=natsorted(glob('docked/build/*/'))
    for i,b in enumerate(builds):
        md.write(b,'docked/equil/{}/'.format(i+1))
    mdx = AcemdLocal()
    mdx.submit(glob('./docked/equil/*/'))
    mdx.wait()


def Produce(run_time=50):
    equils=natsorted(glob('docked/equil/*/'))
    for i,b in enumerate(equils):
        md= Production()
        md.acemd.bincoordinates = 'output.coor'
        md.acemd.extendedsystem  = 'output.xsc'
        md.acemd.binvelocities=None
        md.acemd.binindex=None
        md.acemd.run=str(run_time)+'ns'
        md.temperature = 300
        equils=natsorted(glob('docked/equil/*/'))
        md.write('./docked/equil/{}/'.format(i+1), 'docked/generators/{}/'.format(i+1))

    mdx = AcemdLocal()
    mdx.submit(glob('./docked/generators/*/'))
    mdx.wait()

def adaptive(minsim=6,maxsim=8,numbep=12,dimtica=3,sleeping=14400):
    md = AdaptiveRun()
    md.nmin=minsim
    md.nmax=maxsim
    md.nepochs = numbep
    md.app = AcemdLocal()
    md.generatorspath='./docked/generators/'
    md.datapath='./docked/generators/'
    md.inputpath='./docked/generators/'
    md.dryrun = False 
    md.metricsel1 = 'name CA'
    md.metricsel2 = 'resname MOL and noh'
    md.metrictype = 'contacts'
    md.ticadim = dimtica
    md.updateperiod = sleeping
    md.run()

def analysis(boot=0.8,clusters=1000,merge=5):
    sims = simlist(glob('./filtered/*/'), './filtered/filtered.pdb')
    metr = Metric(sims)
    metr.projection(MetricDistance('protein and name CA', 'resname MOL and noh', metric='contacts'))
    data = metr.project()
    data.fstep = 0.1
    data.plotTrajSizes()
    data.dropTraj()
    tica = TICA(data, 10)
    dataTica = tica.project(3)
    dataBoot = dataTica.bootstrap(0.8)
    dataBoot.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=5) #try with dataTica instead of dataBoot
    model = Model(dataBoot) #try with dataTica
    model.plotTimescales() 

    nframes=input('At which time do time scales converge?')
    ntimescales=input('How many different time scales do you see?')
    model.markovModel(int(nframes)*10, int(ntimescales)) 
    mols = model.getStates()

    kin = Kinetics(model, temperature=298, concentration=0.0037)
    goodmacros=dict()
    mols = model.getStates()
    for i in range(len(mols)):
        for j in range(len(mols)):
            r = kin.getRates(i,j)
            if r.g0eq < -2:
                try:
                    goodmacros[i].append(j)
                except:
                    goodmacros[i]=[j]
    
    curr_max_len=0
    for keys in goodmacros:
        if len(goodmacros[keys])>curr_max_len:
            curr_max_len=len(goodmacros[keys])
            thekey=keys

    retlist=list()
    print('These models contain the best interactions/poses:')
    for sinks in goodmacros[thekey]:
        print(mols[sinks])
        retlist.append(mols[sinks])

    kin.plotRates(rates=('g0eq'))
    kin.plotFluxPathways()
    return retlist


simulate(args.prot, args.ligand, args.rtf,args.params,nbuilds=1,run_time=10,minsim=1,maxsim=2,numbep=4,dimtica=3,sleeping=14400)
