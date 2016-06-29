import os
from htmd import *
from htmd.molecule.util import maxDistance
from htmd.protocols.equilibration_v1 import Equilibration
from htmd.protocols.production_v1 import Production
from natsort import natsorted
import sys
import argparse
import random

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
required=True,
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

parser.add_argument('-c', '--config',
dest='config',
action='store',
default='./parameters.config',
required=False,
help='Parameters configuration file')

parser.add_argument('-mol2', '--mol2',
dest='mol2',
action='store',
default=None,
required=False,
help='mol2 file to generate rtf and prm files')

args = parser.parse_args()

def check_arguments():
    if not args.prot:
        sys.stderr.write("Error: You forget to put the protein file path\n")
        exit(1)
    if args.ligand:
        if args.params and args.rtf:
            if args.mol2:
                sys.stderr.write("Error: You Introduce both options: mol2 and pdb,rtf,prm files. Choose only one option\n")
                exit(1)
            ligand_path = args.ligand
            rtf_path = args.rtf
            params_path = args.params
        else:
            sys.stderr.write("Error: You introduce a ligand pdb file, but rtf and prm files are missing. Introduce them with -rtf and -prm input options \n")
    if not args.ligand or not args.params or not args.rtf:
        if not args.mol2:
            sys.stderr.write("You need to introduce one ligand input options: a mol2 file, or  pdb,rtf and prm files.\n")
            exit(1)
        if args.mol2:
            (ligand_path,rtf_path,params_path)=parameter(args.mol2, netcharge)
    return(ligand_path,rtf_path,params_path)

def parse_config (config_file):
    op_config = open(config_file, "r")
    for line in op_config:
        if line.startswith("nbuilds"):
            nbuilds = line.split("\t")[1].strip()
        if line.startswith("minsim"):
            minsim = line.split("\t")[1].strip()
        if line.startswith("maxsim"):
            maxsim = line.split("\t")[1].strip()
        if line.startswith("run_time"):
            run_time = line.split("\t")[1].strip()
        if line.startswith("numbep"):
            numbep = line.split("\t")[1].strip()
        if line.startswith("dimtica"):
            dimtica = line.split("\t")[1].strip()
        if line.startswith("sleeping"):
            sleeping = line.split("\t")[1].strip()
        if line.startswith("netcharge"):
            netcharge = line.split("\t")[1].strip() 
            print(netcharge)
    return(nbuilds, run_time, minsim, maxsim, numbep, dimtica, sleeping, netcharge)

def parameter(mol2, netcharge):
    molec = Molecule(mol2)
    config = ParameterizationConfig()
    config.FileName = mol2
    molec_name = mol2
    molec_name = molec_name.split(".")[0]
    config.JobName = molec_name.split("/")[-1]+str(random.randint(1,1000))
    config.NetCharge = netcharge
    param = Parameterization(config=config)
    paramfiles = param.getParameters(outdir='./', outname=molec_name)
    ligand_path = molec_name+".pdb"
    params_path = molec_name+".prm"
    rtf_path = molec_name+".rtf"
    print('Parameters needed has been created on your working directory.\n Use: python3 SimulationModule.py -p trypsin.pdb --ligand '+ligand_path + ' -rtf '+rtf_path +' -prm '+ params_path+'\n')
    quit()
    

def dockinit(protein_path, ligand_path):
    prot = Molecule(protein_path) 
    prot.filter('protein or water or resname CA')
    prot.set('segid', 'P', sel='protein and noh')
    prot.set('segid', 'W', sel='water')
    prot.set('segid', 'CA', sel='resname CA')
    D = maxDistance(prot, 'all')
    D = D + 15
    prot.center()
    lig = Molecule(ligand_path)
    poses, scores = dock(prot, lig)
    return (prot, poses, D)

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
        print(path_ligand_rtf,path_ligand_rpm)
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

if __name__ == "__main__":
    if len(glob('./docked/'))!=0:
        sys.stderr.write ('A folder called docked already exists, please change its name to avoid overwritting \n')
        quit()
    (nbuilds, run_time, minsim, maxsim, numbep, dimtica, sleeping, netcharge) = parse_config(args.config)
    (ligand_path,rtf_path,params_path)=check_arguments()
    (prot, poses, D) = dockinit(args.prot, ligand_path)
    sys.stderr.write('\nDocking finished.\n')
    building(prot,poses,D,rtf_path,params_path,nbuilds)
    sys.stderr.write('\nAll systems build.\n')
    Equilibrate()
    sys.stderr.write('\nAll systems equilibrated.Entering production, this could take days of running...\n')
    Produce(run_time)
    sys.stderr.write('\nFinished producing. Starting the adaptive run, this could take days of running...\n')
    adaptive(minsim,maxsim,numbep,dimtica,sleeping)
