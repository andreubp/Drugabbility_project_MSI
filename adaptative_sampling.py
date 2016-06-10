from htmd.protocols.production_v1 import Production

adapt = Production()
#adapt.acemd.show()

adapt.acemd.bincoordinates = '/docked/equil/1/output.coor'
adapt.acemd.extendedsystem  = '/docked/equil/1/output.xsc'
adapt.acemd.binvelocities=None
adapt.acemd.binindex=None
adapt.acemd.run='50ns'
adapt.temperature = 300
adapt.write('./equil', './generators/s1')

#adapt = AdaptiveRun()
#adapt.nmin = 2
#adapt.nmax = 3
#adapt.nepochs = 2
#adapt.ticadim = 0
#adapt.metricsel1 = 'name CA'
#adapt.generatorspath = htmd.home()+'/data/dhfr'
#adapt.app = AcemdLocal()
#adapt.run()
