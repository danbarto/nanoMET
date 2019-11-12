import copy, os, sys
from RootTools.core.Sample import Sample 
import ROOT

## Logging
import logging
logger = logging.getLogger(__name__)

# Data directory
try:    data_directory = sys.modules['__main__'].data_directory
except: from nanoMET.tools.user import data_directory

# Take post processing directory if defined in main module
try:    postProcessing_directory = sys.modules['__main__'].postProcessing_directory
except: postProcessing_directory = '2017_v21/dimuon/'

logger.info("Loading data samples from directory %s", os.path.join(data_directory, postProcessing_directory))

dirs = {}
for (r, version) in [('B',''),('C',''),('D',''),('E',''),('F','')]:
    runTag = 'Run2017' + r + '_31Mar2018' + version
    dirs["DoubleMuon_Run2017"       + r + version ] = ["DoubleMuon_"        + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for r in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + r])

for pd in ['DoubleMuon']:
    merge(pd, 'Run2017',    ['Run2017B', 'Run2017C', 'Run2017D', 'Run2017E', 'Run2017F'])
    merge(pd, 'Run2017BCDE',    ['Run2017B', 'Run2017C', 'Run2017D', 'Run2017E'])

for key in dirs:
    dirs[key] = [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]]


def getSample(pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName])
    sample.lumi = lumi
    return sample

DoubleMuon_Run2017              = getSample('DoubleMuon',       'Run2017',       (41.53)*1000)
DoubleMuon_Run2017F             = getSample('DoubleMuon',      'Run2017F',       (13.5)*1000)
DoubleMuon_Run2017BCDE          = getSample('DoubleMuon',      'Run2017BCDE',    (28.0)*1000)

allSamples_Data25ns = []
allSamples_Data25ns += [DoubleMuon_Run2017, DoubleMuon_Run2017F, DoubleMuon_Run2017BCDE]

Run2017 = Sample.combine("Run2017", [DoubleMuon_Run2017], texName = "Data")
Run2017.lumi = (41.53)*1000

Run2017F = Sample.combine("Run2017F", [DoubleMuon_Run2017F], texName = "Data")
Run2017F.lumi = (13.50)*1000

Run2017BCDE = Sample.combine("Run2017BCDE", [DoubleMuon_Run2017BCDE], texName = "Data")
Run2017BCDE.lumi = (28.0)*1000

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True


