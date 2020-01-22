import copy, os, sys
from RootTools.core.Sample import Sample 
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

# Data directory
try:    data_directory = sys.modules['__main__'].data_directory
except: from nanoMET.tools.user import data_directory

# Take post processing directory if defined in main module
try:    postProcessing_directory = sys.modules['__main__'].postProcessing_directory
except: postProcessing_directory = '2018_v22/dimuon/'

logger.info("Loading data samples from directory %s", os.path.join(data_directory, postProcessing_directory))

dirs = {}
for (r, version) in [('A',''),('B',''),('C',''),('D','')]:
    runTag = 'Run2018' + r + '_17Sep2018' + version
    dirs["DoubleMuon_Run2018"       + r + version ] = ["DoubleMuon_"        + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for r in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + r])

for pd in ['DoubleMuon']:
    merge(pd, 'Run2018',    ['Run2018A', 'Run2018B', 'Run2018C', 'Run2018D'])

for key in dirs:
    dirs[key] = [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]]


def getSample(pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName])
    sample.lumi = lumi
    return sample

DoubleMuon_Run2018              = getSample('DoubleMuon',       'Run2018',       (59.97)*1000)

allSamples_Data25ns = []
allSamples_Data25ns += [DoubleMuon_Run2018]

Run2018 = Sample.combine("Run2018", [DoubleMuon_Run2018], texName = "Data")
Run2018.lumi = (59.97)*1000

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True


