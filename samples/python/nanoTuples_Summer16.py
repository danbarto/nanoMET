import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

dbFile = '/afs/hephy.at/data/dspitzbart01/nanoAOD/DB_Summer16.sql'
dataDirectory = '/afs/hephy.at/data/dspitzbart01/nanoAOD/Fall17/'

# specify a local directory if you want to create (and afterwards automatically use) a local copy of the sample, otherwise use the grid.
DYJetsToLL_M50_LO = Sample.nanoAODfromDAS('DYJetsToLL_M50_LO','/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dbFile=dbFile, xSection=2008.*3)


DY = [
    DYJetsToLL_M50_LO,
    ]

dummy = [
    ]

allSamples = DY + dummy

for s in allSamples:
    s.isData = False

