'''
Post-processing based on nanoAOD tools
'''

import os

# Import tools
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

# Import modules
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer, pufile_data
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lumiWeightProducer import lumiWeightProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer import METSigProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.private.privateModule import privateProducer

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel', action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--samples', action='store', nargs='*', type=str, default=['TTZToLLNuNu_ext'], help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--skim', action='store', nargs='?', type=str, default='dimuon', help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--job', action='store', type=int, default=0, help="Run only jobs i" )
argParser.add_argument('--nJobs', action='store', nargs='?', type=int,default=1, help="Maximum number of simultaneous jobs.")
argParser.add_argument('--prepare', action='store_true', help="Prepare, don't acutally run" )
#argParser.add_argument(
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)

# Import samples

from Samples.nanoAOD.Summer16 import allSamples as Summer16
from Samples.nanoAOD.Summer16_private import allSamples as Summer16_priv

allSamples = Summer16_priv

samples = []
for selectedSample in options.samples:
    for sample in allSamples:
        if selectedSample == sample.name:
            samples.append(sample)
            sample.normalization = float(sample.normalization)

if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample = Sample.combine(sample_name, samples, maxN = None)
    # Clean up
    for s in samples:
        sample.clear()
elif len(samples)==1:
    sample = samples[0]
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

# calculate the lumi scale factor for the weight
targetLumi = 1000.
if sample.isData:
    lumiScaleFactor = 1
else:
    lumiScaleFactor = sample.xSection * targetLumi / float(sample.normalization)

# split into 1 job per file.
len_orig = len(sample.files)
logger.info("Sample has %s files", len_orig)
sample = sample.split( n=len_orig, nSub=options.job)

# Put together skim
isDiMuon        = options.skim.lower().startswith('dimuon')
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')

skimConds = []
if isDiMuon:
    skimConds.append( "Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=2" )
elif isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=2" )
elif isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

cut = '&&'.join(skimConds)

# Main part

directory = "/afs/hephy.at/data/dspitzbart03/nanoAOD/"
output_directory = os.path.join( directory, options.skim, sample.name )

modules = [
    puWeightProducer("auto",pufile_data,"pu_mc","pileup",verbose=False),
    lumiWeightProducer(lumiScaleFactor),
    privateProducer(),
    METSigProducer("Summer16_25nsV1_MC", [1.39,1.26,1.21,1.23,1.28,-0.26,0.62]),
    # MET significance producer
]

p = PostProcessor(output_directory,sample.files,cut=cut, modules=modules)
if not options.prepare:
    p.run()


