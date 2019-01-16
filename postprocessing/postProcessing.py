'''
Post-processing based on nanoAOD tools
'''

import os

# Import tools
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

# Import modules
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer       import puWeightProducer, pufile_data, pufile_mc, pufile_data2017, pufile_data2018
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer            import METSigProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.private.METSigTools           import METSigTools
from PhysicsTools.NanoAODTools.postprocessing.modules.private.lumiWeightProducer    import lumiWeightProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.private.applyJSON             import applyJSON

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel',    action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--samples',     action='store', nargs='*', type=str, default=['TTZToLLNuNu_ext'], help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--skim',        action='store', nargs='?', type=str, default='dimuon', help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--job',         action='store', type=int, default=0, help="Run only jobs i" )
argParser.add_argument('--nJobs',       action='store', nargs='?', type=int,default=1, help="Maximum number of simultaneous jobs.")
argParser.add_argument('--prepare',     action='store_true', help="Prepare, don't acutally run" )
argParser.add_argument('--year',        action='store', default=None, help="Which year? Important for json file.")
argParser.add_argument('--era',         action='store', default="v1", help="Which era/subdirectory?")
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)

# from RootTools
from RootTools.core.standard            import *

# Import samples
year = int(options.year)

allSamples = []
if year == 2016:
    from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as Summer16_private_legacy_v1
    from Samples.nanoAOD.Run2016_17Jul2018_private import allSamples as Run2016_17Jul2018_private
    allSamples = Summer16_private_legacy_v1 + Run2016_17Jul2018_private
elif year == 2017:
    from Samples.nanoAOD.Fall17_private_legacy_v1 import allSamples as Fall17_private_legacy_v1
    from Samples.nanoAOD.Run2017_31Mar2018_private import allSamples as Run2017_31Mar2018_private
    allSamples = Fall17_private_legacy_v1 + Run2017_31Mar2018_private
elif year == 2018:
    from Samples.nanoAOD.Autumn18_private_legacy_v1 import allSamples as Autumn18_private_legacy_v1
    from Samples.nanoAOD.Run2018_17Sep2018_private import allSamples as Run2018_17Sep2018_private
    allSamples = Autumn18_private_legacy_v1 + Run2018_17Sep2018_private

print "Searching for sample %s"%options.samples[0]  

samples = []
for selectedSample in options.samples:
    for sample in allSamples:
        if selectedSample == sample.name:
            samples.append(sample)
            logger.info("Adding sample %s", sample.name)
            logger.info("Sample has normalization %s", sample.normalization)
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
    logger.info("Final normalization is %s", sample.normalization)
elif len(samples)==1:
    sample = samples[0]
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

logger.info("Sample contains %s files", len(sample.files))
sample.files = sorted(sample.files) # in order to avoid some random ordered file list, different in each job

# calculate the lumi scale factor for the weight
targetLumi = 1000.
if sample.isData:
    lumiScaleFactor = 1
else:
    lumiScaleFactor = sample.xSection * targetLumi / float(sample.normalization)

# filebased job splitting
len_orig = len(sample.files)
logger.info("Sample has %s files", len_orig)
json = sample.json # pickup json from sample, as defined in the Sample repository
sample = sample.split( n=options.nJobs, nSub=options.job)

logger.info("Will run over %s files", len(sample.files))

# Put together skim
isDiMuon        = options.skim.lower().startswith('dimuon')
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')

skimConds = []
if isDiMuon:
    skimConds.append( "Sum$(Muon_pt>15&&abs(Muon_eta)<2.5)>=2" )
elif isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=2" )
elif isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

cut = '&&'.join(skimConds)

logger.info("Using selection: %s", cut)

# Main part

directory = "/afs/hephy.at/data/dspitzbart03/nanoSamples/%s_%s/"%(options.year, options.era)
output_directory = os.path.join( directory, options.skim, sample.name )

logger.info("Loading modules.")

if year == 2016:
    puwProducer = puWeightProducer(pufile_mc,pufile_data,"pu_mc","pileup",verbose=False)
elif year == 2017:
    puwProducer = puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)
elif year == 2018:
    puwProducer = puWeightProducer("auto",pufile_data2018,"pu_mc","pileup",verbose=False)


if sample.isData:
    modules = [
        METSigTools(),
        lumiWeightProducer(1, isData=True),
        #METSigProducer("Summer16_25nsV1_MC", [1.39,1.26,1.21,1.23,1.28,-0.26,0.62]),
        applyJSON(json)
        # MET significance producer
    ]

else:
    modules = [
        puwProducer,
        lumiWeightProducer(lumiScaleFactor),
        METSigTools(),
        #METSigProducer("Summer16_25nsV1_MC", [1.39,1.26,1.21,1.23,1.28,-0.26,0.62]),
        applyJSON(None)
        # MET significance producer
    ]

logger.info("Preparing post-processor.")

p = PostProcessor(output_directory,sample.files,cut=cut, modules=modules)
if not options.prepare:
    logger.info("Running...")
    p.run()


