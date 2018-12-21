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
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer       import puWeightProducer, pufile_data, pufile_mc
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer            import METSigProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.private.METSigTools           import METSigTools
from PhysicsTools.NanoAODTools.postprocessing.modules.private.lumiWeightProducer    import lumiWeightProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.private.applyJSON             import applyJSON

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel', action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--samples', action='store', nargs='*', type=str, default=['TTZToLLNuNu_ext'], help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--skim', action='store', nargs='?', type=str, default='dimuon', help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--job', action='store', type=int, default=0, help="Run only jobs i" )
argParser.add_argument('--nJobs', action='store', nargs='?', type=int,default=1, help="Maximum number of simultaneous jobs.")
argParser.add_argument('--prepare', action='store_true', help="Prepare, don't acutally run" )
argParser.add_argument('--year', action='store', default=None, help="Which year? Important for json file.")
#argParser.add_argument(
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)

# Import samples

#from Samples.nanoAOD.Summer16 import allSamples as Summer16
from Samples.nanoAOD.Summer16_METSig_v1 import allSamples as Summer16_METSig_v1
from Samples.nanoAOD.Run2018_26Sep2018_private import allSamples as Run2018_26Sep2018_private
from Samples.nanoAOD.Run2016_17Jul2018_private import allSamples as Run2016_17Jul2018_private

allSamples = Summer16_METSig_v1 + Run2018_26Sep2018_private + Run2016_17Jul2018_private

print "Searching for sample %s"%options.samples[0]  

samples = []
for selectedSample in options.samples:
    for sample in allSamples:
        if selectedSample == sample.name:
            samples.append(sample)
            sample.normalization = float(sample.normalization)

print "Working on sample %s"%sample.name

year = int(options.year)
if year == 2016:
    json_file = "$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
elif year == 2017:
    json_file = "$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
elif year == 2018:
    json_file = "$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt"
else:
    raise ValueError( "Please define a supported year." )

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

directory = "/afs/hephy.at/data/dspitzbart03/nanoSamples/2016_v4/"
output_directory = os.path.join( directory, options.skim, sample.name )

logger.info("Loading modules.")

if sample.isData:
    modules = [
        METSigTools(),
        lumiWeightProducer(1, isData=True),
        #METSigProducer("Summer16_25nsV1_MC", [1.39,1.26,1.21,1.23,1.28,-0.26,0.62]),
        applyJSON(json_file)
        # MET significance producer
    ]

else:
    modules = [
        puWeightProducer(pufile_mc,pufile_data,"pu_mc","pileup",verbose=False),
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


