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
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties       import jetmetUncertaintiesProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib                import jetRecalib
from PhysicsTools.NanoAODTools.postprocessing.modules.private.METSigTools           import METSigTools
from PhysicsTools.NanoAODTools.postprocessing.modules.private.METminProducer        import METminProducer
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

def extractEra(sampleName):
    return sampleName[sampleName.find("Run"):sampleName.find("Run")+len('Run2000A')]

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

if sample.isData:
    era = extractEra(samples[0].name)[-1]
else:
    era = None
print "######### Era %s ########"%era

# calculate the lumi scale factor for the weight
targetLumi = 1000.
if sample.isData:
    lumiScaleFactor = 1
else:
    lumiScaleFactor = sample.xSection * targetLumi / float(sample.normalization)

keepSampleName = sample.name # to mitigate Mateusz change of naming convention

# filebased job splitting
len_orig = len(sample.files)
logger.info("Sample has %s files", len_orig)
json = sample.json # pickup json from sample, as defined in the Sample repository
sample = sample.split( n=options.nJobs, nSub=options.job)
sample.name = keepSampleName

logger.info("Will run over %s files", len(sample.files))
logger.info("Running over files %s", sample.files)

# Put together skim
isDiMuon        = options.skim.lower().startswith('dimuon')
isDiElectron    = options.skim.lower().startswith('dielectron')
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')

skimConds = []
if isDiMuon:
    skimConds.append( "Sum$(Muon_pt>15&&abs(Muon_eta)<2.5)>=2" )
elif isDiElectron:
    skimConds.append( "Sum$(Electron_pt>15&&abs(Electron_eta)<2.5)>=2" )
elif isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=2" )
elif isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

if options.skim.lower().count('met'):
    skimConds.append( "MET_pt>100" )

cut = '&&'.join(skimConds)

logger.info("Using selection: %s", cut)

# Main part

directory = "/afs/hephy.at/data/dspitzbart03/nanoSamples/%s_%s/"%(options.year, options.era)
output_directory = os.path.join( directory, options.skim, sample.name )

logger.info("Loading modules.")

if year == 2016:
    puwProducer = puWeightProducer(pufile_mc,pufile_data,"pu_mc","pileup",verbose=False)
    metSigParamsMC      = [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]
    metSigParamsData    = [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
    JER                 = "Summer16_25nsV1_MC"          if not sample.isData else "Summer16_25nsV1_DATA"
    JERera              = "Fall17_V3"
    if sample.isData:
        if sample.name.count("Run2016B") or sample.name.count("Run2016C") or sample.name.count("Run2016D"):
            JEC         = "Summer16_07Aug2017BCD_V11_DATA"
        elif sample.name.count("Run2016E") or sample.name.count("Run2016F"):
            JEC         = "Summer16_07Aug2017EF_V11_DATA"
        elif sample.name.count("Run2016G") or sample.name.count("Run2016H"):
            JEC         = "Summer16_07Aug2017GH_V11_DATA"
        else:
            raise NotImplementedError ("Don't know what JECs to use for sample %s"%sample.name)
    else:
        JEC             = "Summer16_07Aug2017_V11_MC"
    jetThreshold = 15

elif year == 2017:
    puwProducer = puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)
    # tuning from April 4th, with sumPt threshold of 25
    metSigParamsMC      = [1.7760438537732681, 1.720421230892687, 1.6034765551361112, 1.5336832981702226, 2.0928447254019757, 0.0011228025809342157, 0.7287313412909979]
    metSigParamsData    = [1.518621014453362, 1.611898248687222, 1.5136936762143423, 1.4878342676980971, 1.9192499533282406, -0.0005835026352392627, 0.749718704693196]
    JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
    JERera              = "Fall17_V3"
    if sample.isData:
        if sample.name.count('Run2017B'):
            JEC         = "Fall17_17Nov2017B_V32_DATA"
        elif sample.name.count('Run2017C'):
            JEC         = "Fall17_17Nov2017C_V32_DATA"
        elif sample.name.count('Run2017D') or sample.name.count('Run2017E'):
            JEC         = "Fall17_17Nov2017DE_V32_DATA"
        #elif sample.name.count('Run2017E'):
        #    JEC         = "Fall17_17Nov2017E_V32_DATA"
        elif sample.name.count('Run2017F'):
            JEC         = "Fall17_17Nov2017F_V32_DATA"
        else:
            raise NotImplementedError ("Don't know what JECs to use for sample %s"%sample.name)
    else:
        JEC             = "Fall17_17Nov2017_V32_MC"
    jetThreshold = 25

elif year == 2018:
    puwProducer = puWeightProducer("auto",pufile_data2018,"pu_mc","pileup",verbose=False)
    # tuning from April 4th, with sumPt threshold of 25
    metSigParamsMC      = [1.8430848616315363, 1.8572853766660877, 1.613083160233781, 1.3966398718198898, 1.4831008506492056, 0.0011310724285762122, 0.6929410058142578]
    metSigParamsData    = [1.6231076732985186, 1.615595174619551, 1.4731794897915416, 1.5183631493937553, 2.145670387603659, -0.0001524158603362826, 0.7510574688006575]
    JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
    JERera              = "Fall17_V3"
    if sample.isData:
        if sample.name.count("Run2018"):
            JEC         = "Autumn18_Run%s_V8_DATA"%era
        else:
            raise NotImplementedError ("Don't know what JECs to use for sample %s"%sample.name)
    else:
        JEC             = "Autumn18_V8_MC"
    jetThreshold = 25


unclEnThreshold = 15
metSigParams    = metSigParamsMC                if not sample.isData else metSigParamsData
METCollection   = "METFixEE2017" if year == 2017 else "MET"

jetmetProducer = jetmetUncertaintiesProducer(str(year), JEC, [ "Total" ], jer=JERera, jetType = "AK4PFchs", redoJEC=True, unclEnThreshold=unclEnThreshold, METBranchName=METCollection)
#jetmetProducer = jetmetUncertaintiesProducer(str(year), JEC, [ "Total" ], jer=JERera, jetType = "AK4PFchs", redoJEC=True, unclEnThreshold=unclEnThreshold)


if sample.isData:
    modules = [
        jetRecalib(JEC, unclEnThreshold=unclEnThreshold, METBranchName=METCollection),
        METSigTools(),
        lumiWeightProducer(1, isData=True),
        METSigProducer(JER, metSigParams, useRecorr=True, jetThreshold=jetThreshold, METCollection=METCollection),
        applyJSON(json),
        #METminProducer(isData=True),
        # MET significance producer
    ]

else:
    modules = [
        puwProducer,
        lumiWeightProducer(lumiScaleFactor),
        jetmetProducer,
        METSigTools(),
        METSigProducer(JER, metSigParams, useRecorr=True, calcVariations=True, jetThreshold=jetThreshold, METCollection=METCollection),
        applyJSON(None),
        #METminProducer(isData=False, calcVariations=True),
        # MET significance producer
    ]

logger.info("Preparing post-processor.")

p = PostProcessor(output_directory,sample.files,cut=cut, modules=modules)
if not options.prepare:
    logger.info("Running...")
    p.run()


