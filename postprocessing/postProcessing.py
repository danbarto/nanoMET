'''
Post-processing based on nanoAOD tools
'''

import os
import ROOT

# Import tools
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

# Import modules
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer, pufile_data2016, pufile_mc2016, pufile_data2017, pufile_data2018, pufile_mc2017, pufile_mc2018
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer      import METSigProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2    import *

# private modules
from METSigTools           import METSigTools
from lumiWeightProducer    import lumiWeightProducer
from nanoMET.tools.user    import postprocessing_output_directory

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel',    action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--samples',     action='store', nargs='*', type=str, default=['TTZToLLNuNu_ext'], help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--skim',        action='store', nargs='?', type=str, default='dimuon', help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--job',         action='store', type=int, default=0, help="Run only jobs i" )
argParser.add_argument('--nJobs',       action='store', nargs='?', type=int,default=1, help="Maximum number of simultaneous jobs.")
argParser.add_argument('--prepare',     action='store_true', help="Prepare, don't acutally run" )
argParser.add_argument('--overwrite',   action='store_true', help="Overwrite" )
argParser.add_argument('--year',        action='store', default=None, help="Which year? Important for json file.")
argParser.add_argument('--era',         action='store', default="v1", help="Which era/subdirectory?")
options = argParser.parse_args()

# Logger
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   options.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None)

# from RootTools
from RootTools.core.standard            import *

def nonEmptyFile(f, treeName='Events'):
    logger.info("Checking file: %s"%f)
    rf = ROOT.TFile.Open(f)
    if not rf: return False
    tree = getattr(rf, treeName)
    nonEmpty = True if tree.GetEntries() else False
    if not nonEmpty: logger.info("File is empty")
    rf.Close()
    return nonEmpty

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

logger.info("Searching for sample %s"%options.samples[0])

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
logger.info("######### Era %s ########"%era)

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

directory = os.path.join( postprocessing_output_directory, "%s_%s"%(options.year, options.era) )
output_directory = os.path.join( directory, options.skim, sample.name )

fileNames = [ ('/'.join(x.split('/')[:-1]), x.split('/')[-1]) for x in sample.files if nonEmptyFile(x)  ]

if not options.overwrite:
    allFiles = []
    for f in fileNames:
        outfile = f[1].replace('.root','_Skim.root')
        if os.path.isfile(output_directory+'/'+outfile):
            try:
                a = helpers.checkRootFile(output_directory+'/'+outfile)
                if a:
                    logger.info('Found file %s, skipping.', outfile)
                    continue
                else:
                    logger.info("Found file %s, but need to rerun.", outfile)
                    allFiles.append(f)
            except:
                logger.info("Found file %s, but need to rerun.", outfile)
                allFiles.append(f)
        else:
            logger.info('No files found for %s, will produce them', outfile)
            allFiles.append(f)
else:
    allFiles = fileNames
        


sample.files = [ f[0] + '/' + f[1] for f in allFiles ]

logger.info("Loading modules.")

if year == 2016:
    ## tuning from November 2019, with sumPt threshold of 15
    puwProducer = puWeightProducer(pufile_mc2016,pufile_data2016,"pu_mc","pileup",verbose=False)
    metSigParamsMC      = [1.6789559564013943, 1.543666136735388, 1.4728342034302846, 1.4983602533711493, 1.4758351625239376, 0.008039429222660197, 0.6698834337575063]
    metSigParamsData    = [1.9034557745999647, 1.704569089762286, 1.5854229036413823, 1.4974876665993915, 1.673074548622476, 0.0015993706020479338, 0.6288393591242573]
    JER                 = "Summer16_25nsV1_MC"          if not sample.isData else "Summer16_25nsV1_DATA"
    archive             = '' if not sample.isData else "Summer16_07Aug2017_V11_DATA"
    jetThreshold = 15

elif year == 2017:
    puwProducer = puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)
    ## tuning from November 2019, with sumPt threshold of 15
    metSigParamsMC      = [1.7037614210331564, 1.7166071080686363, 1.6701114323915047, 1.502876236941622, 1.5780611987345947, -0.00012634174329968426, 0.6834329126092852]
    metSigParamsData    = [1.9410724258735805, 1.878895369894863, 1.9122297708165825, 1.7090755971750793, 2.0004413703111146, -0.0001347239857148459, 0.6732736907156339]
    JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
    archive             = '' if not sample.isData else  "Fall17_17Nov2017_V32_DATA"
    jetThreshold = 15

elif year == 2018:
    puwProducer = puWeightProducer("auto",pufile_data2018,"pu_mc","pileup",verbose=False)
    ## tuning from November 2019, with sumPt threshold of 15
    metSigParamsMC      = [1.8549938037723896, 1.6509411315853, 1.636437587406195, 1.4672590328655033, 1.247095140902097, -6.784962823364049e-05, 0.6354413715418223]
    metSigParamsData    = [1.7597963093944353, 1.7529358941285436, 1.6542030082916106, 1.355947946214444, 1.62529299229821, 0.0003583146878367062, 0.6808117480506645]
    JER                 = "Autumn18_V7b_MC" if not sample.isData else "Autumn18_V7b_DATA"
    archive             = '' if not sample.isData else "Autumn18_V19_DATA"
    jetThreshold = 15

if len(metSigParamsMC) == 12 and len(metSigParamsData) == 12:
    pTdependent = True
elif len(metSigParamsMC) == 7 and len(metSigParamsData) == 7:
    pTdependent = False
else:
    raise Exception("Wrong input parameter settings!" )



unclEnThreshold = 15
metSigParams    = metSigParamsMC if not sample.isData else metSigParamsData
METCollection   = "METFixEE2017" if year == 2017 else "MET"

vetoEtaRegion = (2.65, 3.14) if year == 2017 else (10,10)

METBranchName = 'MET' if not year == 2017 else 'METFixEE2017'
JMECorrector = createJMECorrector(isMC=(not sample.isData), dataYear=year, runPeriod=era, jesUncert="Total", jetType = "AK4PFchs", metBranchName=METBranchName, isFastSim=False, applySmearing=False)
modules = [
    JMECorrector()
]


if sample.isData:
    modules += [
        METSigTools(),
        lumiWeightProducer(1, isData=True),
        METSigProducer(JER, metSigParams, useRecorr=True, jetThreshold=jetThreshold, METCollection=METCollection, vetoEtaRegion=vetoEtaRegion, pTdependent=pTdependent),
    ]

else:
    modules += [
        puwProducer,
        lumiWeightProducer(lumiScaleFactor),
        METSigTools(),
        METSigProducer(JER, metSigParams, useRecorr=True, calcVariations=True, jetThreshold=jetThreshold, METCollection=METCollection, vetoEtaRegion=vetoEtaRegion, pTdependent=pTdependent),
    ]

logger.info("Preparing post-processor.")

p = PostProcessor(output_directory,sample.files,cut=cut, modules=modules, jsonInput=(json if sample.isData else None) ) # could make use of hadd to reduce number of files in the future.
if not options.prepare:
    logger.info("Running ... ")
    p.run()

