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
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer       import puWeightProducer, pufile_data2016, pufile_mc2016, pufile_data2017, pufile_data2018, pufile_mc2017, pufile_mc2018
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
argParser.add_argument('--overwrite',   action='store_true', help="Overwrite" )
argParser.add_argument('--year',        action='store', default=None, help="Which year? Important for json file.")
argParser.add_argument('--era',         action='store', default="v1", help="Which era/subdirectory?")
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)

# from RootTools
from RootTools.core.standard            import *

def nonEmptyFile(f, treeName='Events'):
    print "Checking file: %s"%f
    rf = ROOT.TFile.Open(f)
    if not rf: return False
    tree = getattr(rf, treeName)
    nonEmpty = True if tree.GetEntries() else False
    if not nonEmpty: print "File is empty"
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

#logger.info("Sample contains %s files", len(sample.files))
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
#logger.info("Sample has %s files", len_orig)
json = sample.json # pickup json from sample, as defined in the Sample repository
sample = sample.split( n=options.nJobs, nSub=options.job)
sample.name = keepSampleName

logger.info("Will run over %s files", len(sample.files))
#logger.info("Running over files %s", sample.files)

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
    puwProducer = puWeightProducer(pufile_mc2016,pufile_data2016,"pu_mc","pileup",verbose=False)
    metSigParamsMC      = [1.617529475909303, 1.617529475909303, 1.4505983036866312, 1.4505983036866312, 1.411498565372343, 1.411498565372343, 1.4087559908291813, 1.4087559908291813, 1.3633674107893856, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]
    metSigParamsData    = [1.843242937068234, 1.843242937068234, 1.64107911184195,   1.64107911184195,   1.567040591823117, 1.567040591823117, 1.5077143780804294, 1.5077143780804294, 1.614014783345394,  1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
    JER                 = "Summer16_25nsV1_MC"          if not sample.isData else "Summer16_25nsV1_DATA"
    JERera              = "Summer16_25nsV1"
    archive             = ''
    if sample.isData:
        archive = "Summer16_07Aug2017_V11_DATA"
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
    # tuning from June 24th, with sumPt threshold of 15
    metSigParamsMC      = [1.9648214119268503, 1.5343086462230238, 1.9167197601498538, 1.5145044341064964, 1.8069380221985405, 1.3217263662622654, 1.5506294867561126, 1.272977540964842,  1.50742322311234,   1.6542883449796797, -0.0017865650107230548,  0.6593106706741719]
    metSigParamsData    = [2.228118299837604,  1.2420725475347338, 2.227630982417529,  1.256752205787215,  2.0215250734187853, 1.1557507029911258, 1.7350536144535336, 1.1587692458345757, 1.9385081854607988, 1.8726188460472792, -2.6697894266706265e-05, 0.646984812801919]
    JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
    JERera              = "Fall17_V3"
    archive             = ''
    if sample.isData:
        archive = "Fall17_17Nov2017_V32_DATA"
        if sample.name.count('Run2017B'):
            JEC         = "Fall17_17Nov2017B_V32_DATA"
        elif sample.name.count('Run2017C'):
            JEC         = "Fall17_17Nov2017C_V32_DATA"
        elif sample.name.count('Run2017D') or sample.name.count('Run2017E'):
            JEC         = "Fall17_17Nov2017DE_V32_DATA"
        elif sample.name.count('Run2017F'):
            JEC         = "Fall17_17Nov2017F_V32_DATA"
        else:
            raise NotImplementedError ("Don't know what JECs to use for sample %s"%sample.name)
    else:
        JEC             = "Fall17_17Nov2017_V32_MC"
    jetThreshold = 15

elif year == 2018:
    puwProducer = puWeightProducer("auto",pufile_data2018,"pu_mc","pileup",verbose=False)
    ## tuning from October 2019, with sumPt threshold of 25
    metSigParamsMC      = [1.9033100259447273, 1.7374326087706358, 1.6957838005481387, 1.7283187962364615, 1.5868244614361302, 1.526252837049335, 1.3744055574137417, 1.4500298644941831, 1.4796204632654997, 1.4481227819959115, 0.0019899110367503207, 0.6927496536100137]
    metSigParamsData    = [1.7492714572981323, 1.3430198915313956, 1.911348554103867, 1.3718438058490257, 1.5885672661901442, 1.4385903478138795, 1.521901070409261, 1.4522895772008289, 1.8870084799263003, 1.7138750357657668, -1.359837542505224e-06, 0.7434240965656385]
    JER                 = "Autumn18_V1_MC"                if not sample.isData else "Autumn18_V1_DATA"
    JERera              = "Autumn18_V1"
    archive             = ''
    if sample.isData:
        archive = "Autumn18_V19_DATA"
        if sample.name.count("Run2018"):
            JEC         = "Autumn18_Run%s_V19_DATA"%era
        else:
            raise NotImplementedError ("Don't know what JECs to use for sample %s"%sample.name)
    else:
        JEC             = "Autumn18_V19_MC"
    jetThreshold = 25


unclEnThreshold = 15
metSigParams    = metSigParamsMC                if not sample.isData else metSigParamsData
METCollection   = "METFixEE2017" if year == 2017 else "MET"

vetoEtaRegion = (2.65, 3.14) if year == 2017 else (10,10)

print vetoEtaRegion




if sample.isData:
    modules = [
        #jetRecalib(JEC, JEC, METBranchName=METCollection),
        jetmetUncertaintiesProducer(str(year), JEC, archive=archive, jesUncertainties=[ "Total" ], jetType = "AK4PFchs", redoJEC=True, metBranchName=METCollection, isData=True),
        METSigTools(),
        lumiWeightProducer(1, isData=True),
        METSigProducer(JER, metSigParams, useRecorr=True, jetThreshold=jetThreshold, METCollection=METCollection, vetoEtaRegion=vetoEtaRegion),
        #applyJSON(json),
        #METminProducer(isData=True),
    ]

else:
    modules = [
        puwProducer,
        lumiWeightProducer(lumiScaleFactor),
        jetmetUncertaintiesProducer(str(year), JEC, jesUncertainties=[ "Total" ], jetType = "AK4PFchs", redoJEC=True, metBranchName=METCollection),
        #jetmetProducer,
        METSigTools(),
        METSigProducer(JER, metSigParams, useRecorr=True, calcVariations=True, jetThreshold=jetThreshold, METCollection=METCollection, vetoEtaRegion=vetoEtaRegion),
        #applyJSON(None),
        #METminProducer(isData=False, calcVariations=True),
    ]

logger.info("Preparing post-processor.")

p = PostProcessor(output_directory,sample.files,cut=cut, modules=modules, jsonInput=(json if sample.isData else None) ) # could make use of hadd to reduce number of files in the future.
if not options.prepare:
    logger.info("Running ... ")
    p.run()

