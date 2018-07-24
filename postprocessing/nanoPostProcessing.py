#!/usr/bin/env python

'''
Slimmed down post-processing taken from several physics analysis.
Doing several simple tasks that simplfy life afterwards, e.g.
- Calculate PU and lumi weights for MC.
- Apply golden json to data.
- Skim and remove collections that are not needed.
'''

# standard imports
import ROOT
import sys
import os
import copy
import random
import subprocess
import datetime
import shutil

from array import array
from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import nanoMET.tools.user as user

from nanoMET.tools.helpers import closestOSDLMassToMZ, checkRootFile, writeObjToFile, m3, deltaR, bestDRMatchInCollection

targetLumi = 1000 #pb-1 Which lumi to normalize to

logChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',                    action='store',         nargs='?',              choices=logChoices,     default='INFO',                     help="Log level for logging")
    argParser.add_argument('--overwrite',                   action='store_true',                                                                                        help="Overwrite existing output files, bool flag set to True  if used")
    argParser.add_argument('--samples',                     action='store',         nargs='*',  type=str,                           default=['WZTo3LNu'],               help="List of samples to be post-processed, given as CMG component name")
    argParser.add_argument('--eventsPerJob',                action='store',         nargs='?',  type=int,                           default=30000000,                   help="Maximum number of events per job (Approximate!).") # mul by 100
    argParser.add_argument('--nJobs',                       action='store',         nargs='?',  type=int,                           default=1,                          help="Maximum number of simultaneous jobs.")
    argParser.add_argument('--job',                         action='store',                     type=int,                           default=0,                          help="Run only job i")
    argParser.add_argument('--minNJobs',                    action='store',         nargs='?',  type=int,                           default=1,                          help="Minimum number of simultaneous jobs.")
    argParser.add_argument('--fileBasedSplitting',          action='store_true',                                                                                        help="Split njobs according to files")
    argParser.add_argument('--dataDir',                     action='store',         nargs='?',  type=str,                           default="/a/b/c",                   help="Name of the directory where the input data is stored (for samples read from Heppy).")
    argParser.add_argument('--targetDir',                   action='store',         nargs='?',  type=str,                           default=user.postprocessing_output_directory, help="Name of the directory the post-processed files will be saved")
    argParser.add_argument('--processingEra',               action='store',         nargs='?',  type=str,                           default='v1',                       help="Name of the processing era")
    argParser.add_argument('--skim',                        action='store',         nargs='?',  type=str,                           default='dilepTiny',                help="Skim conditions to be applied for post-processing")
    argParser.add_argument('--small',                       action='store_true',                                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used")
    argParser.add_argument('--year',                        action='store',                     type=int,   choices=[2016,2017],    required = True,                    help="Which year?")

    return argParser

options = get_parser().parse_args()

# Logging
import nanoMET.tools.logger as logger
logFile = '/tmp/%s_%s_%s_njob%s.txt'%(options.skim, '_'.join(options.samples), os.environ['USER'], str(0 if options.nJobs==1 else options.job) )
logger  = logger.get_logger(options.logLevel, logFile = logFile)

import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None )


# Flags 
isSingleLep =   options.skim.lower().startswith('singlelep')
isDiLep     =   options.skim.lower().startswith('dilep')
isSmall     =   options.skim.lower().count('small')

# Skim condition
skimConds = []
if isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )
elif isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=2" )

#Samples: Load samples
maxN = 2 if options.small else None

# import data and MC
#from nanoMET.samples.helpers import fromNanoSample
if options.year == 2016:
    from nanoMET.samples.nanoTuples_Summer16 import *
elif options.year == 2017:
    from nanoMET.samples.nanoTuples_Fall17 import *
else:
    raise NotImplementedError

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
            samples.append(sample)

#samples = [ fromNanoSample(s, maxN = maxN) for s in options.samples ]

isData = False not in [s.isData for s in samples]
isMC = True not in [s.isData for s in samples]

# Check that all samples which are concatenated have the same x-section.
assert isData or len(set([s.xSection for s in samples]))==1, "Not all samples have the same xSection: %s !"%(",".join([s.name for s in samples]))
assert isMC or len(samples)==1, "Don't concatenate data samples"

xSection = samples[0].xSection if isMC else None

if options.fileBasedSplitting:
    len_orig = len(sample.files)
    sample = sample.split( n=options.nJobs, nSub=options.job)
    logger.info( "fileBasedSplitting: Run over %i/%i files for job %i/%i."%(len(sample.files), len_orig, options.job, options.nJobs))
    logger.debug( "fileBasedSplitting: Files to be run over:\n%s", "\n".join(sample.files) )

if isMC:
    from nanoMET.tools.puReweighting import getReweightingFunction
    mcProfile = "Summer16"
    # nTrueIntReweighting
    nTrueInt_puRW        = getReweightingFunction(data="PU_2016_36000_XSecCentral", mc=mcProfile)
    nTrueInt_puRWDown    = getReweightingFunction(data="PU_2016_36000_XSecDown",    mc=mcProfile)
    nTrueInt_puRWUp      = getReweightingFunction(data="PU_2016_36000_XSecUp",      mc=mcProfile)

## top pt reweighting to be added later
#

# file and directory handling
directory = os.path.join(options.targetDir, options.processingEra) 
output_directory = os.path.join( directory, options.skim, sample.name )

if os.path.exists(output_directory) and options.overwrite:
    if options.nJobs > 1:
        logger.warning( "NOT removing directory %s because nJobs = %i", output_directory, options.nJobs )
    else:
        logger.info( "Output directory %s exists. Deleting.", output_directory )
        shutil.rmtree(output_directory)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(output_directory)
    logger.info( "Created output directory %s.", output_directory )
except:
    pass


### branch handling
#branches to be kept for data and MC
branchKeepStrings_DATAMC = [\
    "run", "luminosityBlock", "event", "fixedGridRhoFastjetAll", "PV_npvs", "PV_npvsGood",
    "MET_pt", "MET_phi","MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY", "MET_sumEt", "CaloMET_phi", "CaloMET_pt", "CaloMET_sumEt", "MET_covXX", "MET_covXY", "MET_covYY", "MET_significance",
    "RawMET_phi", "RawMET_pt", "RawMET_sumEt",
    "Flag_*","HLT_*",
    "nJet", "Jet_*",
    "nElectron", "Electron_*",
    "nPhoton", "Photon_*",
    "nMuon", "Muon_*",
    "nTau", "Tau_*",
]

#branches to be kept for MC samples only
branchKeepStrings_MC = [\
    #"nTrueInt", 
    "genWeight", "LHE_HTIncoming",
    "LHEWeight_originalXWGTUP","nLHEPdfWeight","LHEPdfWeight","nLHEScaleWeight","LHEScaleWeight",
    "nGenPart","GenPart_*",
    "GenMET_pt", "GenMET_phi",
]

#branches to be kept for data only
branchKeepStrings_DATA = [ ]


if sample.isData:
    lumiScaleFactor=None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA
    from FWCore.PythonUtilities.LumiList import LumiList
    # Apply golden JSON
    if options.year == 2016:
        sample.json = '$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    elif options.year == 2017:
        sample.json = '$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    else:
        sample.json = '$CMSSW_BASE/src/nanoMET/tools/data/json/Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt'
    
    lumiList = LumiList(os.path.expandvars(sample.json))
    logger.info( "Loaded json %s", sample.json )
else:
    lumiScaleFactor = xSection*targetLumi/float(sample.normalization) if xSection is not None else None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC

#print lumiScaleFactor, sample.normalization

read_variables = map(TreeVariable.fromString, ['MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] )
new_variables = [ 'weight/F']
if isMC:
    read_variables+= map(TreeVariable.fromString, ['genWeight/F', 'PV_npvsGood/I'] )
if isData: new_variables.extend( ['jsonPassed/I','isData/I'] )


logger.info("Defining reader")

# Define a reader
reader = sample.treeReader( \
    variables = read_variables ,
    selectionString = "&&".join(skimConds)
    )

def filler( event ):
    # shortcut
    r = reader.event
    event.isData = s.isData

    if isMC:
        event.weight = lumiScaleFactor*r.genWeight if lumiScaleFactor is not None else 1
    elif isData:
        event.weight = 1

    if isData:
        event.jsonPassed  = lumiList.contains(r.run, r.luminosityBlock)

    if isMC:
        event.reweightPU        = nTrueInt_puRW     ( r.PV_npvsGood ) 
        event.reweightPUDown    = nTrueInt_puRWDown ( r.PV_npvsGood ) 
        event.reweightPUUp      = nTrueInt_puRWUp   ( r.PV_npvsGood ) 


# Create a maker. Maker class will be compiled. This instance will be used as a parent in the loop
logger.info("Compling maker class")
treeMaker_parent = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in new_variables ],
    treeName = "Events"
)


# Split input in ranges
if options.nJobs>1 and not options.fileBasedSplitting:
    eventRanges = reader.getEventRanges( nJobs = options.nJobs )
else:
    eventRanges = reader.getEventRanges( maxNEvents = options.eventsPerJob, minJobs = options.minNJobs )

logger.info( "Splitting into %i ranges of %i events on average. FileBasedSplitting: %s. Job number %s",  
        len(eventRanges), 
        (eventRanges[-1][1] - eventRanges[0][0])/len(eventRanges), 
        'Yes' if options.fileBasedSplitting else 'No',
        options.job)

#Define all jobs
jobs = [(i, eventRanges[i]) for i in range(len(eventRanges))]

filename, ext = os.path.splitext( os.path.join(output_directory, sample.name + '.root') )

if options.fileBasedSplitting and len(eventRanges)>1:
    raise RuntimeError("Using fileBasedSplitting but have more than one event range!")

clonedEvents = 0
convertedEvents = 0
outputLumiList = {}

for ievtRange, eventRange in enumerate( eventRanges ):

    
    if not options.fileBasedSplitting and options.nJobs>1:
        if ievtRange != options.job: continue

    logger.info( "Processing range %i/%i from %i to %i which are %i events.",  ievtRange, len(eventRanges), eventRange[0], eventRange[1], eventRange[1]-eventRange[0] )

    # Check whether file exists
    fileNumber = options.job if options.job is not None else 0
    outfilename = filename+'_'+str(fileNumber)+ext
    if os.path.isfile(outfilename):
        logger.info( "Output file %s found.", outfilename)
        if not checkRootFile(outfilename, checkForObjects=["Events"]):
            logger.info( "File %s is broken. Overwriting.", outfilename)
        elif not options.overwrite:
            logger.info( "Skipping.")
            continue
        else:
            logger.info( "Overwriting.")

    tmp_directory = ROOT.gDirectory
    outputfile = ROOT.TFile.Open(outfilename, 'recreate')
    tmp_directory.cd()

    if options.small: 
        logger.info("Running 'small'. Not more than 10000 events") 
        nMaxEvents = eventRange[1]-eventRange[0]
        eventRange = ( eventRange[0], eventRange[0] +  min( [nMaxEvents, 10000] ) )

    # Set the reader to the event range
    reader.setEventRange( eventRange )

    clonedTree = reader.cloneTree( branchKeepStrings, newTreename = "Events", rootfile = outputfile )
    clonedEvents += clonedTree.GetEntries()
    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()

    while reader.run():
        maker.run()
        if isData:
            if maker.event.jsonPassed_:
                if reader.event.run not in outputLumiList.keys():
                    outputLumiList[reader.event.run] = set([reader.event.luminosityBlock])
                else:
                    if reader.event.luminosityBlock not in outputLumiList[reader.event.run]:
                        outputLumiList[reader.event.run].add(reader.event.luminosityBlock)

    convertedEvents += maker.tree.GetEntries()
    maker.tree.Write()
    outputfile.Close()
    logger.info( "Written %s", outfilename)

  # Destroy the TTree
    maker.clear()


logger.info( "Converted %i events of %i, cloned %i",  convertedEvents, reader.nEvents , clonedEvents )

# Storing JSON file of processed events
if isData:
    jsonFile = filename+'_%s.json'%(0 if options.nJobs==1 else options.job)
    LumiList( runsAndLumis = outputLumiList ).writeJSON(jsonFile)
    logger.info( "Written JSON file %s", jsonFile )

logger.info("Copying log file to %s", output_directory )
copyLog = subprocess.call(['cp', logFile, output_directory] )
if copyLog:
    logger.info( "Copying log from %s to %s failed", logFile, output_directory)
else:
    logger.info( "Successfully copied log file" )
    os.remove(logFile)
    logger.info( "Removed temporary log file" )


