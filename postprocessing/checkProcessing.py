import os
# from RootTools
from RootTools.core.standard            import *

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel',    action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--samples',     action='store', nargs='*', type=str, default=['TTZToLLNuNu_ext'], help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--skim',        action='store', nargs='?', type=str, default='dimuon', help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--era',         action='store', default="2016_v1", help="Which era/subdirectory?")
#argParser.add_argument(
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)


#from Samples.nanoAOD.Summer16 import allSamples as Summer16
from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as Summer16_METSig_v1
from Samples.nanoAOD.Run2018_26Sep2018_private import allSamples as Run2018_26Sep2018_private
from Samples.nanoAOD.Run2016_17Jul2018_private import allSamples as Run2016_17Jul2018_private

allSamples = Summer16_METSig_v1 + Run2018_26Sep2018_private + Run2016_17Jul2018_private

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


directory = "/afs/hephy.at/data/dspitzbart03/nanoSamples/%s/"%options.era
output_directory = os.path.join( directory, options.skim, sample.name )

logger.info("Starting checks.")

import glob

files = glob.glob("%s/*.root"%output_directory)

if len(files) == len(sample.files):
    logger.info("Number of files is correct.")

nCorrupted = 0
for f in files:
    if not helpers.checkRootFile(f):
        logger.info("File corrupted? %s", f)
        nCorrupted += 1

if nCorrupted == 0:
    logger.info("All files seem fine.")

logger.info("Done.")

