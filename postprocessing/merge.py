import os

# argparser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")
argParser.add_argument('--logLevel',    action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging" )
argParser.add_argument('--dir',         action='store', nargs='*', type=str, help="Directory with root files" )
argParser.add_argument('--nFiles',      action='store', nargs='*', type=int, default=100, help="How many root files to merge." )
options = argParser.parse_args()

import nanoMET.tools.logger as _logger
logger = _logger.get_logger(options.logLevel, logFile = None)

# from RootTools
from RootTools.core.standard            import *

import glob

def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in xrange(0, len(l), n))

filePath = str(options.dir[0])+'/*.root'

files = glob.glob(filePath)

for i,c in enumerate(chunks(files, options.nFiles)):
    
    outFile = "%s/merged_nanoAOD_%s.root"%(options.dir[0], i)
    
    nCorrupted = 0
    for f in c:
        if not helpers.checkRootFile(f):
            logger.info("File corrupted? %s", f)
            nCorrupted += 1

    logger.info("Number of corrupted files %s", nCorrupted)

    inFiles = ' '.join(c)
    
    print inFiles
    os.system("haddnano.py %s %s"%(outFile, inFiles))

