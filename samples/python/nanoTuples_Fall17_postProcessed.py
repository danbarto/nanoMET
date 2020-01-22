import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT
from nanoMET.samples.color import color

# Logging
import logging
logger = logging.getLogger(__name__)

# Data directory
try:
    data_directory = sys.modules['__main__'].data_directory
except:
    from nanoMET.tools.user import data_directory as user_data_directory
    data_directory = user_data_directory 

# Take post processing directory if defined in main module
try:
  import sys
  postProcessing_directory = sys.modules['__main__'].postProcessing_directory
except:
  postProcessing_directory = "2017_v22/dimuon/"

logger.info("Loading MC samples from directory %s", os.path.join(data_directory, postProcessing_directory))

dirs = {}

dirs['DY_LO']           = ["DYJetsToLL_M50_LO"]
dirs['Top']             = ["TTLep_pow", "T_tWch_ext", "TBar_tWch_ext"]
dirs['diboson']         = ['VVTo2L2Nu', 'WZTo3LNu_amcatnlo', 'WZTo1L3Nu']
dirs['triboson']        = ['WWZ_4F', 'WZZ', 'ZZZ']

directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

DY_LO_17            = Sample.fromDirectory(name="DY_LO",            treeName="Events", isData=False, color=color.DY,            texName="DY (LO)",                  directory=directories['DY_LO'])
Top_17              = Sample.fromDirectory(name="Top",              treeName="Events", isData=False, color=color.TTJets,        texName="t(#bar{t})",               directory=directories['Top'])
diboson_17          = Sample.fromDirectory(name="diboson",          treeName="Events", isData=False, color=color.diboson,       texName="diboson",                  directory=directories['diboson'])
rare_17             = Sample.fromDirectory(name="rare",             treeName="Events", isData=False, color=color.rare,          texName="rare",                     directory=directories['triboson'])
