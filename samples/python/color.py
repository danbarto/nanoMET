# Standard imports 
import os
import ROOT

# RootTools
from RootTools.core.standard import *

# Logging
import logging
logger = logging.getLogger(__name__)

from nanoMET.samples.helpers import singleton

@singleton
class color():
    pass

color.data          = ROOT.kBlack
color.DY            = ROOT.kGreen-6
color.TTJets        = ROOT.kCyan-6
color.diboson       = ROOT.kRed-7
color.TTZ           = ROOT.kMagenta-5
color.WJets         = ROOT.kBlue-6
color.rare          = ROOT.kBlue-2
