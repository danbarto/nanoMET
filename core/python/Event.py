'''

'''

import math

from nanoMET.core.METSignificance import METSignificance

MS = METSignificance()

etabins = [0.8,1.3,1.9,2.5,100]

def getBin(abseta):
    for i, a in enumerate(etabins):
        if abseta < a:
            return int(i)
            break

def cartesian(pt, phi):
    return (pt*math.cos(phi), pt*math.sin(phi))

class Event:
    def __init__(self, event):
        #self.MS = METSignificance()
        MS.getJER(event)

        self.nJet       = event.nJet
        self.Jet_pt     = [ x for x in event.Jet_pt     if x > 0 ]
        self.Jet_eta    = [ x for x in event.Jet_eta    if not math.isnan(x) ]
        self.Jet_etabin = [ getBin(abs(x)) for x in self.Jet_eta ]
        self.Jet_phi    = [ x for x in event.Jet_phi    if not math.isnan(x) ]
        self.Jet_dpt    = [ x for x in event.Jet_dpt ]
        self.Jet_dphi   = [ x for x in event.Jet_dphi ]

        self.MET_pt             = event.MET_pt
        self.MET_phi            = event.MET_phi
        self.MET_sumEt          = event.MET_sumEt
        self.MET_significance   = event.MET_significance
        


        self.fixedGridRhoFastjetAll = event.fixedGridRhoFastjetAll
        self.weight = event.weight

    def METSignificance(self, args):

        MS.calculate(self, args)

        try:
            LL = self.weight * ( self.MET_sig + math.log(self.det) )
        except:
            LL = self.weight * self.MET_sig

        return LL


