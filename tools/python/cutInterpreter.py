''' Class to interpret string based cuts
'''

import logging
logger = logging.getLogger(__name__)

jetSelection    = "nJetGood"
bJetSelectionM  = "nBTag"

mIsoWP = { "VT":5, "T":4, "M":3 , "L":2 , "VL":1, 0:"None" }

special_cuts = {
    # ("multiIsoVT":        "(1)", 
    "diMuonTune":       "(Sum$(Muon_pt>25&&Muon_isGoodMuon)==2)",
    "diMuon":           "(Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0)",
    "diMuon16":         "(Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>25&&Muon_isGoodMuon)>0)",
    "diMuon1718":       "(Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>30&&Muon_isGoodMuon)>0)",
    "looseLeptonVeto":  "(Sum$(Electron_pt>15&&abs(Electron_eta)<2.4&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>15&&abs(Muon_eta)<2.4&&Muon_pfRelIso03_all<0.4) )==2",
    "tuneElectronVeto": "(Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0)",
    "lepSel":           "l1_pt>30 && l2_pt>20",
    "badJetSrEVeto":    "Sum$(Jet_neEmEF*Jet_pt*cosh(Jet_eta)*(2.5<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt<50))<200",
    "noEEJets":         "Sum$(Jet_pt>15&&abs(Jet_eta)>2.65&&abs(Jet_eta)<3.14)==0",
    "badJetSrEVetoV2":  "Sum$(Jet_neEmEF*Jet_pt*cosh(Jet_eta)*(2.5<abs(Jet_eta)&&abs(Jet_eta)<3))<(-(3./7)*Sum$(Jet_pt*cosh(Jet_eta)*(2.5<abs(Jet_eta)&&abs(Jet_eta)<3))+430)",
    "dPhiEEMet":        "(Sum$((-cos(Jet_phi-MET_phi)*Jet_pt*Jet_neEmEF)*(cos(Jet_phi-MET_phi)<cos(2*pi/3.))*(2.5<abs(Jet_eta)&&abs(Jet_eta)<3))/MET_pt)<0.2",
    "dPhiEEMetV2":      "(Sum$((-cos(Jet_phi-MET_phi)*Jet_pt*Jet_neEmEF*cosh(Jet_eta))*(cos(Jet_phi-MET_phi)<cos(2*pi/3.))*(2.5<abs(Jet_eta)&&abs(Jet_eta)<3))/MET_pt)<0.5",
    "incl":             "(1)",
    "allZ":             "(1)",
    "onZ":              "abs(dl_mass-91.1876)<10",
    "offZ":             "abs(dl_mass-91.1876)>10",
    "BadEEJetVeto":     "Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>30))==0",
    "BadEEJetVetoLow":  "Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>15))==0",
  }

continous_variables = [ ("met", "MET_pt"), ("mll", "dl_mass"), ("nPV", "PV_npvsGood"), ("nSoftJet", "Sum$(Jet_pt>15&&Jet_pt<30)"), ("sumPt", "MET_sumPt") ]
#discrete_variables  = [ ("njet", "Sum$(Jet_pt>30&&Jet_jetId>0&&Jet_cleanmask>0&&abs(Jet_eta)<2.4)"), ("btag", "Sum$(Jet_pt>30&&Jet_jetId>0&&Jet_cleanmask>0&&abs(Jet_eta)<2.4&&Jet_btagDeepB>0.4941)")]
discrete_variables  = [ ("nJet", "Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)"), ("nCleanJet", "Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4&&Jet_cleanmask)"), ("btagA", "Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4&&Jet_btagDeepB>0.4184&&Jet_cleanmask)"), ("nsoftJet", "Sum$(Jet_pt<30)")]

class cutInterpreter:
    ''' Translate var100to200-var2p etc.
    '''

    @staticmethod
    def translate_cut_to_string( string ):

        if string.startswith("multiIso"):
            str_ = mIsoWP[ string.replace('multiIso','') ]
            return "l1_mIsoWP>%i&&l2_mIsoWP>%i" % (str_, str_)
        elif string.startswith("relIso"):
           iso = float( string.replace('relIso','') )
           return "l1_relIso03<%3.2f&&l2_relIso03<%3.2f"%( iso, iso )
        # special cuts
        if string in special_cuts.keys(): return special_cuts[string]

        # continous Variables
        for var, tree_var in continous_variables:
            if string.startswith( var ):
                num_str = string[len( var ):].replace("to","To").split("To")
                upper = None
                lower = None
                if len(num_str)==2:
                    lower, upper = num_str
                elif len(num_str)==1:
                    lower = num_str[0]
                else:
                    raise ValueError( "Can't interpret string %s" % string )
                res_string = []
                if lower: res_string.append( tree_var+">="+lower )
                if upper: res_string.append( tree_var+"<"+upper )
                return "&&".join( res_string )

        # discrete Variables
        for var, tree_var in discrete_variables:
            logger.debug("Reading discrete cut %s as %s"%(var, tree_var))
            if string.startswith( var ):
                # So far no njet2To5
                if string[len( var ):].replace("to","To").count("To"):
                    raise NotImplementedError( "Can't interpret string with 'to' for discrete variable: %s. You just volunteered." % string )

                num_str = string[len( var ):]
                # logger.debug("Num string is %s"%(num_str))
                # var1p -> tree_var >= 1
                if num_str[-1] == 'p' and len(num_str)==2:
                    # logger.debug("Using cut string %s"%(tree_var+">="+num_str[0]))
                    return tree_var+">="+num_str[0]
                # var123->tree_var==1||tree_var==2||tree_var==3
                else:
                    vls = [ tree_var+"=="+c for c in num_str ]
                    if len(vls)==1:
                      # logger.debug("Using cut string %s"%vls[0])
                      return vls[0]
                    else:
                      # logger.debug("Using cut string %s"%'('+'||'.join(vls)+')')
                      return '('+'||'.join(vls)+')'
        raise ValueError( "Can't interpret string %s. All cuts %s" % (string,  ", ".join( [ c[0] for c in continous_variables + discrete_variables] +  special_cuts.keys() ) ) )

    @staticmethod
    def cutString( cut, select = [""], ignore = []):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        cuts = cut.split('-')
        # require selected
        cuts = filter( lambda c: any( sel in c for sel in select ), cuts )
        # ignore
        cuts = filter( lambda c: not any( ign in c for ign in ignore ), cuts )

        cutString = "&&".join( map( cutInterpreter.translate_cut_to_string, cuts ) )


        return cutString
    
    @staticmethod
    def cutList ( cut, select = [""], ignore = []):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        cuts = cut.split('-')
        # require selected
        cuts = filter( lambda c: any( sel in c for sel in select ), cuts )
        # ignore
        cuts = filter( lambda c: not any( ign in c for ign in ignore ), cuts )
        return [ cutInterpreter.translate_cut_to_string(cut) for cut in cuts ] 
        #return  "&&".join( map( cutInterpreter.translate_cut_to_string, cuts ) )

if __name__ == "__main__":
    print cutInterpreter.cutString("njet2p-onZ")
    print
    print cutInterpreter.cutList("njet2p-onZ")
