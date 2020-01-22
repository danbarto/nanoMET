# MetFilter Analysis Recommendations according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# Flag_BadChargedCandidateFilter is not recommended anymore, under review
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2

def getFilterCut( year, isData=False, isFastSim=False, skipBadChargedCandidate=True, skipBadPFMuon=False ):
    if year == 2016:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)

    elif year == 2017:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
#        filters            += [ "Flag_ecalBadCalibReducedMINIAODFilter" ]     # ECAL bad calibration filter update (needs to be re-run on miniAOD)
        filters            += [ "Flag_ecalBadCalibFilter" ]                   # current replacement for Flag_ecalBadCalibReducedMINIAODFilter -> change to Flag_ecalBadCalibFilterv2
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)
        if isFastSim:
            filters        += ["Flag_ecalBadCalibFilter"]
        else:
            filters        += ["Flag_ecalBadCalibFilterV2"]

    elif year == 2018:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
#        filters            += [ "Flag_ecalBadCalibReducedMINIAODFilter" ]     # ECAL bad calibration filter update (needs to be re-run on miniAOD)
        filters            += [ "Flag_ecalBadCalibFilter" ]                   # current replacement for Flag_ecalBadCalibReducedMINIAODFilter -> change to Flag_ecalBadCalibFilterv2
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)
        if isFastSim:
            filters        += ["Flag_ecalBadCalibFilter"]
        else:
            filters        += ["Flag_ecalBadCalibFilterV2"]

    else:
        raise NotImplementedError( "No MET filter found for year %i" %year )

    return "&&".join(filters)

