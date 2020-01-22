import os

from tunes import tuneParams
tuneEras = tuneParams.keys()
tuneEras.sort()
commands = []

for tune in tuneEras:
    print tune
    if "2017" in str(tune):
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ-BadEEJetVeto             --plot_directory ppv22_v1"%tune )
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ-nCleanJet0-BadEEJetVeto  --plot_directory ppv22_v1"%tune )
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ-nCleanJet1p-BadEEJetVeto --plot_directory ppv22_v1"%tune )
    else:
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ             --plot_directory ppv22_v1"%tune )
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ-nCleanJet0  --plot_directory ppv22_v1"%tune )
        commands.append( "python analysisPlots.py --calcMETSig --tuneEra %s --selection diMuonTune-tuneElectronVeto-onZ-nCleanJet1p --plot_directory ppv22_v1"%tune )


with open( "runAllPlots.sh", "w" ) as f:
    for com in commands:
        f.write(com+"\n")
