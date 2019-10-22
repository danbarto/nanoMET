# nanoMET
Repository to use nanoAOD tuples to tune MET Significance, and produce validation plots

Recipe:

# CMSSW
```
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
git cms-init
```

# get all the packages
```
cd $CMSSW_BASE/src
git clone https://github.com/HephyAnalysisSW/Samples.git
git clone https://github.com/HephyAnalysisSW/RootTools.git
git clone https://github.com/HephyAnalysisSW/Analysis.git
```

# nanoAOD-tools
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
git remote add metsig https://github.com/danbarto/nanoAOD-tools.git
git fetch metsig
git checkout METSigProducer
```

# nanoMET
```
cd $CMSSW_BASE/src
git clone -b moreParams https://github.com/danbarto/nanoMET.git
scram b -j 8
```

