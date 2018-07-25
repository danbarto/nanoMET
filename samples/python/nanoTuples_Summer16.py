import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

dbFile = '/afs/hephy.at/data/dspitzbart01/nanoAOD/DB_Summer16.sql'
dataDirectory = '/afs/hephy.at/data/dspitzbart01/nanoAOD/Summer16/'

# specify a local directory if you want to create (and afterwards automatically use) a local copy of the sample, otherwise use the grid.

## DY
DYJetsToLL_M50_LO   = Sample.nanoAODfromDAS('DYJetsToLL_M50_LO','/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dbFile=dbFile, xSection=2008.*3)
DYJetsToLL_M50      = Sample.nanoAODfromDAS('DYJetsToLL_M50','/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM', dbFile=dbFile, xSection=2008.*3)



DY = [
    DYJetsToLL_M50_LO,
    DYJetsToLL_M50,
    ]

## top
TTLep_pow       = Sample.nanoAODfromDAS("TTLep_pow", "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM", dbFile=dbFile, xSection=831.76*((3*0.108)**2) )

top = [
    TTLep_pow,
    ]

## di/multiboson
WW      = Sample.nanoAODfromDAS("WW", "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM", dbFile=dbFile, xSection=63.21 * 1.82)
WW_ext  = Sample.nanoAODfromDAS("WW_ext", "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM", dbFile=dbFile, xSection=63.21 * 1.82)
WZ      = Sample.nanoAODfromDAS("WZ", "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM", dbFile=dbFile, xSection=47.13)
WZ_ext  = Sample.nanoAODfromDAS("WZ_ext", "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM", dbFile=dbFile, xSection=47.13)
ZZ      = Sample.nanoAODfromDAS("ZZ", "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM", dbFile=dbFile, xSection=16.523)
ZZ_ext  = Sample.nanoAODfromDAS("ZZ_ext", "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM", dbFile=dbFile, xSection=16.523)


boson = [
    WW,WW_ext,
    WZ,WZ_ext,
    ZZ,ZZ_ext,
    ]

## W+jets
WJetsToLNu      = Sample.nanoAODfromDAS("WJetsToLNu", "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM", dbFile=dbFile, xSection=3* 20508.9)
WJetsToLNu_ext  = Sample.nanoAODfromDAS("WJetsToLNu_ext", "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM", dbFile=dbFile, xSection=3* 20508.9)

wjets = [
    WJetsToLNu,
    WJetsToLNu_ext,
    ]

## rare

rare = [
    ]


## other

other = [
    ]

allSamples = DY + top + boson + wjets + rare + other

for s in allSamples:
    s.isData = False

