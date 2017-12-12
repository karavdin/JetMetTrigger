import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList

from JetMetTrigger.Analyzer.hltJetMETNtuple_cfi import *

process = cms.Process("hltJetMET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if runOnData:
  process.GlobalTag.globaltag = '92X_dataRun2_HLT_v7'
else:
  process.GlobalTag.globaltag = '92X_upgrade2017_TSG_For90XSamples_V2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(

	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/1A34E997-3F76-E711-A45C-002590D0B004.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/244A3B9E-3F76-E711-9855-0025905B85EE.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/28A17841-4076-E711-97C5-0CC47A4D75F0.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/28B8C208-4076-E711-A562-0242AC130002.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/80A744DB-3F76-E711-85C5-0CC47A6C1054.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/DAF0BB7C-9D75-E711-9D64-0025905A48D0.root',

	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/180ED238-3371-E711-AA34-02163E01A413.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/64FCD9F0-2471-E711-89D1-02163E019C9F.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/7C064643-3A71-E711-A956-02163E0119C5.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/90370974-2971-E711-94F0-02163E011E64.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/B2EB96A9-2B71-E711-AD25-02163E019DC2.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C0CE5CF7-2671-E711-A65F-02163E014218.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C81001BD-2C71-E711-BFAA-02163E01418D.root',
	'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/E03B5A6E-2F71-E711-B4D5-02163E01A4AD.root',
 
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/183EA4EB-3B6D-E711-81F3-02163E0135E4.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/249027CF-3F6D-E711-8781-02163E014115.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/24AECC6E-506D-E711-981F-02163E01A3F5.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/3A86E91B-436D-E711-99E1-02163E0143F0.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/44B6E402-456D-E711-B9DE-02163E01A6EC.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/76BF158B-416D-E711-9E06-02163E01A1C8.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/AA793AED-7B6D-E711-B793-02163E012BA6.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/C2951904-486D-E711-B0BD-02163E01A35F.root',
	

 
    
                                )
)

configureJetMetNtuple(process)
process.ntuple = cms.Path(process.JetMetNtupleSequence)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)

