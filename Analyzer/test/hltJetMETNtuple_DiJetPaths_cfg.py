import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList


isData = True
#isData = False
runOnData=isData #data/MC switch

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')

process = cms.Process("hltJetMET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if runOnData:
  process.GlobalTag.globaltag = '101X_upgrade2018_realistic_v7'
else:
  process.GlobalTag.globaltag = '101X_upgrade2018_realistic_v7'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(

	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/1A34E997-3F76-E711-A45C-002590D0B004.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/244A3B9E-3F76-E711-9855-0025905B85EE.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/28A17841-4076-E711-97C5-0CC47A4D75F0.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/28B8C208-4076-E711-A562-0242AC130002.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/80A744DB-3F76-E711-85C5-0CC47A6C1054.root',
	#'/store/mc/RunIISummer17MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v7-v1/110000/DAF0BB7C-9D75-E711-9D64-0025905A48D0.root',

	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/180ED238-3371-E711-AA34-02163E01A413.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/64FCD9F0-2471-E711-89D1-02163E019C9F.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/7C064643-3A71-E711-A956-02163E0119C5.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/90370974-2971-E711-94F0-02163E011E64.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/B2EB96A9-2B71-E711-AD25-02163E019DC2.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C0CE5CF7-2671-E711-A65F-02163E014218.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C81001BD-2C71-E711-BFAA-02163E01418D.root',
	#'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/E03B5A6E-2F71-E711-B4D5-02163E01A4AD.root',
 
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/183EA4EB-3B6D-E711-81F3-02163E0135E4.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/249027CF-3F6D-E711-8781-02163E014115.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/24AECC6E-506D-E711-981F-02163E01A3F5.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/3A86E91B-436D-E711-99E1-02163E0143F0.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/44B6E402-456D-E711-B9DE-02163E01A6EC.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/76BF158B-416D-E711-9E06-02163E01A1C8.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/AA793AED-7B6D-E711-B793-02163E012BA6.root',
	#'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/329/00000/C2951904-486D-E711-B0BD-02163E01A35F.root',
	

    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/105643A1-328F-E711-943D-02163E014641.root',
    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/1AE1B9B1-2F8F-E711-A3F6-02163E01441F.root',
    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/247DE91B-318F-E711-BC5C-02163E012AFE.root'
    '/store/data/Run2018B/JetHT/MINIAOD/PromptReco-v1/000/317/392/00000/FCD59F71-FF69-E811-A297-FA163EEA1CE0.root',
 
    
                                )
)

# Lumi filter ==========================================================================================

###process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_2017_HCAL_DCS_GOOD_post_JEC_bugfix.txt').getVLuminosityBlockRange()

# Electron ID ==========================================================================================

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff'
		 ]


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

### END Electron ID ====================================================================================


process.hltJetMetNtuple = cms.EDAnalyzer('HLTJetMETNtupleProducer',
                                         runJets = cms.bool(False),
                                         PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                         MetCollectionTag = cms.InputTag('slimmedMETs'),
					 applyMETFilters = cms.bool(options.applyMETFilters),
					 BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
					 BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
					 MuonCollectionTag = cms.InputTag('slimmedMuons'),
                                         ElectronCollectionTag = cms.InputTag('slimmedElectrons'),
                                         PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                         HLTPFJetCollectionTag = cms.InputTag('hltAK4PFJetsCorrected'),
                                         HLTCaloJetCollectionTag = cms.InputTag('hltAK4CaloJetsCorrected'),
                                         #eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                         #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                         #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                         #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                         #eleMvaNonTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                         #eleMvaNonTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                         #eleMvaTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
                                         #eleMvaTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
                                         #mvaNonTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                         #mvaNonTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                                         #mvaTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                                         #mvaTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
                                         eleMvaSpring16WPMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                         eleMvaSpring16WPTightMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                         mvaSpring16ValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                         mvaSpring16CategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                                         eleSummer16VetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                         eleSummer16LooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                         eleSummer16MediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                         eleSummer16TightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
					 hltprocess = cms.InputTag('TriggerResults::HLT'),
                                         triggerObjects = cms.InputTag("slimmedPatTrigger"),
                                         triggerPaths = cms.untracked.vstring('HLT_PFMET200_NotCleaned_v',					 
                                                                              'HLT_PFMET200_HBHECleaned_v',
									      'HLT_PFMET250_HBHECleaned_v',
									      'HLT_PFMET300_HBHECleaned_v',
                                                                              'HLT_PFMET200_HBHE_BeamHaloCleaned_v',
									      'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v', 
									      'HLT_PFMET110_PFMHT110_IDTight_v', 
                                                                              'HLT_Ele25_WPTight_Gsf_v',
                                                                              'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele27_WPTight_Gsf_v',
                                                                              'HLT_Ele27_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele30_WPTight_Gsf_v',
                                                                              'HLT_Ele30_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele32_WPTight_Gsf_v',
                                                                              'HLT_Ele32_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_PFJet40_v',
                                                                              'HLT_PFJet60_v',
                                                                              'HLT_PFJet80_v',
                                                                              'HLT_PFJet140_v',
                                                                              'HLT_PFJet200_v',
                                                                              'HLT_PFJet260_v',
                                                                              'HLT_PFJet320_v',
                                                                              'HLT_PFJet400_v',
                                                                              'HLT_PFJet450_v',
                                                                              'HLT_PFJet500_v',
                                                                              'HLT_PFJet550_v',
                                                                              'HLT_PFJetFwd40_v', 
                                                                              'HLT_PFJetFwd60_v',
                                                                              'HLT_PFJetFwd80_v',
                                                                              'HLT_PFJetFwd140_v',
                                                                              'HLT_PFJetFwd200_v',
                                                                              'HLT_PFJetFwd260_v',
                                                                              'HLT_PFJetFwd320_v',
                                                                              'HLT_PFJetFwd400_v',
                                                                              'HLT_PFJetFwd450_v',
                                                                              'HLT_PFJetFwd500_v'
                                                                              )
                                         )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)
process.p = cms.Path(process.egmGsfElectronIDSequence*
		     process.BadChargedCandidateFilter*
		     process.BadPFMuonFilter* 
		     process.hltJetMetNtuple
		     )
