#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
//#include "CMS_lumi.C"

void LegendSettings(TLegend *leg, int ncolumn){
  leg->SetNColumns(ncolumn);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
}

TGraphAsymmErrors* getMetEfficiency(std::string pathName_, TString inputFileName_)
{
  //float pTbins_[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 
  //210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 325., 350., 375., 400., 450., 600., 775., 1000.};
  float pTbins_[] = {0., 10., 20., 30., 40., 50., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 
		     130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400.};

  //TH1F* histMETAll = new TH1F("METAll", "METAll", 100, 0., 1000.);
  TH1F* histMETAll = new TH1F("METAll", "METAll", 38, pTbins_);
  histMETAll->Sumw2(); histMETAll->Reset();
  //TH1F* histMETPass = new TH1F("METPass", "METPass", 100, 0., 1000.);
  TH1F* histMETPass = new TH1F("METPass", "METPass", 38, pTbins_);
  histMETPass->Sumw2(); histMETPass->Reset();

  TFile *file = new TFile(inputFileName_);
  TTree* tree = (TTree*)file->Get("hltJetMetNtuple/tree");

  float metPt_, metPhi_;
  float caloMetPt_, caloMetPhi_;
  UInt_t passedL1MET_, passedHLTCaloMET_, passedHLTCaloMETClean_;
  std::vector<std::string>* triggerResults_ = new std::vector<std::string>();
  std::vector<float>* muonPt_ = new std::vector<float>();
  std::vector<float>* muonPx_ = new std::vector<float>();
  std::vector<float>* muonPy_ = new std::vector<float>();
  std::vector<float>* muonEta_ = new std::vector<float>();
  std::vector<float>* muonDxy_ = new std::vector<float>();
  std::vector<float>* muonDz_ = new std::vector<float>();
  std::vector<float>* muonR04SumChargedHadronPt_ = new std::vector<float>();
  std::vector<float>* muonR04SumNeutralHadronEt_ = new std::vector<float>();
  std::vector<float>* muonR04SumPhotonEt_ = new std::vector<float>();
  std::vector<float>* muonR04SumPUPt_ = new std::vector<float>();
  std::vector<bool>* muonIsICHEPMedium_ = new std::vector<bool>();
  std::vector<bool>* muonIsGlobal_ = new std::vector<bool>();
  std::vector<bool>* muonIsTracker_ = new std::vector<bool>();
  std::vector<bool>* muonIsPF_ = new std::vector<bool>();

  std::vector<float>* elecPt_ = new std::vector<float>(); 
  std::vector<float>* elecPx_ = new std::vector<float>();
  std::vector<float>* elecPy_ = new std::vector<float>();   
  std::vector<float>* elecEta_ = new std::vector<float>();
  std::vector<float>* elecDxy_ = new std::vector<float>();
  std::vector<float>* elecDz_ = new std::vector<float>();
  std::vector<bool>* elec_pass_conversion_ = new std::vector<bool>();
  std::vector<float>* elecR03SumChargedHadronPt_ = new std::vector<float>();
  std::vector<float>* elecR03SumNeutralHadronEt_ = new std::vector<float>();
  std::vector<float>* elecR03SumPhotonEt_ = new std::vector<float>();
  std::vector<float>* elecR03SumPUPt_ = new std::vector<float>();
  std::vector<bool>* elec_mva_tight_Spring16_v1_ = new std::vector<bool>();
  std::vector<bool>* elec_mva_medium_Spring16_v1_ = new std::vector<bool>();
  std::vector<bool>* elec_cutId_tight_Summer16_ = new std::vector<bool>();
  
//  std::vector<long>* run_ = new std::vector<long>();
//  std::vector<long>* event_ = new std::vector<long>(); 
//  std::vector<long>* lumi_ = new std::vector<long>(); 


  tree->SetBranchAddress("metPt", &metPt_);
  tree->SetBranchAddress("metPhi", &metPhi_);
  tree->SetBranchAddress("caloMetPt", &caloMetPt_);
  tree->SetBranchAddress("caloMetPhi", &caloMetPhi_);
  tree->SetBranchAddress("triggerResults", &triggerResults_);
  tree->SetBranchAddress("passedL1MET", &passedL1MET_);
  tree->SetBranchAddress("passedHLTCaloMET", &passedHLTCaloMET_);
  tree->SetBranchAddress("passedHLTCaloMETClean", &passedHLTCaloMETClean_);
  /*tree->SetBranchAddress("muonPt", &muonPt_);
  tree->SetBranchAddress("muonPx", &muonPx_);
  tree->SetBranchAddress("muonPy", &muonPy_);
  tree->SetBranchAddress("muonEta", &muonEta_);
  tree->SetBranchAddress("muonDxy", &muonDxy_);
  tree->SetBranchAddress("muonDz", &muonDz_);
  tree->SetBranchAddress("muonR04SumChargedHadronPt", &muonR04SumChargedHadronPt_);
  tree->SetBranchAddress("muonR04SumNeutralHadronEt", &muonR04SumNeutralHadronEt_);
  tree->SetBranchAddress("muonR04SumPhotonEt", &muonR04SumPhotonEt_);
  tree->SetBranchAddress("muonR04SumPUPt", &muonR04SumPUPt_);
  tree->SetBranchAddress("muonIsICHEPMedium", &muonIsICHEPMedium_);
  tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal_);
  tree->SetBranchAddress("muonIsTracker", &muonIsTracker_);
  tree->SetBranchAddress("muonIsPF", &muonIsPF_);
  
  
  tree->SetBranchAddress("elecPt", &elecPt_);
  tree->SetBranchAddress("elecPx", &elecPx_);
  tree->SetBranchAddress("elecPy", &elecPy_);
  tree->SetBranchAddress("elecEta", &elecEta_);
  tree->SetBranchAddress("elecDxy", &elecDxy_);
  tree->SetBranchAddress("elecDz", &elecDz_);
  tree->SetBranchAddress("elec_pass_conversion", &elec_pass_conversion_);
  tree->SetBranchAddress("elecR03SumChargedHadronPt", &elecR03SumChargedHadronPt_);
  tree->SetBranchAddress("elecR03SumNeutralHadronEt", &elecR03SumNeutralHadronEt_);
  tree->SetBranchAddress("elecR03SumPhotonEt", &elecR03SumPhotonEt_);
  tree->SetBranchAddress("elecR03SumPUPt", &elecR03SumPUPt_);
  tree->SetBranchAddress("elec_mva_tight_Spring16_v1", &elec_mva_tight_Spring16_v1_);
  tree->SetBranchAddress("elec_mva_medium_Spring16_v1", &elec_mva_medium_Spring16_v1_);
  tree->SetBranchAddress("elec_cutId_tight_Summer16", &elec_cutId_tight_Summer16_);
  */
//  tree->SetBranchAddress("run", &run_);
//  tree->SetBranchAddress("event", &event_);
//  tree->SetBranchAddress("lumi", &lumi_);

  int nEntries = tree->GetEntries();

  for (int n = 0; n < nEntries; n++){
   
    
    tree->GetEntry(n);
    /*
    //check singleElectron trigger
    bool passSingleElectron = false;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      if((*triggerResults_)[i].find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos ||
        (*triggerResults_)[i].find("HLT_Ele27_eta2p1_WPTight_Gsf_v") != std::string::npos)passSingleElectron = true;
    }
    if(!passSingleElectron) continue;
    
    //select events with at least one electron
    std::vector<float>selElecPx_; selElecPx_.clear();
    std::vector<float>selElecPy_; selElecPy_.clear();    
    for(size_t ie = 0; ie < elecPt_->size(); ie++){
     if((*elecPt_)[ie] > 30. && fabs((*elecEta_)[ie]) < 2.5 && 
	(!(fabs((*elecEta_)[ie]) > 1.4442 && fabs((*elecEta_)[ie]) < 1.5660)) && 
	fabs((*elecDxy_)[ie]) < 0.05 && fabs((*elecDz_)[ie]) < 0.1 && (*elec_pass_conversion_)[ie] == true) {
	float iso = ( (*elecR03SumChargedHadronPt_)[ie] + std::max(0., (*elecR03SumNeutralHadronEt_)[ie] + (*elecR03SumPhotonEt_)[ie] - (0.5 * (*elecR03SumPUPt_)[ie])) )/(*elecPt_)[ie];
	if(iso < 0.15){
	  //if((*elec_mva_tight_Spring16_v1_)[ie] == true) {
	  if((*elec_mva_medium_Spring16_v1_)[ie] == true) {
	  //if((*elec_cutId_tight_Summer16_)[ie] == true) {
     	   selElecPx_.push_back((*elecPx_)[ie]);
	   selElecPy_.push_back((*elecPy_)[ie]);
	  }
	 //}
	}
     }
    }
    if(selElecPx_.size() == 0) continue;
    
    //find selected muons
    std::vector<float>rejectMuon_; rejectMuon_.clear();
    std::vector<float>selMuonPx_; selMuonPx_.clear();
    std::vector<float>selMuonPy_; selMuonPy_.clear();
    for(size_t im = 0; im < muonPt_->size(); im++){
      if((*muonPt_)[im] > 10. && (*muonIsGlobal_)[im] == true){
	rejectMuon_.push_back((*muonPx_)[im]);
	float iso = ( (*muonR04SumChargedHadronPt_)[im] + std::max(0., (*muonR04SumNeutralHadronEt_)[im] + (*muonR04SumPhotonEt_)[im] - (0.5 * (*muonR04SumPUPt_)[im])) )/(*muonPt_)[im];
	if((*muonPt_)[im] > 30. && fabs((*muonEta_)[im]) < 2.4 && fabs((*muonDxy_)[im]) < 0.05 && fabs((*muonDz_)[im]) < 0.1 && iso < 0.15){
	  if((*muonIsICHEPMedium_)[im] == true) {
	  selMuonPx_.push_back((*muonPx_)[im]);
	  selMuonPy_.push_back((*muonPy_)[im]);
	}
       }
      }
    }

    if(rejectMuon_.size() > 0) continue;
    TLorentzVector newMetP4_;
    newMetP4_.SetPtEtaPhiM(metPt_, 0, metPhi_, 0);
    float newMetPx_ = newMetP4_.Px();
    float newMetPy_ = newMetP4_.Py();
    //Remove muons from MET to make it similar to CaloMet
    //for(size_t im = 0; im < selMuonPx_.size(); im++){
    //  newMetPx_ += selMuonPx_[im];
    //  newMetPy_ += selMuonPy_[im];
    //}
    float newMetPt_ = sqrt(newMetPx_*newMetPx_ + newMetPy_*newMetPy_);
    */
    /*
    if(passedL1MET_ && passedHLTCaloMET_ && passedHLTCaloMETClean_)
      histMETAll->Fill(metPt_);
    bool passTrigger = false;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      if((*triggerResults_)[i].find(pathName_) != std::string::npos)passTrigger = true;
    }
    if(passTrigger)histMETPass->Fill(metPt_);
    //if(selMuonPx_.size() > 0 && metPt_ > 100) std::cout<<" old met "<<metPt_<<" new met "<<newMetPt_<<std::endl;
    */


    histMETAll->Fill(caloMetPt_);                                                                                                                                                                          
    bool passTrigger = false;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      if((*triggerResults_)[i].find(pathName_) != std::string::npos)passTrigger = true;
    }
    if(passTrigger)histMETPass->Fill(caloMetPt_);
  }


    
  TGraphAsymmErrors* histMETEff = new TGraphAsymmErrors(histMETPass, histMETAll);
  
  return histMETEff;
}

void plotMetEfficiency(std::vector<std::string> pathNames, std::vector<TString> filenames, std::vector<TString> legends, TString outFileName_, TString label_string)
{
  gROOT->ProcessLine(".L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  //gROOT->LoadMacro("CMS_lumi.C");
  
  //writeExtraText = true;       // if extra text
  //extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  //lumi_sqrtS = "2016, 13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  int colors[] = { 1, 2, 4, 6, 7, 9};
  int style[] = {20, 21, 22, 23, 24, 25};

  std::vector<TGraphAsymmErrors*> histMETEffs;
  histMETEffs.resize(filenames.size());
  for(size_t ifile = 0; ifile < filenames.size(); ifile++){
    histMETEffs[ifile] = (TGraphAsymmErrors*)getMetEfficiency(pathNames[ifile], filenames[ifile]);
    histMETEffs[ifile]->SetMarkerStyle(style[ifile]);
    histMETEffs[ifile]->SetMarkerColor(colors[ifile]);
    histMETEffs[ifile]->SetMarkerSize(1.0);
    histMETEffs[ifile]->SetLineColor(colors[ifile]);
  }


  TCanvas* c1 = new TCanvas();
  c1->SetGridx();
  c1->SetGridy();
  c1->cd();
  
  //TH1F* histogram_base = new TH1F("histogram_base", "", 100, 1., 999.);
  TH1F* histogram_base = new TH1F("histogram_base", "", 40, 1., 401.);
  histogram_base->SetTitle("");
  histogram_base->SetStats(false);
  histogram_base->SetMinimum(0.0);
  histogram_base->SetMaximum(1.2);
  histogram_base->GetXaxis()->SetTitle("Offline Calo E_{T}^{miss} (GeV)");
  histogram_base->GetYaxis()->SetTitle("Efficiency");
  histogram_base->Draw("hist");
  
  for(size_t ih = 0; ih < histMETEffs.size(); ih++){ 
    histMETEffs[ih]->Draw("Psame");
  }

  TLegend *leg1 = new TLegend(0.5,0.6,0.9,(0.6+0.05*(histMETEffs.size())));
  LegendSettings(leg1, 1);
  for(size_t ih = 0; ih < histMETEffs.size(); ih++){ 
    leg1->AddEntry(histMETEffs[ih], legends[ih], "lp");
  }
  leg1->Draw();

  TPaveText* label = new TPaveText(0.2, 0.85, 0.65, 0.9, "brNDC");
  label->AddText(label_string);
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(0.03);
  label->SetTextFont(42);
  label->Draw();

  //CMS_lumi( c1, iPeriod, 11 );
  
  c1->Update();
  c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+".png");
  c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+".pdf");
  c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+".jpg");
}

void plotMetEfficiencyAll()
{

  std::vector<std::string> pathNames_PFMET;
  pathNames_PFMET.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned_v");
  pathNames_PFMET.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned_v");
  pathNames_PFMET.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned_v");

  /*
  vector<TString> filenames_PFMET;
  filenames_PFMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiOFF.root");
  filenames_PFMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiON.root");
  filenames_PFMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_2_100X_upgrade2018_realistic_v10.root");

  vector<TString>legends_PFMET;
  legends_PFMET.push_back("10_0_0, mahi off");
  legends_PFMET.push_back("10_0_0, mahi on");
  legends_PFMET.push_back("10_0_2, m3+mahi");

  plotMetEfficiency(pathNames_PFMET, filenames_PFMET, legends_PFMET, "ADDMonojet_HLT_PFMET200_HBHE_BeamHaloCleaned_OnlyPFatHLT", "HLT_PFMET200_HBHE_BeamHaloCleaned");

  std::vector<std::string> pathNames_PFMETTypeOne;
  pathNames_PFMETTypeOne.push_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v");
  pathNames_PFMETTypeOne.push_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v");
  pathNames_PFMETTypeOne.push_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v");
  
  vector<TString> filenames_PFMETTypeOne;
  filenames_PFMETTypeOne.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiOFF.root");
  filenames_PFMETTypeOne.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiON.root");
  filenames_PFMETTypeOne.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_2_100X_upgrade2018_realistic_v10.root");

  vector<TString>legends_PFMETTypeOne;
  legends_PFMETTypeOne.push_back("10_0_0, mahi off");
  legends_PFMETTypeOne.push_back("10_0_0, mahi on");
  legends_PFMETTypeOne.push_back("10_0_2, m3+mahi");

  plotMetEfficiency(pathNames_PFMETTypeOne, filenames_PFMETTypeOne, legends_PFMETTypeOne, "ADDMonojet_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_OnlyPFatHLT", "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned");
  */

  std::vector<std::string> pathNames_CaloMET;
  pathNames_CaloMET.push_back("HLT_CaloMET100_NotCleaned_v");
  pathNames_CaloMET.push_back("HLT_CaloMET100_NotCleaned_v");
  pathNames_CaloMET.push_back("HLT_CaloMET100_NotCleaned_v");


  vector<TString> filenames_CaloMET;
  filenames_CaloMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiOFF.root");
  filenames_CaloMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiON.root");
  filenames_CaloMET.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_2_100X_upgrade2018_realistic_v10.root");

  vector<TString>legends_CaloMET;
  legends_CaloMET.push_back("10_0_0, mahi off");
  legends_CaloMET.push_back("10_0_0, mahi on");
  legends_CaloMET.push_back("10_0_2, m3+mahi");

  plotMetEfficiency(pathNames_CaloMET, filenames_CaloMET, legends_CaloMET, "ADDMonojet_HLT_CaloMET100_NotCleaned", "HLT_CaloMET100_NotCleaned");
    
  std::vector<std::string> pathNames_CaloMETClean;
  pathNames_CaloMETClean.push_back("HLT_CaloMET100_HBHECleaned_v");
  pathNames_CaloMETClean.push_back("HLT_CaloMET100_HBHECleaned_v");
  pathNames_CaloMETClean.push_back("HLT_CaloMET100_HBHECleaned_v");


  vector<TString> filenames_CaloMETClean;
  filenames_CaloMETClean.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiOFF.root");
  filenames_CaloMETClean.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_0-100X_upgrade2018_realistic_v6_mahiON.root");
  filenames_CaloMETClean.push_back("../test/hltJetMetNtuple_RelValADDMonojet_10_0_2_100X_upgrade2018_realistic_v10.root");

  vector<TString>legends_CaloMETClean;
  legends_CaloMETClean.push_back("10_0_0, mahi off");
  legends_CaloMETClean.push_back("10_0_0, mahi on");
  legends_CaloMETClean.push_back("10_0_2, m3+mahi");

  plotMetEfficiency(pathNames_CaloMETClean, filenames_CaloMETClean, legends_CaloMETClean, "ADDMonojet_HLT_CaloMET100_HBHECleaned", "HLT_CaloMET100_HBHECleaned"); 
}
