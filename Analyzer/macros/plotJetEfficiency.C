#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
//#include "CMS_lumi.C"

double Pi = 3.14159265359;
void LegendSettings(TLegend *leg, int ncolumn){
  leg->SetNColumns(ncolumn);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  float deta = fabs(eta1 - eta2);
  float dphi = fabs(phi1 - phi2);
  if(dphi > Pi)dphi = 2*Pi - dphi;
  float dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

TGraphAsymmErrors* getJetEfficiency(float hltThreshold_, float etaMin_, float etaMax_, TString inputFileName_)
{
  float pTbins_[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 
		     210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 325., 350., 375., 400., 450., 600., 775., 1000.};

  TH1F* histJetAll = new TH1F("JetAll", "JetAll", 100, 0., 500.);
  //TH1F* histJetAll = new TH1F("JetAll", "JetAll", 38, pTbins_);
  histJetAll->Sumw2(); histJetAll->Reset();
  TH1F* histJetPass = new TH1F("JetPass", "JetPass", 100, 0., 500.);
  //TH1F* histJetPass = new TH1F("JetPass", "JetPass", 38, pTbins_);
  histJetPass->Sumw2(); histJetPass->Reset();

  TFile *file = new TFile(inputFileName_);
  TTree* tree = (TTree*)file->Get("hltJetMetNtuple/tree");

  std::vector<float>* pfjetPt_ = new std::vector<float>();
  std::vector<float>* pfjetEta_ = new std::vector<float>();
  std::vector<float>* pfjetPhi_ = new std::vector<float>();
  std::vector<float>* hltpfjetPt_ = new std::vector<float>();
  std::vector<float>* hltpfjetEta_ = new std::vector<float>();
  std::vector<float>* hltpfjetPhi_ = new std::vector<float>();
  //std::vector<float>* hltcalojetPt_ = new std::vector<float>();
  //std::vector<float>* hltcalojetEta_ = new std::vector<float>();
  //std::vector<float>* hltcalojetPhi_ = new std::vector<float>();

  std::vector<std::string>* triggerResults_ = new std::vector<std::string>();
  std::vector<float>* muonPt_ = new std::vector<float>();
  std::vector<float>* muonEta_ = new std::vector<float>();
  std::vector<float>* muonPhi_ = new std::vector<float>();
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
  std::vector<float>* elecEta_ = new std::vector<float>();
  std::vector<float>* elecPhi_ = new std::vector<float>();
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
  
  tree->SetBranchAddress("pfjetPt", &pfjetPt_);
  tree->SetBranchAddress("pfjetEta", &pfjetEta_);
  tree->SetBranchAddress("pfjetPhi", &pfjetPhi_);
  tree->SetBranchAddress("hltpfjetPt", &hltpfjetPt_);
  tree->SetBranchAddress("hltpfjetEta", &hltpfjetEta_);
  tree->SetBranchAddress("hltpfjetPhi", &hltpfjetPhi_);
  //tree->SetBranchAddress("hltcalojetPt", &hltcalojetPt_);
  //tree->SetBranchAddress("hltcalojetEta", &hltcalojetEta_);
  //tree->SetBranchAddress("hltcalojetPhi", &hltcalojetPhi_);
  
  tree->SetBranchAddress("triggerResults", &triggerResults_);
  tree->SetBranchAddress("muonPt", &muonPt_);
  tree->SetBranchAddress("muonEta", &muonEta_);
  tree->SetBranchAddress("muonPhi", &muonPhi_);
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
  tree->SetBranchAddress("elecEta", &elecEta_);
  tree->SetBranchAddress("elecPhi", &elecPhi_);
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
  
  int nEntries = tree->GetEntries();
  std::cout<<" no. of events "<<nEntries<<std::endl;

  for (int n = 0; n < nEntries; n++){
    tree->GetEntry(n);

    //check singleElectron trigger
    bool passSingleElectron = false;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      if((*triggerResults_)[i].find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos ||
        (*triggerResults_)[i].find("HLT_Ele27_eta2p1_WPTight_Gsf_v") != std::string::npos)passSingleElectron = true;
    }
    //if(!passSingleElectron) continue;
    
    //select events with at least one electron
    std::vector<float>selElecEta_; selElecEta_.clear();
    std::vector<float>selElecPhi_; selElecPhi_.clear();    
    for(size_t ie = 0; ie < elecPt_->size(); ie++){
     if((*elecPt_)[ie] > 30. && fabs((*elecEta_)[ie]) < 2.5 && 
	(!(fabs((*elecEta_)[ie]) > 1.4442 && fabs((*elecEta_)[ie]) < 1.5660)) && 
	fabs((*elecDxy_)[ie]) < 0.05 && fabs((*elecDz_)[ie]) < 0.1 && (*elec_pass_conversion_)[ie] == true) {
	float iso = ( (*elecR03SumChargedHadronPt_)[ie] + std::max(0., (*elecR03SumNeutralHadronEt_)[ie] + (*elecR03SumPhotonEt_)[ie] - (0.5 * (*elecR03SumPUPt_)[ie])) )/(*elecPt_)[ie];
	if(iso < 0.15){
	  //if((*elec_mva_tight_Spring16_v1_)[ie] == true) {
	  if((*elec_mva_medium_Spring16_v1_)[ie] == true) {
	  //if((*elec_cutId_tight_Summer16_)[ie] == true) {
     	   selElecEta_.push_back((*elecEta_)[ie]);
	   selElecPhi_.push_back((*elecPhi_)[ie]);
	  }
	}
     }
    }
    
    //find selected muons
    std::vector<float>selMuonEta_; selMuonEta_.clear();
    std::vector<float>selMuonPhi_; selMuonPhi_.clear();
    for(size_t im = 0; im < muonPt_->size(); im++){
      if((*muonPt_)[im] > 10. && (*muonIsGlobal_)[im] == true){
	float iso = ( (*muonR04SumChargedHadronPt_)[im] + std::max(0., (*muonR04SumNeutralHadronEt_)[im] + (*muonR04SumPhotonEt_)[im] - (0.5 * (*muonR04SumPUPt_)[im])) )/(*muonPt_)[im];
	if((*muonPt_)[im] > 25. && fabs((*muonEta_)[im]) < 2.4 && fabs((*muonDxy_)[im]) < 0.05 && fabs((*muonDz_)[im]) < 0.1 && iso < 0.15){
	  if((*muonIsICHEPMedium_)[im] == true) {
	  selMuonEta_.push_back((*muonEta_)[im]);
	  selMuonPhi_.push_back((*muonPhi_)[im]);
	}
       }
      }
    }

    //Require at least one lepton
    //if(selElecEta_.size()+selMuonEta_.size() <= 0) continue;
    
    //Loop over pfJets and find the highest pT jet in event 
    float highestPt_ = 0; int jet_index = -1;
    for(size_t ij = 0; ij < pfjetPt_->size(); ij++){
      if((*pfjetPt_)[ij] > highestPt_ && fabs((*pfjetEta_)[ij]) > etaMin_ && fabs((*pfjetEta_)[ij]) < etaMax_){
	//clean from leptons
	bool match_mu = false; bool match_elec = false;
	for(size_t im = 0; im < selMuonEta_.size(); im++){
	  if(deltaR((*pfjetEta_)[ij], (*pfjetPhi_)[ij], selMuonEta_[im], selMuonPhi_[im]) < 0.3)match_mu = true;
	}
	for(size_t ie = 0; ie < selElecEta_.size(); ie++){
          if(deltaR((*pfjetEta_)[ij], (*pfjetPhi_)[ij], selElecEta_[ie], selElecPhi_[ie]) < 0.3)match_elec = true;
        }
	if(match_mu || match_elec) continue;

	highestPt_ = (*pfjetPt_)[ij];
	jet_index = ij;
      }
    }
    
    //Go to next event if there are no jets
    if(highestPt_ <= 0) continue;
    histJetAll->Fill(highestPt_);
    //std::cout<<" highest Jet pT "<<highestPt_<<endl;
    
    //Find matching to hltJet passing HLT threshold
    bool Found_Match = false;
    for(size_t ij = 0; ij < hltpfjetPt_->size(); ij++){
      if((*hltpfjetPt_)[ij] > hltThreshold_){
	if(deltaR((*hltpfjetEta_)[ij], (*hltpfjetPhi_)[ij], (*pfjetEta_)[jet_index], (*pfjetPhi_)[jet_index]) < 0.3)Found_Match = true;
      }
    }

    //Fill the numerator
    if(Found_Match)histJetPass->Fill(highestPt_);
    
  }
  std::cout<<" Total Jets "<<histJetAll->Integral()<<" HLT passed "<<histJetPass->Integral()<<std::endl;
  TGraphAsymmErrors* histJetEff = new TGraphAsymmErrors(histJetPass, histJetAll);
  
  return histJetEff;
}

void plotJetEfficiency(std::vector<float> hltThresholds_, vector<float> etaMins_, vector<float> etaMaxs_, std::vector<TString> filenames, std::vector<TString> legends, TString outFileName_, TString label_string)
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

  std::vector<TGraphAsymmErrors*> histJetEffs;
  histJetEffs.resize(filenames.size());
  for(size_t ifile = 0; ifile < filenames.size(); ifile++){
    histJetEffs[ifile] = (TGraphAsymmErrors*)getJetEfficiency(hltThresholds_[ifile], etaMins_[ifile], etaMaxs_[ifile], filenames[ifile]);
    histJetEffs[ifile]->SetMarkerStyle(style[ifile]);
    histJetEffs[ifile]->SetMarkerColor(colors[ifile]);
    histJetEffs[ifile]->SetMarkerSize(1.0);
    histJetEffs[ifile]->SetLineColor(colors[ifile]);
  }


  TCanvas* c1 = new TCanvas();
  c1->SetGridx();
  c1->SetGridy();
  c1->cd();
  
  TH1F* histogram_base = new TH1F("histogram_base", "", 100, 0., 500.);
  histogram_base->SetTitle("");
  histogram_base->SetStats(false);
  histogram_base->SetMinimum(0.0);
  histogram_base->SetMaximum(1.2);
  histogram_base->GetXaxis()->SetTitle("Offline Jet p_{T} (GeV)");
  histogram_base->GetYaxis()->SetTitle("Efficiency");
  histogram_base->Draw("hist");
  
  for(size_t ih = 0; ih < histJetEffs.size(); ih++){ 
    histJetEffs[ih]->Draw("Psame");
  }

  TLegend *leg1 = new TLegend(0.5,0.2,0.9,(0.2+0.05*(histJetEffs.size())));
  LegendSettings(leg1, 1);
  for(size_t ih = 0; ih < histJetEffs.size(); ih++){ 
    leg1->AddEntry(histJetEffs[ih], legends[ih], "lp");
  }
  leg1->Draw();
/*
  TPaveText* label = new TPaveText(0.2, 0.85, 0.65, 0.9, "brNDC");
  label->AddText(label_string);
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(0.03);
  label->SetTextFont(42);
  label->Draw();
*/  
  //CMS_lumi( c1, iPeriod, 11 );
  
  c1->Update();
  c1->SaveAs("plots_test/Efficiency_"+TString(outFileName_)+"_sel.png");
  c1->SaveAs("plots_test/Efficiency_"+TString(outFileName_)+"_sel.pdf");
  c1->SaveAs("plots_test/Efficiency_"+TString(outFileName_)+"_sel.jpg");
}

void plotJetEfficiencyAll()
{

  std::vector<float> hltThresholds;
  hltThresholds.push_back(40.);
  hltThresholds.push_back(60.);
  hltThresholds.push_back(80.);
  hltThresholds.push_back(120.);
  hltThresholds.push_back(200.);
  
  std::vector<float> etaMins;
  etaMins.push_back(0.);
  etaMins.push_back(0.);
  etaMins.push_back(0.);
  etaMins.push_back(0.);
  etaMins.push_back(0.);

  std::vector<float> etaMaxs;
  etaMaxs.push_back(5.);
  etaMaxs.push_back(5.);
  etaMaxs.push_back(5.);
  etaMaxs.push_back(5.);
  etaMaxs.push_back(5.);

  vector<TString> filenames_data;
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");

  vector<TString>legends_data;
  legends_data.push_back("HLT_PFJet40");
  legends_data.push_back("HLT_PFJet60");
  legends_data.push_back("HLT_PFJet80");
  legends_data.push_back("HLT_PFJet120");
  legends_data.push_back("HLT_PFJet200");

  //plotJetEfficiency(hltThresholds, etaMins, etaMaxs, filenames_data, legends_data, "HLT_PFJet_Test", "");
  
  hltThresholds.clear();
  hltThresholds.push_back(40.);
  hltThresholds.push_back(40.);
  hltThresholds.push_back(40.);
  hltThresholds.push_back(40.);
  hltThresholds.push_back(40.);
  hltThresholds.push_back(40.);

  etaMins.clear();
  etaMins.push_back(0.);
  etaMins.push_back(1.3);
  etaMins.push_back(2.0);
  etaMins.push_back(2.5);
  etaMins.push_back(2.75);
  etaMins.push_back(3.0);

  etaMaxs.clear();
  etaMaxs.push_back(1.3);
  etaMaxs.push_back(2.0);
  etaMaxs.push_back(2.5);
  etaMaxs.push_back(2.75);
  etaMaxs.push_back(3.0);
  etaMaxs.push_back(5.0);

  filenames_data.clear();
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");
  filenames_data.push_back("/tmp/mdjordje/hltJetMetNtuple_171.root");

  legends_data.clear();
  legends_data.push_back("HLT_PFJet40, 0 < #eta < 1.3");
  legends_data.push_back("HLT_PFJet40, 1.3 < #eta < 2.0");
  legends_data.push_back("HLT_PFJet40, 2.0 < #eta < 2.5");
  legends_data.push_back("HLT_PFJet40, 2.5 < #eta < 2.75");
  legends_data.push_back("HLT_PFJet40, 2.75 < #eta < 3.0");
  legends_data.push_back("HLT_PFJet40, 3.0 < #eta < 5.0");

  plotJetEfficiency(hltThresholds, etaMins, etaMaxs, filenames_data, legends_data, "HLT_PFJet40_Eta", ""); 
  
}
