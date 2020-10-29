
#include<TApplication.h>
#include<TFile.h>
#include<TMath.h>
#include<TROOT.h>
#include<TSystem.h>
#include<TTree.h>
#include<TVector2.h>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraph.h>
#include<TH2F.h>
#include<TLatex.h>
#include<TLegend.h>
#include<TLine.h>
#include<TPaveStats.h>
#include<TStyle.h>

#include<algorithm>
#include<chrono>
#include<ctime>
#include<dirent.h>
#include<fstream>
#include<iostream>
#include<mutex>
#include<string>
#include<sys/stat.h>
#include<sys/types.h>
#include<thread>
#include<unistd.h>
#include<vector>

#ifdef __MAKECINT__ 
#pragma link C++ class vector<vector<double> >+;  
#pragma link C++ class vector<vector<float> >+; 
#pragma link C++ class vector<vector<int> >+;  
#pragma link C++ class vector<vector<bool> >+; 
#pragma link C++ class vector<vector<map<int,float>> >+; 
#endif
    

/* Histograms */
   vector<string> binning = {"","_EB","_EE"};

   TH1F *h_goodSC_R9_new[3];
   TH1F *h_goodSC_SigmaIetaIeta_new[3];
   TH1F *h_goodSC_SigmaIphiIphi_new[3];
   TH1F *h_goodSC_fullR9_new[3];
   TH1F *h_goodSC_fullSigmaIetaIeta_new[3];
   TH1F *h_goodSC_fullSigmaIphiIphi_new[3];
   TH1F *h_goodSC_etaWidth_new[3];
   TH1F *h_goodSC_phiWidth_new[3];
   TH1F *h_fakeSC_R9_new[3];
   TH1F *h_fakeSC_SigmaIetaIeta_new[3];
   TH1F *h_fakeSC_SigmaIphiIphi_new[3];
   TH1F *h_fakeSC_fullR9_new[3];
   TH1F *h_fakeSC_fullSigmaIetaIeta_new[3];
   TH1F *h_fakeSC_fullSigmaIphiIphi_new[3];
   TH1F *h_fakeSC_etaWidth_new[3];
   TH1F *h_fakeSC_phiWidth_new[3];

   TH1F *h_goodSC_eta_new[3];
   TH1F *h_goodSC_phi_new[3];
   TH1F *h_goodSC_et_new[3];
   TH1F *h_goodSC_energy_new[3];
   
   TH1F *h_fakeSC_eta_new[3];
   TH1F *h_fakeSC_phi_new[3];
   TH1F *h_fakeSC_et_new[3];
   TH1F *h_fakeSC_energy_new[3];

   TH1F *h_goodFakeSC_deta_new[3];
   TH1F *h_goodFakeSC_dphi_new[3];
   TH1F *h_goodFakeSC_detrel_new[3];
   TH1F *h_goodFakeSC_denergyrel_new[3];

   TH1F *h_goodSC_R9_old[3];
   TH1F *h_goodSC_SigmaIetaIeta_old[3];
   TH1F *h_goodSC_SigmaIphiIphi_old[3];
   TH1F *h_goodSC_fullR9_old[3];
   TH1F *h_goodSC_fullSigmaIetaIeta_old[3];
   TH1F *h_goodSC_fullSigmaIphiIphi_old[3];
   TH1F *h_goodSC_etaWidth_old[3];
   TH1F *h_goodSC_phiWidth_old[3];
   TH1F *h_fakeSC_R9_old[3];
   TH1F *h_fakeSC_SigmaIetaIeta_old[3];
   TH1F *h_fakeSC_SigmaIphiIphi_old[3];
   TH1F *h_fakeSC_fullR9_old[3];
   TH1F *h_fakeSC_fullSigmaIetaIeta_old[3];
   TH1F *h_fakeSC_fullSigmaIphiIphi_old[3];
   TH1F *h_fakeSC_etaWidth_old[3];
   TH1F *h_fakeSC_phiWidth_old[3];

   TH1F *h_goodSC_eta_old[3];
   TH1F *h_goodSC_phi_old[3];
   TH1F *h_goodSC_et_old[3];
   TH1F *h_goodSC_energy_old[3];
   
   TH1F *h_fakeSC_eta_old[3];
   TH1F *h_fakeSC_phi_old[3];
   TH1F *h_fakeSC_et_old[3];
   TH1F *h_fakeSC_energy_old[3];

   TH1F *h_goodFakeSC_deta_old[3];
   TH1F *h_goodFakeSC_dphi_old[3];
   TH1F *h_goodFakeSC_detrel_old[3];
   TH1F *h_goodFakeSC_denergyrel_old[3];

    //integrated
   TH1F *h_goodSC_R9_new_integrated[3][3];
   TH1F *h_goodSC_SigmaIetaIeta_new_integrated[3][3];
   TH1F *h_goodSC_SigmaIphiIphi_new_integrated[3][3];
   TH1F *h_goodSC_fullR9_new_integrated[3][3];
   TH1F *h_goodSC_fullSigmaIetaIeta_new_integrated[3][3];
   TH1F *h_goodSC_fullSigmaIphiIphi_new_integrated[3][3];
   TH1F *h_goodSC_etaWidth_new_integrated[3][3];
   TH1F *h_goodSC_phiWidth_new_integrated[3][3];
   TH1F *h_fakeSC_R9_new_integrated[3][3];
   TH1F *h_fakeSC_SigmaIetaIeta_new_integrated[3][3];
   TH1F *h_fakeSC_SigmaIphiIphi_new_integrated[3][3];
   TH1F *h_fakeSC_fullR9_new_integrated[3][3];
   TH1F *h_fakeSC_fullSigmaIetaIeta_new_integrated[3][3];
   TH1F *h_fakeSC_fullSigmaIphiIphi_new_integrated[3][3];
   TH1F *h_fakeSC_etaWidth_new_integrated[3][3];
   TH1F *h_fakeSC_phiWidth_new_integrated[3][3];

   TH1F *h_goodSC_eta_new_integrated[3][3];
   TH1F *h_goodSC_phi_new_integrated[3][3];
   TH1F *h_goodSC_et_new_integrated[3][3];
   TH1F *h_goodSC_energy_new_integrated[3][3];
   
   TH1F *h_fakeSC_eta_new_integrated[3][3];
   TH1F *h_fakeSC_phi_new_integrated[3][3];
   TH1F *h_fakeSC_et_new_integrated[3][3];
   TH1F *h_fakeSC_energy_new_integrated[3][3];

   TH1F *h_goodFakeSC_deta_new_integrated[3][3];
   TH1F *h_goodFakeSC_dphi_new_integrated[3][3];
   TH1F *h_goodFakeSC_detrel_new_integrated[3][3];
   TH1F *h_goodFakeSC_denergyrel_new_integrated[3][3];

   TH1F *h_goodSC_R9_old_integrated[3][3];
   TH1F *h_goodSC_SigmaIetaIeta_old_integrated[3][3];
   TH1F *h_goodSC_SigmaIphiIphi_old_integrated[3][3];
   TH1F *h_goodSC_fullR9_old_integrated[3][3];
   TH1F *h_goodSC_fullSigmaIetaIeta_old_integrated[3][3];
   TH1F *h_goodSC_fullSigmaIphiIphi_old_integrated[3][3];
   TH1F *h_goodSC_etaWidth_old_integrated[3][3];
   TH1F *h_goodSC_phiWidth_old_integrated[3][3];
   TH1F *h_fakeSC_R9_old_integrated[3][3];
   TH1F *h_fakeSC_SigmaIetaIeta_old_integrated[3][3];
   TH1F *h_fakeSC_SigmaIphiIphi_old_integrated[3][3];
   TH1F *h_fakeSC_fullR9_old_integrated[3][3];
   TH1F *h_fakeSC_fullSigmaIetaIeta_old_integrated[3][3];
   TH1F *h_fakeSC_fullSigmaIphiIphi_old_integrated[3][3];
   TH1F *h_fakeSC_etaWidth_old_integrated[3][3];
   TH1F *h_fakeSC_phiWidth_old_integrated[3][3];

   TH1F *h_goodSC_eta_old_integrated[3][3];
   TH1F *h_goodSC_phi_old_integrated[3][3];
   TH1F *h_goodSC_et_old_integrated[3][3];
   TH1F *h_goodSC_energy_old_integrated[3][3];
   
   TH1F *h_fakeSC_eta_old_integrated[3][3];
   TH1F *h_fakeSC_phi_old_integrated[3][3];
   TH1F *h_fakeSC_et_old_integrated[3][3];
   TH1F *h_fakeSC_energy_old_integrated[3][3];

   TH1F *h_goodFakeSC_deta_old_integrated[3][3];
   TH1F *h_goodFakeSC_dphi_old_integrated[3][3];
   TH1F *h_goodFakeSC_detrel_old_integrated[3][3];
   TH1F *h_goodFakeSC_denergyrel_old_integrated[3][3];

int _20to40_events,
    _20toInf_events,
    _40toInf_events;
   
float _20to40_xs = 232.9,
      _20toInf_xs = 3186.0,
      _40toInf_xs = 878.1; //pb

string cur_time(){
//returns the current time, for logging
    std::time_t tt = std::time(NULL);
    std::string s = std::ctime(&tt);
    return s.substr(0, s.size()-1);
}

void ReadInfile(string SuperClusterRef, string SuperClusterVal){
    cout<<"Reading input file"<<endl;


    TFile *hist_infile_SCRef;
    TFile *hist_infile_SCVal;

    string ref_1 = "", ref_2 = "", val_1 = "", val_2 = "";
    if(SuperClusterRef == "SingleSet"){
        hist_infile_SCRef = TFile::Open("./SCAlgoPlots_SingleSetParams.root");
        ref_1 = "DeepSC";
        ref_2 = "new";
    }
    else if(SuperClusterRef == "Mustache"){
        hist_infile_SCRef = TFile::Open("./SCAlgoPlots_SingleSetParams.root");
        ref_1 = "SC";
        ref_2 = "old";
    }
    if(SuperClusterVal == "SingleSet"){
        hist_infile_SCVal = TFile::Open("./SCAlgoPlots_SingleSetParams.root");
        val_1 = "DeepSC";
        val_2 = "new";
    }
    else if(SuperClusterVal == "LocalSets"){
        hist_infile_SCVal = TFile::Open("./SCAlgoPlots_LocalSetParams.root");
        val_1 = "DeepSC";
        val_2 = "new";
    }

    for(int i=0; i<3; i++){

        hist_infile_SCVal->GetObject(("h_goodSC_R9"+val_2+binning[i]).c_str(), h_goodSC_R9_new[i]);h_goodSC_R9_new[i]->SetName(("h_goodSC_R9_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_SigmaIetaIeta"+val_2+binning[i]).c_str(), h_goodSC_SigmaIetaIeta_new[i]);h_goodSC_SigmaIetaIeta_new[i]->SetName(("h_goodSC_SigmaIetaIeta_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_SigmaIphiIphi"+val_2+binning[i]).c_str(), h_goodSC_SigmaIphiIphi_new[i]);h_goodSC_SigmaIphiIphi_new[i]->SetName(("h_goodSC_SigmaIphiIphi_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_fullR9"+val_2+binning[i]).c_str(), h_goodSC_fullR9_new[i]);h_goodSC_fullR9_new[i]->SetName(("h_goodSC_fullR9_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_fullSigmaIetaIeta"+val_2+binning[i]).c_str(), h_goodSC_fullSigmaIetaIeta_new[i]);h_goodSC_fullSigmaIetaIeta_new[i]->SetName(("h_goodSC_fullSigmaIetaIeta_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_fullSigmaIphiIphi"+val_2+binning[i]).c_str(), h_goodSC_fullSigmaIphiIphi_new[i]);h_goodSC_fullSigmaIphiIphi_new[i]->SetName(("h_goodSC_fullSigmaIphiIphi_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_etaWidth"+val_2+binning[i]).c_str(), h_goodSC_etaWidth_new[i]);h_goodSC_etaWidth_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_phiWidth"+val_2+binning[i]).c_str(), h_goodSC_phiWidth_new[i]);h_goodSC_phiWidth_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_R9"+val_2+binning[i]).c_str(), h_fakeSC_R9_new[i]);h_fakeSC_R9_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_SigmaIetaIeta"+val_2+binning[i]).c_str(), h_fakeSC_SigmaIetaIeta_new[i]);h_fakeSC_SigmaIetaIeta_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_SigmaIphiIphi"+val_2+binning[i]).c_str(), h_fakeSC_SigmaIphiIphi_new[i]);h_fakeSC_SigmaIphiIphi_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_fullR9"+val_2+binning[i]).c_str(), h_fakeSC_fullR9_new[i]);h_fakeSC_fullR9_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_fullSigmaIetaIeta"+val_2+binning[i]).c_str(), h_fakeSC_fullSigmaIetaIeta_new[i]);h_fakeSC_fullSigmaIetaIeta_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_fullSigmaIphiIphi"+val_2+binning[i]).c_str(), h_fakeSC_fullSigmaIphiIphi_new[i]);h_fakeSC_fullSigmaIphiIphi_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_etaWidth"+val_2+binning[i]).c_str(), h_fakeSC_etaWidth_new[i]);h_fakeSC_etaWidth_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_phiWidth"+val_2+binning[i]).c_str(), h_fakeSC_phiWidth_new[i]);h_fakeSC_phiWidth_new[i]->SetName(("_new"+binning[i]).c_str());

        hist_infile_SCVal->GetObject(("h_goodSC_eta"+val_2+binning[i]).c_str(), h_goodSC_eta_new[i]);h_goodSC_eta_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_phi"+val_2+binning[i]).c_str(), h_goodSC_phi_new[i]);h_goodSC_phi_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_et"+val_2+binning[i]).c_str(), h_goodSC_et_new[i]);h_goodSC_et_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodSC_kinenergy"+val_2+binning[i]).c_str(), h_goodSC_energy_new[i]);h_goodSC_energy_new[i]->SetName(("_new"+binning[i]).c_str());

        hist_infile_SCVal->GetObject(("h_fakeSC_eta"+val_2+binning[i]).c_str(), h_fakeSC_eta_new[i]);h_fakeSC_eta_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_phi"+val_2+binning[i]).c_str(), h_fakeSC_phi_new[i]);h_fakeSC_phi_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_et"+val_2+binning[i]).c_str(), h_fakeSC_et_new[i]);h_fakeSC_et_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_fakeSC_kinenergy"+val_2+binning[i]).c_str(), h_fakeSC_energy_new[i]);h_fakeSC_energy_new[i]->SetName(("_new"+binning[i]).c_str());

        hist_infile_SCVal->GetObject(("h_goodFakeSC_deta"+val_2+binning[i]).c_str(), h_goodFakeSC_deta_new[i]);h_goodFakeSC_deta_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodFakeSC_dphi"+val_2+binning[i]).c_str(), h_goodFakeSC_dphi_new[i]);h_goodFakeSC_dphi_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodFakeSC_detrel"+val_2+binning[i]).c_str(), h_goodFakeSC_detrel_new[i]);h_goodFakeSC_detrel_new[i]->SetName(("_new"+binning[i]).c_str());
        hist_infile_SCVal->GetObject(("h_goodFakeSC_denergyrel"+val_2+binning[i]).c_str(), h_goodFakeSC_denergyrel_new[i]);h_goodFakeSC_denergyrel_new[i]->SetName(("_new"+binning[i]).c_str());

        hist_infile_SCRef->GetObject(("h_goodSC_R9"+ref_2+binning[i]).c_str(), h_goodSC_R9_old[i]);h_goodSC_R9_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_SigmaIetaIeta"+ref_2+binning[i]).c_str(), h_goodSC_SigmaIetaIeta_old[i]);h_goodSC_SigmaIetaIeta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_SigmaIphiIphi"+ref_2+binning[i]).c_str(), h_goodSC_SigmaIphiIphi_old[i]);h_goodSC_SigmaIphiIphi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_fullR9"+ref_2+binning[i]).c_str(), h_goodSC_fullR9_old[i]);h_goodSC_fullR9_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_fullSigmaIetaIeta"+ref_2+binning[i]).c_str(), h_goodSC_fullSigmaIetaIeta_old[i]);h_goodSC_fullSigmaIetaIeta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_fullSigmaIphiIphi"+ref_2+binning[i]).c_str(), h_goodSC_fullSigmaIphiIphi_old[i]);h_goodSC_fullSigmaIphiIphi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_etaWidth"+ref_2+binning[i]).c_str(), h_goodSC_etaWidth_old[i]);h_goodSC_etaWidth_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_phiWidth"+ref_2+binning[i]).c_str(), h_goodSC_phiWidth_old[i]);h_goodSC_phiWidth_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_R9"+ref_2+binning[i]).c_str(), h_fakeSC_R9_old[i]);h_fakeSC_R9_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_SigmaIetaIeta"+ref_2+binning[i]).c_str(), h_fakeSC_SigmaIetaIeta_old[i]);h_fakeSC_SigmaIetaIeta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_SigmaIphiIphi"+ref_2+binning[i]).c_str(), h_fakeSC_SigmaIphiIphi_old[i]);h_fakeSC_SigmaIphiIphi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_fullR9"+ref_2+binning[i]).c_str(), h_fakeSC_fullR9_old[i]);h_fakeSC_fullR9_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_fullSigmaIetaIeta"+ref_2+binning[i]).c_str(), h_fakeSC_fullSigmaIetaIeta_old[i]);h_fakeSC_fullSigmaIetaIeta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_fullSigmaIphiIphi"+ref_2+binning[i]).c_str(), h_fakeSC_fullSigmaIphiIphi_old[i]);h_fakeSC_fullSigmaIphiIphi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_etaWidth"+ref_2+binning[i]).c_str(), h_fakeSC_etaWidth_old[i]);h_fakeSC_etaWidth_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_phiWidth"+ref_2+binning[i]).c_str(), h_fakeSC_phiWidth_old[i]);h_fakeSC_phiWidth_old[i]->SetName(("_old"+binning[i]).c_str());

        hist_infile_SCRef->GetObject(("h_goodSC_eta"+ref_2+binning[i]).c_str(), h_goodSC_eta_old[i]);h_goodSC_eta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_phi"+ref_2+binning[i]).c_str(), h_goodSC_phi_old[i]);h_goodSC_phi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_et"+ref_2+binning[i]).c_str(), h_goodSC_et_old[i]);h_goodSC_et_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodSC_kinenergy"+ref_2+binning[i]).c_str(), h_goodSC_energy_old[i]);h_goodSC_energy_old[i]->SetName(("_old"+binning[i]).c_str());

        hist_infile_SCRef->GetObject(("h_fakeSC_eta"+ref_2+binning[i]).c_str(), h_fakeSC_eta_old[i]);h_fakeSC_eta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_phi"+ref_2+binning[i]).c_str(), h_fakeSC_phi_old[i]);h_fakeSC_phi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_et"+ref_2+binning[i]).c_str(), h_fakeSC_et_old[i]);h_fakeSC_et_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_fakeSC_kinenergy"+ref_2+binning[i]).c_str(), h_fakeSC_energy_old[i]);h_fakeSC_energy_old[i]->SetName(("_old"+binning[i]).c_str());

        hist_infile_SCRef->GetObject(("h_goodFakeSC_deta"+ref_2+binning[i]).c_str(), h_goodFakeSC_deta_old[i]);h_goodFakeSC_deta_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodFakeSC_dphi"+ref_2+binning[i]).c_str(), h_goodFakeSC_dphi_old[i]);h_goodFakeSC_dphi_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodFakeSC_detrel"+ref_2+binning[i]).c_str(), h_goodFakeSC_detrel_old[i]);h_goodFakeSC_detrel_old[i]->SetName(("_old"+binning[i]).c_str());
        hist_infile_SCRef->GetObject(("h_goodFakeSC_denergyrel"+ref_2+binning[i]).c_str(), h_goodFakeSC_denergyrel_old[i]);h_goodFakeSC_denergyrel_old[i]->SetName(("_old"+binning[i]).c_str());
    }

    cout<<"Closed input file"<<endl;
}

void ReadInfiles_Integrated(string SuperClusterRef, string SuperClusterVal)
{
    cout << "Reading input file" << endl;

    //1: Pt20to40
    //2: Pt20toInf
    //3: Pt40toInf

    TFile *hist_infile_SCRef[3];
    TFile *hist_infile_SCVal[3];

    string ref_1 = "", ref_2 = "", val_1 = "", val_2 = "";
    if (SuperClusterRef == "SingleSet")
    {
        hist_infile_SCRef[0] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20to40_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCRef[1] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCRef[2] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-40toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        ref_1 = "DeepSC";
        ref_2 = "_new";
    }
    else if (SuperClusterRef == "Mustache")
    {
        hist_infile_SCRef[0] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20to40_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCRef[1] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCRef[2] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-40toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        ref_1 = "SC";
        ref_2 = "_old";
    }
    if (SuperClusterVal == "SingleSet")
    {
        hist_infile_SCVal[0] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20to40_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCVal[1] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        hist_infile_SCVal[2] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-40toInf_OptPFRHT_OptMustSingle_ETWeight.root");
        val_1 = "DeepSC";
        val_2 = "_new";
    }
    else if (SuperClusterVal == "LocalSets")
    {
        hist_infile_SCVal[0] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20to40_OptPFRHT_OptMustLocalParams.root");
        hist_infile_SCVal[1] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-20toInf_OptPFRHT_OptMustLocalParams.root");
        hist_infile_SCVal[2] = TFile::Open("plotFiles/GJetValidation_Plots_Pt-40toInf_OptPFRHT_OptMustLocalParams.root");
        val_1 = "DeepSC";
        val_2 = "_new";
    }

    vector<string> ptBins = {"_pt20to40","_pt20toInf","_pt40toInf"};

    for (int ptIdx = 0; ptIdx < 3; ptIdx++)
    {
        for (int i = 0; i < 3; i++)
        {
            cout<< ptIdx<<"\t"<<i<<endl;

            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_R9" + ref_2 + binning[i]).c_str(), h_goodSC_R9_old_integrated[ptIdx][i]); h_goodSC_R9_old_integrated[ptIdx][i]->SetName(("h_goodSC_R9_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_SigmaIetaIeta" + ref_2 + binning[i]).c_str(), h_goodSC_SigmaIetaIeta_old_integrated[ptIdx][i]);h_goodSC_SigmaIetaIeta_old_integrated[ptIdx][i]->SetName(("h_goodSC_SigmaIetaIeta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_SigmaIphiIphi" + ref_2 + binning[i]).c_str(), h_goodSC_SigmaIphiIphi_old_integrated[ptIdx][i]);h_goodSC_SigmaIphiIphi_old_integrated[ptIdx][i]->SetName(("h_goodSC_SigmaIphiIphi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_fullR9" + ref_2 + binning[i]).c_str(), h_goodSC_fullR9_old_integrated[ptIdx][i]);h_goodSC_fullR9_old_integrated[ptIdx][i]->SetName(("h_goodSC_fullR9_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_fullSigmaIetaIeta" + ref_2 + binning[i]).c_str(), h_goodSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]);h_goodSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]->SetName(("h_goodSC_fullSigmaIetaIeta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_fullSigmaIphiIphi" + ref_2 + binning[i]).c_str(), h_goodSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]);h_goodSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]->SetName(("h_goodSC_fullSigmaIphiIphi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_etaWidth" + ref_2 + binning[i]).c_str(), h_goodSC_etaWidth_old_integrated[ptIdx][i]);h_goodSC_etaWidth_old_integrated[ptIdx][i]->SetName(("h_goodSC_etaWidth_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_phiWidth" + ref_2 + binning[i]).c_str(), h_goodSC_phiWidth_old_integrated[ptIdx][i]);h_goodSC_phiWidth_old_integrated[ptIdx][i]->SetName(("h_goodSC_phiWidth_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_R9" + ref_2 + binning[i]).c_str(), h_fakeSC_R9_old_integrated[ptIdx][i]);h_fakeSC_R9_old_integrated[ptIdx][i]->SetName(("h_fakeSC_R9_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_SigmaIetaIeta" + ref_2 + binning[i]).c_str(), h_fakeSC_SigmaIetaIeta_old_integrated[ptIdx][i]);h_fakeSC_SigmaIetaIeta_old_integrated[ptIdx][i]->SetName(("h_fakeSC_SigmaIetaIeta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_SigmaIphiIphi" + ref_2 + binning[i]).c_str(), h_fakeSC_SigmaIphiIphi_old_integrated[ptIdx][i]);h_fakeSC_SigmaIphiIphi_old_integrated[ptIdx][i]->SetName(("h_fakeSC_SigmaIphiIphi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_fullR9" + ref_2 + binning[i]).c_str(), h_fakeSC_fullR9_old_integrated[ptIdx][i]);h_fakeSC_fullR9_old_integrated[ptIdx][i]->SetName(("h_fakeSC_fullR9_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_fullSigmaIetaIeta" + ref_2 + binning[i]).c_str(), h_fakeSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]);h_fakeSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]->SetName(("h_fakeSC_fullSigmaIetaIeta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_fullSigmaIphiIphi" + ref_2 + binning[i]).c_str(), h_fakeSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]);h_fakeSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]->SetName(("h_fakeSC_fullSigmaIphiIphi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_etaWidth" + ref_2 + binning[i]).c_str(), h_fakeSC_etaWidth_old_integrated[ptIdx][i]);h_fakeSC_etaWidth_old_integrated[ptIdx][i]->SetName(("h_fakeSC_etaWidth_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_phiWidth" + ref_2 + binning[i]).c_str(), h_fakeSC_phiWidth_old_integrated[ptIdx][i]);h_fakeSC_phiWidth_old_integrated[ptIdx][i]->SetName(("h_fakeSC_phiWidth_old" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_eta" + ref_2 + binning[i]).c_str(), h_goodSC_eta_old_integrated[ptIdx][i]);h_goodSC_eta_old_integrated[ptIdx][i]->SetName(("h_goodSC_eta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_phi" + ref_2 + binning[i]).c_str(), h_goodSC_phi_old_integrated[ptIdx][i]);h_goodSC_phi_old_integrated[ptIdx][i]->SetName(("h_goodSC_phi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_et" + ref_2 + binning[i]).c_str(), h_goodSC_et_old_integrated[ptIdx][i]);h_goodSC_et_old_integrated[ptIdx][i]->SetName(("h_goodSC_et_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodSC_kinenergy" + ref_2 + binning[i]).c_str(), h_goodSC_energy_old_integrated[ptIdx][i]);h_goodSC_energy_old_integrated[ptIdx][i]->SetName(("h_goodSC_energy_old" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_eta" + ref_2 + binning[i]).c_str(), h_fakeSC_eta_old_integrated[ptIdx][i]);h_fakeSC_eta_old_integrated[ptIdx][i]->SetName(("h_fakeSC_eta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_phi" + ref_2 + binning[i]).c_str(), h_fakeSC_phi_old_integrated[ptIdx][i]);h_fakeSC_phi_old_integrated[ptIdx][i]->SetName(("h_fakeSC_phi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_et" + ref_2 + binning[i]).c_str(), h_fakeSC_et_old_integrated[ptIdx][i]);h_fakeSC_et_old_integrated[ptIdx][i]->SetName(("h_fakeSC_et_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_fakeSC_kinenergy" + ref_2 + binning[i]).c_str(), h_fakeSC_energy_old_integrated[ptIdx][i]);h_fakeSC_energy_old_integrated[ptIdx][i]->SetName(("h_fakeSC_energy_old" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCRef[ptIdx]->GetObject(("h_goodFakeSC_deta" + ref_2 + binning[i]).c_str(), h_goodFakeSC_deta_old_integrated[ptIdx][i]);h_goodFakeSC_deta_old_integrated[ptIdx][i]->SetName(("h_goodFakeSC_deta_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodFakeSC_dphi" + ref_2 + binning[i]).c_str(), h_goodFakeSC_dphi_old_integrated[ptIdx][i]);h_goodFakeSC_dphi_old_integrated[ptIdx][i]->SetName(("h_goodFakeSC_dphi_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodFakeSC_detrel" + ref_2 + binning[i]).c_str(), h_goodFakeSC_detrel_old_integrated[ptIdx][i]);h_goodFakeSC_detrel_old_integrated[ptIdx][i]->SetName(("h_goodFakeSC_detrel_old" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCRef[ptIdx]->GetObject(("h_goodFakeSC_denergyrel" + ref_2 + binning[i]).c_str(), h_goodFakeSC_denergyrel_old_integrated[ptIdx][i]);h_goodFakeSC_denergyrel_old_integrated[ptIdx][i]->SetName(("h_goodFakeSC_denergyrel_old" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_R9" + val_2 + binning[i]).c_str(), h_goodSC_R9_new_integrated[ptIdx][i]); h_goodSC_R9_new_integrated[ptIdx][i] ->SetName(("h_goodSC_R9_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_SigmaIetaIeta" + val_2 + binning[i]).c_str(), h_goodSC_SigmaIetaIeta_new_integrated[ptIdx][i]); h_goodSC_SigmaIetaIeta_new_integrated[ptIdx][i]->SetName(("h_goodSC_SigmaIetaIeta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_SigmaIphiIphi" + val_2 + binning[i]).c_str(), h_goodSC_SigmaIphiIphi_new_integrated[ptIdx][i]); h_goodSC_SigmaIphiIphi_new_integrated[ptIdx][i]->SetName(("h_goodSC_SigmaIphiIphi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_fullR9" + val_2 + binning[i]).c_str(), h_goodSC_fullR9_new_integrated[ptIdx][i]); h_goodSC_fullR9_new_integrated[ptIdx][i]->SetName(("h_goodSC_fullR9_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_fullSigmaIetaIeta" + val_2 + binning[i]).c_str(), h_goodSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]); h_goodSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]->SetName(("h_goodSC_fullSigmaIetaIeta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_fullSigmaIphiIphi" + val_2 + binning[i]).c_str(), h_goodSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]); h_goodSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]->SetName(("h_goodSC_fullSigmaIphiIphi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_etaWidth" + val_2 + binning[i]).c_str(), h_goodSC_etaWidth_new_integrated[ptIdx][i]); h_goodSC_etaWidth_new_integrated[ptIdx][i]->SetName(("h_goodSC_etaWidth_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_phiWidth" + val_2 + binning[i]).c_str(), h_goodSC_phiWidth_new_integrated[ptIdx][i]); h_goodSC_phiWidth_new_integrated[ptIdx][i]->SetName(("h_goodSC_phiWidth_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_R9" + val_2 + binning[i]).c_str(), h_fakeSC_R9_new_integrated[ptIdx][i]); h_fakeSC_R9_new_integrated[ptIdx][i]->SetName(("h_fakeSC_R9_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_SigmaIetaIeta" + val_2 + binning[i]).c_str(), h_fakeSC_SigmaIetaIeta_new_integrated[ptIdx][i]); h_fakeSC_SigmaIetaIeta_new_integrated[ptIdx][i]->SetName(("h_fakeSC_SigmaIetaIeta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_SigmaIphiIphi" + val_2 + binning[i]).c_str(), h_fakeSC_SigmaIphiIphi_new_integrated[ptIdx][i]); h_fakeSC_SigmaIphiIphi_new_integrated[ptIdx][i]->SetName(("h_fakeSC_SigmaIphiIphi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_fullR9" + val_2 + binning[i]).c_str(), h_fakeSC_fullR9_new_integrated[ptIdx][i]); h_fakeSC_fullR9_new_integrated[ptIdx][i]->SetName(("h_fakeSC_fullR9_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_fullSigmaIetaIeta" + val_2 + binning[i]).c_str(), h_fakeSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]); h_fakeSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]->SetName(("h_fakeSC_fullSigmaIetaIeta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_fullSigmaIphiIphi" + val_2 + binning[i]).c_str(), h_fakeSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]); h_fakeSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]->SetName(("h_fakeSC_fullSigmaIphiIphi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_etaWidth" + val_2 + binning[i]).c_str(), h_fakeSC_etaWidth_new_integrated[ptIdx][i]); h_fakeSC_etaWidth_new_integrated[ptIdx][i]->SetName(("h_fakeSC_etaWidth_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_phiWidth" + val_2 + binning[i]).c_str(), h_fakeSC_phiWidth_new_integrated[ptIdx][i]); h_fakeSC_phiWidth_new_integrated[ptIdx][i]->SetName(("h_fakeSC_phiWidth_new" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_eta" + val_2 + binning[i]).c_str(), h_goodSC_eta_new_integrated[ptIdx][i]); h_goodSC_eta_new_integrated[ptIdx][i]->SetName(("h_goodSC_eta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_phi" + val_2 + binning[i]).c_str(), h_goodSC_phi_new_integrated[ptIdx][i]); h_goodSC_phi_new_integrated[ptIdx][i]->SetName(("h_goodSC_phi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_et" + val_2 + binning[i]).c_str(), h_goodSC_et_new_integrated[ptIdx][i]); h_goodSC_et_new_integrated[ptIdx][i]->SetName(("h_goodSC_et_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodSC_kinenergy" + val_2 + binning[i]).c_str(), h_goodSC_energy_new_integrated[ptIdx][i]); h_goodSC_energy_new_integrated[ptIdx][i]->SetName(("h_goodSC_energy_new" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_eta" + val_2 + binning[i]).c_str(), h_fakeSC_eta_new_integrated[ptIdx][i]); h_fakeSC_eta_new_integrated[ptIdx][i]->SetName(("h_fakeSC_eta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_phi" + val_2 + binning[i]).c_str(), h_fakeSC_phi_new_integrated[ptIdx][i]); h_fakeSC_phi_new_integrated[ptIdx][i]->SetName(("h_fakeSC_phi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_et" + val_2 + binning[i]).c_str(), h_fakeSC_et_new_integrated[ptIdx][i]);h_fakeSC_et_new_integrated[ptIdx][i]->SetName(("h_fakeSC_et_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_fakeSC_kinenergy" + val_2 + binning[i]).c_str(), h_fakeSC_energy_new_integrated[ptIdx][i]);h_fakeSC_energy_new_integrated[ptIdx][i]->SetName(("h_fakeSC_energy_new" + binning[i]+ptBins[ptIdx]).c_str());

            hist_infile_SCVal[ptIdx]->GetObject(("h_goodFakeSC_deta" + val_2 + binning[i]).c_str(), h_goodFakeSC_deta_new_integrated[ptIdx][i]);h_goodFakeSC_deta_new_integrated[ptIdx][i]->SetName(("h_goodFakeSC_deta_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodFakeSC_dphi" + val_2 + binning[i]).c_str(), h_goodFakeSC_dphi_new_integrated[ptIdx][i]);h_goodFakeSC_dphi_new_integrated[ptIdx][i]->SetName(("h_goodFakeSC_dphi_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodFakeSC_detrel" + val_2 + binning[i]).c_str(), h_goodFakeSC_detrel_new_integrated[ptIdx][i]);h_goodFakeSC_detrel_new_integrated[ptIdx][i]->SetName(("h_goodFakeSC_detrel_new" + binning[i]+ptBins[ptIdx]).c_str());
            hist_infile_SCVal[ptIdx]->GetObject(("h_goodFakeSC_denergyrel" + val_2 + binning[i]).c_str(), h_goodFakeSC_denergyrel_new_integrated[ptIdx][i]);h_goodFakeSC_denergyrel_new_integrated[ptIdx][i]->SetName(("h_goodFakeSC_denergyrel_new" + binning[i]+ptBins[ptIdx]).c_str());

        }
    }

    double weight_pt20to40 = 232.9 * 5000.0 / 505000, 
          weight_pt20toInf = 3186.0 * 5000.0 / 506800, 
          weight_pt40toInf = 878.1 * 5000.0 / 502100;
    vector<double> weights = {weight_pt20to40,weight_pt20toInf,weight_pt40toInf};

    //Add files together - weighted
    for (int i = 0; i < 3; i++)
    {
        for (int ptIdx = 0; ptIdx < 3; ptIdx++)
        {

            h_goodSC_R9_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_SigmaIetaIeta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_SigmaIphiIphi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullR9_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_etaWidth_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_phiWidth_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_R9_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_SigmaIetaIeta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_SigmaIphiIphi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullR9_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullSigmaIetaIeta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullSigmaIphiIphi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_etaWidth_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_phiWidth_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_goodSC_eta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_phi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_et_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_energy_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_fakeSC_eta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_phi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_et_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_energy_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_goodFakeSC_deta_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_dphi_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_detrel_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_denergyrel_new_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_goodSC_R9_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_SigmaIetaIeta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_SigmaIphiIphi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullR9_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_etaWidth_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_phiWidth_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_R9_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_SigmaIetaIeta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_SigmaIphiIphi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullR9_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullSigmaIetaIeta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_fullSigmaIphiIphi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_etaWidth_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_phiWidth_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_goodSC_eta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_phi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_et_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodSC_energy_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_fakeSC_eta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_phi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_et_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_fakeSC_energy_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);

            h_goodFakeSC_deta_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_dphi_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_detrel_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
            h_goodFakeSC_denergyrel_old_integrated[ptIdx][i]->Scale(weights[ptIdx]);
        }

            h_goodSC_R9_new[i] = (TH1F*)h_goodSC_R9_new_integrated[0][i]->Clone(); h_goodSC_R9_new[i] ->Add(h_goodSC_R9_new_integrated[1][i]); h_goodSC_R9_new[i] ->Add(h_goodSC_R9_new_integrated[2][i]);
            h_goodSC_SigmaIetaIeta_new[i] = (TH1F*)h_goodSC_SigmaIetaIeta_new_integrated[0][i]->Clone();  h_goodSC_SigmaIetaIeta_new[i]->Add(h_goodSC_SigmaIetaIeta_new_integrated[1][i]); h_goodSC_SigmaIetaIeta_new[i]->Add(h_goodSC_SigmaIetaIeta_new_integrated[2][i]);
            h_goodSC_SigmaIphiIphi_new[i] = (TH1F*)h_goodSC_SigmaIphiIphi_new_integrated[0][i]->Clone();  h_goodSC_SigmaIphiIphi_new[i]->Add(h_goodSC_SigmaIphiIphi_new_integrated[1][i]); h_goodSC_SigmaIphiIphi_new[i]->Add(h_goodSC_SigmaIphiIphi_new_integrated[2][i]);
            h_goodSC_fullR9_new[i] = (TH1F*)h_goodSC_fullR9_new_integrated[0][i]->Clone();  h_goodSC_fullR9_new[i]->Add(h_goodSC_fullR9_new_integrated[1][i]); h_goodSC_fullR9_new[i]->Add(h_goodSC_fullR9_new_integrated[2][i]);
            h_goodSC_fullSigmaIetaIeta_new[i] = (TH1F*)h_goodSC_fullSigmaIetaIeta_new_integrated[0][i]->Clone();  h_goodSC_fullSigmaIetaIeta_new[i]->Add(h_goodSC_fullSigmaIetaIeta_new_integrated[1][i]); h_goodSC_fullSigmaIetaIeta_new[i]->Add(h_goodSC_fullSigmaIetaIeta_new_integrated[2][i]);
            h_goodSC_fullSigmaIphiIphi_new[i] = (TH1F*)h_goodSC_fullSigmaIphiIphi_new_integrated[0][i]->Clone();  h_goodSC_fullSigmaIphiIphi_new[i]->Add(h_goodSC_fullSigmaIphiIphi_new_integrated[1][i]); h_goodSC_fullSigmaIphiIphi_new[i]->Add(h_goodSC_fullSigmaIphiIphi_new_integrated[2][i]);
            h_goodSC_etaWidth_new[i] = (TH1F*)h_goodSC_etaWidth_new_integrated[0][i]->Clone();  h_goodSC_etaWidth_new[i]->Add(h_goodSC_etaWidth_new_integrated[1][i]); h_goodSC_etaWidth_new[i]->Add(h_goodSC_etaWidth_new_integrated[2][i]);
            h_goodSC_phiWidth_new[i] = (TH1F*)h_goodSC_phiWidth_new_integrated[0][i]->Clone();  h_goodSC_phiWidth_new[i]->Add(h_goodSC_phiWidth_new_integrated[1][i]); h_goodSC_phiWidth_new[i]->Add(h_goodSC_phiWidth_new_integrated[2][i]);
            h_fakeSC_R9_new[i] = (TH1F*)h_fakeSC_R9_new_integrated[0][i]->Clone();  h_fakeSC_R9_new[i]->Add(h_fakeSC_R9_new_integrated[1][i]); h_fakeSC_R9_new[i]->Add(h_fakeSC_R9_new_integrated[2][i]);
            h_fakeSC_SigmaIetaIeta_new[i] = (TH1F*)h_fakeSC_SigmaIetaIeta_new_integrated[0][i]->Clone();  h_fakeSC_SigmaIetaIeta_new[i]->Add(h_fakeSC_SigmaIetaIeta_new_integrated[1][i]); h_fakeSC_SigmaIetaIeta_new[i]->Add(h_fakeSC_SigmaIetaIeta_new_integrated[2][i]);
            h_fakeSC_SigmaIphiIphi_new[i] = (TH1F*)h_fakeSC_SigmaIphiIphi_new_integrated[0][i]->Clone();  h_fakeSC_SigmaIphiIphi_new[i]->Add(h_fakeSC_SigmaIphiIphi_new_integrated[1][i]); h_fakeSC_SigmaIphiIphi_new[i]->Add(h_fakeSC_SigmaIphiIphi_new_integrated[2][i]);
            h_fakeSC_fullR9_new[i] = (TH1F*)h_fakeSC_fullR9_new_integrated[0][i]->Clone();  h_fakeSC_fullR9_new[i]->Add(h_fakeSC_fullR9_new_integrated[1][i]); h_fakeSC_fullR9_new[i]->Add(h_fakeSC_fullR9_new_integrated[2][i]);
            h_fakeSC_fullSigmaIetaIeta_new[i] = (TH1F*)h_fakeSC_fullSigmaIetaIeta_new_integrated[0][i]->Clone();  h_fakeSC_fullSigmaIetaIeta_new[i]->Add(h_fakeSC_fullSigmaIetaIeta_new_integrated[1][i]); h_fakeSC_fullSigmaIetaIeta_new[i]->Add(h_fakeSC_fullSigmaIetaIeta_new_integrated[2][i]);
            h_fakeSC_fullSigmaIphiIphi_new[i] = (TH1F*)h_fakeSC_fullSigmaIphiIphi_new_integrated[0][i]->Clone();  h_fakeSC_fullSigmaIphiIphi_new[i]->Add(h_fakeSC_fullSigmaIphiIphi_new_integrated[1][i]); h_fakeSC_fullSigmaIphiIphi_new[i]->Add(h_fakeSC_fullSigmaIphiIphi_new_integrated[2][i]);
            h_fakeSC_etaWidth_new[i] = (TH1F*)h_fakeSC_etaWidth_new_integrated[0][i]->Clone();  h_fakeSC_etaWidth_new[i]->Add(h_fakeSC_etaWidth_new_integrated[1][i]); h_fakeSC_etaWidth_new[i]->Add(h_fakeSC_etaWidth_new_integrated[2][i]);
            h_fakeSC_phiWidth_new[i] = (TH1F*)h_fakeSC_phiWidth_new_integrated[0][i]->Clone();  h_fakeSC_phiWidth_new[i]->Add(h_fakeSC_phiWidth_new_integrated[1][i]); h_fakeSC_phiWidth_new[i]->Add(h_fakeSC_phiWidth_new_integrated[2][i]);

            h_goodSC_eta_new[i] = (TH1F*)h_goodSC_eta_new_integrated[0][i]->Clone();  h_goodSC_eta_new[i]->Add(h_goodSC_eta_new_integrated[1][i]); h_goodSC_eta_new[i]->Add(h_goodSC_eta_new_integrated[2][i]);
            h_goodSC_phi_new[i] = (TH1F*)h_goodSC_phi_new_integrated[0][i]->Clone();   h_goodSC_phi_new[i]->Add(h_goodSC_phi_new_integrated[1][i]);  h_goodSC_phi_new[i]->Add(h_goodSC_phi_new_integrated[2][i]);
            h_goodSC_et_new[i] = (TH1F*)h_goodSC_et_new_integrated[0][i]->Clone();  h_goodSC_et_new[i]->Add(h_goodSC_et_new_integrated[1][i]); h_goodSC_et_new[i]->Add(h_goodSC_et_new_integrated[2][i]);
            h_goodSC_energy_new[i] = (TH1F*)h_goodSC_energy_new_integrated[0][i]->Clone();  h_goodSC_energy_new[i]->Add(h_goodSC_energy_new_integrated[1][i]); h_goodSC_energy_new[i]->Add(h_goodSC_energy_new_integrated[2][i]);

            h_fakeSC_eta_new[i] = (TH1F*)h_fakeSC_eta_new_integrated[0][i]->Clone();  h_fakeSC_eta_new[i]->Add(h_fakeSC_eta_new_integrated[1][i]); h_fakeSC_eta_new[i]->Add(h_fakeSC_eta_new_integrated[2][i]);
            h_fakeSC_phi_new[i] = (TH1F*)h_fakeSC_phi_new_integrated[0][i]->Clone();  h_fakeSC_phi_new[i]->Add(h_fakeSC_phi_new_integrated[1][i]); h_fakeSC_phi_new[i]->Add(h_fakeSC_phi_new_integrated[2][i]);
            h_fakeSC_et_new[i] = (TH1F*)h_fakeSC_et_new_integrated[0][i]->Clone();  h_fakeSC_et_new[i]->Add(h_fakeSC_et_new_integrated[1][i]); h_fakeSC_et_new[i]->Add(h_fakeSC_et_new_integrated[2][i]);
            h_fakeSC_energy_new[i] = (TH1F*)h_fakeSC_energy_new_integrated[0][i]->Clone();  h_fakeSC_energy_new[i]->Add(h_fakeSC_energy_new_integrated[1][i]); h_fakeSC_energy_new[i]->Add(h_fakeSC_energy_new_integrated[2][i]);

            h_goodFakeSC_deta_new[i] = (TH1F*)h_goodFakeSC_deta_new_integrated[0][i]->Clone();  h_goodFakeSC_deta_new[i]->Add(h_goodFakeSC_deta_new_integrated[1][i]); h_goodFakeSC_deta_new[i]->Add(h_goodFakeSC_deta_new_integrated[2][i]);
            h_goodFakeSC_dphi_new[i] = (TH1F*)h_goodFakeSC_deta_new_integrated[0][i]->Clone();  h_goodFakeSC_dphi_new[i]->Add(h_goodFakeSC_deta_new_integrated[1][i]); h_goodFakeSC_dphi_new[i]->Add(h_goodFakeSC_deta_new_integrated[2][i]);
            h_goodFakeSC_detrel_new[i] = (TH1F*)h_goodFakeSC_deta_new_integrated[0][i]->Clone();  h_goodFakeSC_detrel_new[i]->Add(h_goodFakeSC_deta_new_integrated[1][i]); h_goodFakeSC_detrel_new[i]->Add(h_goodFakeSC_deta_new_integrated[2][i]);
            h_goodFakeSC_denergyrel_new[i] = (TH1F*)h_goodFakeSC_deta_new_integrated[0][i]->Clone();  h_goodFakeSC_denergyrel_new[i]->Add(h_goodFakeSC_deta_new_integrated[1][i]); h_goodFakeSC_denergyrel_new[i]->Add(h_goodFakeSC_deta_new_integrated[2][i]);

            h_goodSC_R9_old[i] = (TH1F*)h_goodSC_R9_old_integrated[0][i]->Clone(); h_goodSC_R9_old[i] ->Add(h_goodSC_R9_old_integrated[1][i]); h_goodSC_R9_old[i] ->Add(h_goodSC_R9_old_integrated[2][i]);
            h_goodSC_SigmaIetaIeta_old[i] = (TH1F*)h_goodSC_SigmaIetaIeta_old_integrated[0][i]->Clone();  h_goodSC_SigmaIetaIeta_old[i]->Add(h_goodSC_SigmaIetaIeta_old_integrated[1][i]); h_goodSC_SigmaIetaIeta_old[i]->Add(h_goodSC_SigmaIetaIeta_old_integrated[2][i]);
            h_goodSC_SigmaIphiIphi_old[i] = (TH1F*)h_goodSC_SigmaIphiIphi_old_integrated[0][i]->Clone();  h_goodSC_SigmaIphiIphi_old[i]->Add(h_goodSC_SigmaIphiIphi_old_integrated[1][i]); h_goodSC_SigmaIphiIphi_old[i]->Add(h_goodSC_SigmaIphiIphi_old_integrated[2][i]);
            h_goodSC_fullR9_old[i] = (TH1F*)h_goodSC_fullR9_old_integrated[0][i]->Clone();  h_goodSC_fullR9_old[i]->Add(h_goodSC_fullR9_old_integrated[1][i]); h_goodSC_fullR9_old[i]->Add(h_goodSC_fullR9_old_integrated[2][i]);
            h_goodSC_fullSigmaIetaIeta_old[i] = (TH1F*)h_goodSC_fullSigmaIetaIeta_old_integrated[0][i]->Clone();  h_goodSC_fullSigmaIetaIeta_old[i]->Add(h_goodSC_fullSigmaIetaIeta_old_integrated[1][i]); h_goodSC_fullSigmaIetaIeta_old[i]->Add(h_goodSC_fullSigmaIetaIeta_old_integrated[2][i]);
            h_goodSC_fullSigmaIphiIphi_old[i] = (TH1F*)h_goodSC_fullSigmaIphiIphi_old_integrated[0][i]->Clone();  h_goodSC_fullSigmaIphiIphi_old[i]->Add(h_goodSC_fullSigmaIphiIphi_old_integrated[1][i]); h_goodSC_fullSigmaIphiIphi_old[i]->Add(h_goodSC_fullSigmaIphiIphi_old_integrated[2][i]);
            h_goodSC_etaWidth_old[i] = (TH1F*)h_goodSC_etaWidth_old_integrated[0][i]->Clone();  h_goodSC_etaWidth_old[i]->Add(h_goodSC_etaWidth_old_integrated[1][i]); h_goodSC_etaWidth_old[i]->Add(h_goodSC_etaWidth_old_integrated[2][i]);
            h_goodSC_phiWidth_old[i] = (TH1F*)h_goodSC_phiWidth_old_integrated[0][i]->Clone();  h_goodSC_phiWidth_old[i]->Add(h_goodSC_phiWidth_old_integrated[1][i]); h_goodSC_phiWidth_old[i]->Add(h_goodSC_phiWidth_old_integrated[2][i]);
            h_fakeSC_R9_old[i] = (TH1F*)h_fakeSC_R9_old_integrated[0][i]->Clone();  h_fakeSC_R9_old[i]->Add(h_fakeSC_R9_old_integrated[1][i]); h_fakeSC_R9_old[i]->Add(h_fakeSC_R9_old_integrated[2][i]);
            h_fakeSC_SigmaIetaIeta_old[i] = (TH1F*)h_fakeSC_SigmaIetaIeta_old_integrated[0][i]->Clone();  h_fakeSC_SigmaIetaIeta_old[i]->Add(h_fakeSC_SigmaIetaIeta_old_integrated[1][i]); h_fakeSC_SigmaIetaIeta_old[i]->Add(h_fakeSC_SigmaIetaIeta_old_integrated[2][i]);
            h_fakeSC_SigmaIphiIphi_old[i] = (TH1F*)h_fakeSC_SigmaIphiIphi_old_integrated[0][i]->Clone();  h_fakeSC_SigmaIphiIphi_old[i]->Add(h_fakeSC_SigmaIphiIphi_old_integrated[1][i]); h_fakeSC_SigmaIphiIphi_old[i]->Add(h_fakeSC_SigmaIphiIphi_old_integrated[2][i]);
            h_fakeSC_fullR9_old[i] = (TH1F*)h_fakeSC_fullR9_old_integrated[0][i]->Clone();  h_fakeSC_fullR9_old[i]->Add(h_fakeSC_fullR9_old_integrated[1][i]); h_fakeSC_fullR9_old[i]->Add(h_fakeSC_fullR9_old_integrated[2][i]);
            h_fakeSC_fullSigmaIetaIeta_old[i] = (TH1F*)h_fakeSC_fullSigmaIetaIeta_old_integrated[0][i]->Clone();  h_fakeSC_fullSigmaIetaIeta_old[i]->Add(h_fakeSC_fullSigmaIetaIeta_old_integrated[1][i]); h_fakeSC_fullSigmaIetaIeta_old[i]->Add(h_fakeSC_fullSigmaIetaIeta_old_integrated[2][i]);
            h_fakeSC_fullSigmaIphiIphi_old[i] = (TH1F*)h_fakeSC_fullSigmaIphiIphi_old_integrated[0][i]->Clone();  h_fakeSC_fullSigmaIphiIphi_old[i]->Add(h_fakeSC_fullSigmaIphiIphi_old_integrated[1][i]); h_fakeSC_fullSigmaIphiIphi_old[i]->Add(h_fakeSC_fullSigmaIphiIphi_old_integrated[2][i]);
            h_fakeSC_etaWidth_old[i] = (TH1F*)h_fakeSC_etaWidth_old_integrated[0][i]->Clone();  h_fakeSC_etaWidth_old[i]->Add(h_fakeSC_etaWidth_old_integrated[1][i]); h_fakeSC_etaWidth_old[i]->Add(h_fakeSC_etaWidth_old_integrated[2][i]);
            h_fakeSC_phiWidth_old[i] = (TH1F*)h_fakeSC_phiWidth_old_integrated[0][i]->Clone();  h_fakeSC_phiWidth_old[i]->Add(h_fakeSC_phiWidth_old_integrated[1][i]); h_fakeSC_phiWidth_old[i]->Add(h_fakeSC_phiWidth_old_integrated[2][i]);

            h_goodSC_eta_old[i] = (TH1F*)h_goodSC_eta_old_integrated[0][i]->Clone();  h_goodSC_eta_old[i]->Add(h_goodSC_eta_old_integrated[1][i]); h_goodSC_eta_old[i]->Add(h_goodSC_eta_old_integrated[2][i]);
            h_goodSC_phi_old[i] = (TH1F*)h_goodSC_phi_old_integrated[0][i]->Clone();   h_goodSC_phi_old[i]->Add(h_goodSC_phi_old_integrated[1][i]);  h_goodSC_phi_old[i]->Add(h_goodSC_phi_old_integrated[2][i]);
            h_goodSC_et_old[i] = (TH1F*)h_goodSC_et_old_integrated[0][i]->Clone();  h_goodSC_et_old[i]->Add(h_goodSC_et_old_integrated[1][i]); h_goodSC_et_old[i]->Add(h_goodSC_et_old_integrated[2][i]);
            h_goodSC_energy_old[i] = (TH1F*)h_goodSC_energy_old_integrated[0][i]->Clone();  h_goodSC_energy_old[i]->Add(h_goodSC_energy_old_integrated[1][i]); h_goodSC_energy_old[i]->Add(h_goodSC_energy_old_integrated[2][i]);

            h_fakeSC_eta_old[i] = (TH1F*)h_fakeSC_eta_old_integrated[0][i]->Clone();  h_fakeSC_eta_old[i]->Add(h_fakeSC_eta_old_integrated[1][i]); h_fakeSC_eta_old[i]->Add(h_fakeSC_eta_old_integrated[2][i]);
            h_fakeSC_phi_old[i] = (TH1F*)h_fakeSC_phi_old_integrated[0][i]->Clone();  h_fakeSC_phi_old[i]->Add(h_fakeSC_phi_old_integrated[1][i]); h_fakeSC_phi_old[i]->Add(h_fakeSC_phi_old_integrated[2][i]);
            h_fakeSC_et_old[i] = (TH1F*)h_fakeSC_et_old_integrated[0][i]->Clone();  h_fakeSC_et_old[i]->Add(h_fakeSC_et_old_integrated[1][i]); h_fakeSC_et_old[i]->Add(h_fakeSC_et_old_integrated[2][i]);
            h_fakeSC_energy_old[i] = (TH1F*)h_fakeSC_energy_old_integrated[0][i]->Clone();  h_fakeSC_energy_old[i]->Add(h_fakeSC_energy_old_integrated[1][i]); h_fakeSC_energy_old[i]->Add(h_fakeSC_energy_old_integrated[2][i]);

            h_goodFakeSC_deta_old[i] = (TH1F*)h_goodFakeSC_deta_old_integrated[0][i]->Clone();  h_goodFakeSC_deta_old[i]->Add(h_goodFakeSC_deta_old_integrated[1][i]); h_goodFakeSC_deta_old[i]->Add(h_goodFakeSC_deta_old_integrated[2][i]);
            h_goodFakeSC_dphi_old[i] = (TH1F*)h_goodFakeSC_deta_old_integrated[0][i]->Clone();  h_goodFakeSC_dphi_old[i]->Add(h_goodFakeSC_deta_old_integrated[1][i]); h_goodFakeSC_dphi_old[i]->Add(h_goodFakeSC_deta_old_integrated[2][i]);
            h_goodFakeSC_detrel_old[i] = (TH1F*)h_goodFakeSC_deta_old_integrated[0][i]->Clone();  h_goodFakeSC_detrel_old[i]->Add(h_goodFakeSC_deta_old_integrated[1][i]); h_goodFakeSC_detrel_old[i]->Add(h_goodFakeSC_deta_old_integrated[2][i]);
            h_goodFakeSC_denergyrel_old[i] = (TH1F*)h_goodFakeSC_deta_old_integrated[0][i]->Clone();  h_goodFakeSC_denergyrel_old[i]->Add(h_goodFakeSC_deta_old_integrated[1][i]); h_goodFakeSC_denergyrel_old[i]->Add(h_goodFakeSC_deta_old_integrated[2][i]);

    }

    cout << "Closed input file" << endl;
}

void drawHisto(TH1F* h_old, TH1F* h_new, std::string x_label, std::string drawType, std::string Name, bool log, bool fit=false, std::string fitFunc_="cruijff", std::string refLegend="Mustache", std::string valLegend="DeepSC", string outputDir_="")
{

   h_old->Scale(1./h_old->GetEntries());
   h_new->Scale(1./h_new->GetEntries());

   TFitResultPtr frp_old;
   TFitResultPtr frp_new;
   TF1* doubleCB_old = new TF1();
   TF1* doubleCB_new = new TF1();

   if(fit){
      TH1F* h_old_clone = (TH1F*)h_old->Clone();  
      TH1F* h_new_clone = (TH1F*)h_new->Clone();    

      //doubleCB_old = fitHisto(h_old_clone, fitFunc_);
      doubleCB_old->SetLineColor(kRed+1);

      //doubleCB_new = fitHisto(h_new_clone, fitFunc_);
      doubleCB_new->SetLineColor(kBlue+1);
   }

   //h_old->SetTitle("");
   //h_new->SetTitle("");

   std::vector<float> maxima;
   maxima.resize(4);
   maxima[0] = h_old->GetMaximum(); 
   maxima[1] = h_new->GetMaximum();
   if(fit) maxima[2] = doubleCB_old->GetMaximum();
   if(fit) maxima[3] = doubleCB_new->GetMaximum(); 
   std::sort(maxima.begin(),maxima.end());  
   
   h_old->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   h_old->SetLineColor(kRed+1);
   h_old->SetMarkerColor(kRed+1);
   h_old->SetLineWidth(2);
   h_old->GetXaxis()->SetTitle(x_label.c_str());

   h_new->SetLineColor(kBlue+1);
   h_new->SetMarkerColor(kBlue+1);
   h_new->SetLineWidth(2);

   /*float maximum = -1.;
   if( h_old -> GetMaximum() > h_new -> GetMaximum()) maximum=h_old -> GetMaximum();
   else maximum=h_new -> GetMaximum();  
   h_old -> SetMaximum( 1.1*maximum );*/
    
   TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);

   TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();
   
   TCanvas* c = new TCanvas();
 
   /*TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   if(log) cUp->SetLogy();
   cDown->Draw();
     
   cUp->cd();*/
   if(!fit){
      TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.04); 
      legend -> AddEntry(h_old,refLegend.c_str(),"L");
      legend -> AddEntry(h_new,valLegend.c_str(),"L");
      h_old->Draw(drawType.c_str());
      gPad -> Update();
      st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
      st_old->SetX1NDC(0.82); //new x start position
      st_old->SetX2NDC(0.99); //new x end position
      st_old->SetY1NDC(0.82); //new y start position
      st_old->SetY2NDC(0.94); //new y end position
      st_old->SetTextColor(kRed+1);
      st_old->Draw("sames");
      h_new->Draw(std::string(drawType+",sames").c_str());
      gPad -> Update();
      st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
      st_new->SetX1NDC(0.82); //new x start position
      st_new->SetX2NDC(0.99); //new x end position
      st_new->SetY1NDC(0.68); //new y start position
      st_new->SetY2NDC(0.80); //new y end position
      st_new->SetTextColor(kBlue+1);
      st_new->Draw("sames");
      legend -> Draw("same");
   }else{
      TLegend* legend = new TLegend(0.57, 0.77, 0.72, 0.89);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.05);
      legend -> AddEntry(h_old,std::string(refLegend+": "+to_string(doubleCB_old->GetParameter(0))+" +/- "+to_string((doubleCB_old->GetParameter(2)+doubleCB_old->GetParameter(1))/2.)).c_str(),"L");
      legend -> AddEntry(h_new,std::string(valLegend+": "+to_string(doubleCB_new->GetParameter(0))+" +/- "+to_string((doubleCB_new->GetParameter(2)+doubleCB_new->GetParameter(1))/2.)).c_str(),"L");
      h_old->Draw(std::string(drawType).c_str());    
      h_new->Draw(std::string(drawType+",sames").c_str());
      doubleCB_old->Draw("L,same");
      doubleCB_new->Draw("L,same");
      legend -> Draw("same");   
      gPad -> Update();
      st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
      st_old->SetX1NDC(0.); //new x start position
      st_old->SetX2NDC(0.); //new x end position
      st_old->SetY1NDC(0.); //new y start position
      st_old->SetY2NDC(0.); //new y end position
      st_old->SetTextColor(kBlack);
      st_old->Draw("sames");
      gPad -> Update();
      st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
      st_new->SetX1NDC(0.); //new x start position
      st_new->SetX2NDC(0.); //new x end position
      st_new->SetY1NDC(0.); //new y start position
      st_new->SetY2NDC(0.); //new y end position
      st_new->SetTextColor(kBlack);
      st_new->Draw("sames");
   }
   //cDown->cd();
    /*
   TH1F* histo_ratio=(TH1F*)h_new->Clone("histo_ratio");
   if(histo_ratio->GetSumw2N()<=0) histo_ratio->Sumw2();
   histo_ratio->Divide(h_old);
    
   histo_ratio -> GetXaxis() -> SetTitle(x_label.c_str());
   histo_ratio -> GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
   histo_ratio -> SetMaximum(1.2);
   histo_ratio -> SetMinimum(0.8);
   histo_ratio -> SetMarkerColor(kBlack);
   histo_ratio -> SetMarkerSize(0.5);
   histo_ratio ->SetTitle("");
   histo_ratio -> GetXaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetXaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleOffset(0.7);
   histo_ratio -> Draw("e");
   TF1* f_const = new TF1("f_1", "[0]",histo_ratio->GetBinCenter(1)-histo_ratio->GetBinWidth(1)/2, histo_ratio->GetBinCenter(histo_ratio->GetNbinsX())+histo_ratio->GetBinWidth(histo_ratio->GetNbinsX())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   f_const -> SetLineWidth(2);
   f_const -> Draw("same");
    
   gPad -> Update();
   st_ratio= (TPaveStats*)(histo_ratio->GetListOfFunctions()->FindObject("stats"));
   st_ratio->SetX1NDC(0.); //new x start position
   st_ratio->SetX2NDC(0.); //new x end position
   st_ratio->SetY1NDC(0.); //new y start position
   st_ratio->SetY2NDC(0.); //new y end position
   st_ratio->SetTextColor(kBlack);
   st_ratio->Draw("sames");*/

   if(!log){
      c->SaveAs(std::string(outputDir_+Name+".png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+".pdf").c_str(),"pdf");	
   }else{
      c->SaveAs(std::string(outputDir_+Name+"_logY.png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+"_logY.pdf").c_str(),"pdf");
   }

   //delete histo_ratio;
}

void drawHisto(TH1F* h_good_old, TH1F* h_good_new, TH1F* h_fake_old, TH1F* h_fake_new, std::string x_label, std::string drawType, std::string Name, bool log, 
              std::string goodOldLegend, std::string goodNewLegend,std::string fakeOldLegend, std::string fakeNewLegend, string outputDir_="")
{

   h_good_old->Scale(1./h_good_old->GetEntries());
   h_good_new->Scale(1./h_good_new->GetEntries());
   h_fake_old->Scale(1./h_fake_old->GetEntries());
   h_fake_new->Scale(1./h_fake_new->GetEntries());

   std::vector<float> maxima;
   maxima.resize(4);
   maxima[0] = h_good_old->GetMaximum(); 
   maxima[1] = h_good_new->GetMaximum();
   maxima[2] = h_fake_old->GetMaximum(); 
   maxima[3] = h_fake_new->GetMaximum();
   std::sort(maxima.begin(),maxima.end());  
   
   h_good_old->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   h_good_old->SetLineColor(kRed+1);
   h_good_old->SetMarkerColor(kRed+1);
   h_good_old->SetLineWidth(2);
   h_good_old->GetXaxis()->SetTitle(x_label.c_str());

   h_good_new->SetLineColor(kMagenta-7);
   h_good_new->SetLineStyle(10);
   h_good_new->SetMarkerColor(kMagenta-7);
   h_good_new->SetLineWidth(2);

   h_fake_old->SetLineColor(kBlue+1);
   h_fake_old->SetMarkerColor(kBlue+1);
   h_fake_old->SetLineWidth(2);

   h_fake_new->SetLineColor(kCyan+1);
   h_fake_new->SetLineStyle(10);
   h_fake_new->SetMarkerColor(kCyan+1);
   h_fake_new->SetLineWidth(2);

    
   TLegend* legend = new TLegend(0.60, 0.75, 0.95, 0.87);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(40);  
   legend -> SetTextSize(0.04);

   /*TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_old_2 = new TPaveStats();
   TPaveStats* st_new_2 = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();*/
   
   TCanvas* c = new TCanvas();
   gStyle->SetOptStat(0);

   legend -> AddEntry(h_good_old,goodOldLegend.c_str(),"L");
   legend -> AddEntry(h_fake_old,fakeOldLegend.c_str(),"L");
   legend -> AddEntry(h_good_new,goodNewLegend.c_str(),"L");
   legend -> AddEntry(h_fake_new,fakeNewLegend.c_str(),"L");
   h_good_old->Draw(drawType.c_str());
   gPad -> Update();
   /*st_old= (TPaveStats*)(h_good_old->GetListOfFunctions()->FindObject("stats"));
   st_old->SetX1NDC(0.82); //new x start position
   st_old->SetX2NDC(0.99); //new x end position
   st_old->SetY1NDC(0.82); //new y start position
   st_old->SetY2NDC(0.94); //new y end position
   st_old->SetTextColor(kRed+1);*/
   //st_old->Draw("sames");
   h_good_new->Draw(std::string(drawType+",sames").c_str());
   gPad -> Update();
   /*st_new= (TPaveStats*)(h_good_new->GetListOfFunctions()->FindObject("stats"));
   st_new->SetX1NDC(0.82); //new x start position
   st_new->SetX2NDC(0.99); //new x end position
   st_new->SetY1NDC(0.68); //new y start position
   st_new->SetY2NDC(0.80); //new y end position
   st_new->SetTextColor(kBlue+1);*/
   //st_new->Draw("sames");
   h_fake_old->Draw(std::string(drawType+",sames").c_str());
   h_fake_new->Draw(std::string(drawType+",sames").c_str());
   legend -> Draw("same");


   if(!log){
      c->SaveAs(std::string(outputDir_+Name+".png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+".pdf").c_str(),"pdf");	
   }else{
      c->SaveAs(std::string(outputDir_+Name+"_logY.png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+"_logY.pdf").c_str(),"pdf");
   }

   //delete histo_ratio;
}

void drawHisto_ratio(TH1F* h_good_old, TH1F* h_good_new, TH1F* h_fake_old, TH1F* h_fake_new, std::string x_label, std::string drawType, std::string Name, bool log, 
              std::string goodOldLegend, std::string goodNewLegend,std::string fakeOldLegend, std::string fakeNewLegend,
              double ratio_min=0.8, double ratio_max=1.2, string SuperClusterRef="Mustache", string SuperClusterVal="Retuned Mustache", string outputDir_="")
{

   h_good_old->Scale(1./h_good_old->GetEntries());
   h_good_new->Scale(1./h_good_new->GetEntries());
   h_fake_old->Scale(1./h_fake_old->GetEntries());
   h_fake_new->Scale(1./h_fake_new->GetEntries());

   std::vector<float> maxima;
   maxima.resize(4);
   maxima[0] = h_good_old->GetMaximum(); 
   maxima[1] = h_good_new->GetMaximum();
   maxima[2] = h_fake_old->GetMaximum(); 
   maxima[3] = h_fake_new->GetMaximum();
   std::sort(maxima.begin(),maxima.end());  
   
   h_good_old->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   h_good_old->SetLineColor(kRed+1);
   h_good_old->SetMarkerColor(kRed+1);
   h_good_old->SetLineWidth(2);
   h_good_old->GetXaxis()->SetTitle(x_label.c_str());
   h_good_old->GetXaxis()->SetLabelSize(0.5);

   h_good_new->SetLineColor(kMagenta-7);
   h_good_new->SetLineStyle(10);
   h_good_new->SetMarkerColor(kMagenta-7);
   h_good_new->SetLineWidth(2);

   h_fake_old->SetLineColor(kBlue+1);
   h_fake_old->SetMarkerColor(kBlue+1);
   h_fake_old->SetLineWidth(2);

   h_fake_new->SetLineColor(kCyan+1);
   h_fake_new->SetLineStyle(10);
   h_fake_new->SetMarkerColor(kCyan+1);
   h_fake_new->SetLineWidth(2);

    
   TLegend* legend = new TLegend(0.60, 0.7, 0.95, 0.87);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.05);

   /*TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_old_2 = new TPaveStats();
   TPaveStats* st_new_2 = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();*/
   
   TCanvas* c = new TCanvas();
 
   TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   if(log) cUp->SetLogy();
   cDown->Draw();
     
   cUp->cd();
   gStyle->SetOptStat(0);

   legend -> AddEntry(h_good_old,goodOldLegend.c_str(),"L");
   legend -> AddEntry(h_fake_old,fakeOldLegend.c_str(),"L");
   legend -> AddEntry(h_good_new,goodNewLegend.c_str(),"L");
   legend -> AddEntry(h_fake_new,fakeNewLegend.c_str(),"L");
   h_good_old->Draw(drawType.c_str());
   gPad -> Update();
   /*st_old= (TPaveStats*)(h_good_old->GetListOfFunctions()->FindObject("stats"));
   st_old->SetX1NDC(0.82); //new x start position
   st_old->SetX2NDC(0.99); //new x end position
   st_old->SetY1NDC(0.82); //new y start position
   st_old->SetY2NDC(0.94); //new y end position
   st_old->SetTextColor(kRed+1);*/
   //st_old->Draw("sames");
   h_good_new->Draw(std::string(drawType+",sames").c_str());
   gPad -> Update();
   /*st_new= (TPaveStats*)(h_good_new->GetListOfFunctions()->FindObject("stats"));
   st_new->SetX1NDC(0.82); //new x start position
   st_new->SetX2NDC(0.99); //new x end position
   st_new->SetY1NDC(0.68); //new y start position
   st_new->SetY2NDC(0.80); //new y end position
   st_new->SetTextColor(kBlue+1);*/
   //st_new->Draw("sames");
   h_fake_old->Draw(std::string(drawType+",sames").c_str());
   h_fake_new->Draw(std::string(drawType+",sames").c_str());
   legend -> Draw("same");

    
   cDown->cd();
   TH1F* histo_ratio_good=(TH1F*)h_good_new->Clone("histo_ratio_good");
   TH1F* histo_ratio_fake=(TH1F*)h_fake_new->Clone("histo_ratio_fake");
   TPaveStats* st_ratio = new TPaveStats();
   if(histo_ratio_good->GetSumw2N()<=0) histo_ratio_good->Sumw2();
   if(histo_ratio_fake->GetSumw2N()<=0) histo_ratio_fake->Sumw2();
   histo_ratio_good->Divide(h_good_old);
   histo_ratio_fake->Divide(h_fake_old);

       
   TLegend* ratio_legend = new TLegend(0.80, 0.80, 0.92, 0.97);
   ratio_legend -> SetFillColor(kWhite);
   ratio_legend -> SetFillStyle(1000);
   ratio_legend -> SetLineWidth(0);
   ratio_legend -> SetLineColor(kWhite);
   ratio_legend -> SetTextFont(42);  
   ratio_legend -> SetTextSize(0.1);
    
   histo_ratio_good -> GetXaxis() -> SetTitle(x_label.c_str());
   //histo_ratio_good -> GetXaxis() -> SetLabelSize(0.35);
   histo_ratio_good -> GetYaxis() -> SetTitle(std::string(SuperClusterVal +"/"+SuperClusterRef).c_str());
   histo_ratio_good -> SetMaximum(ratio_max);
   histo_ratio_good -> SetMinimum(ratio_min);
   histo_ratio_good -> SetMarkerColor(kRed);
   histo_ratio_good -> SetLineColor(kRed);
   histo_ratio_good -> SetLineStyle(1);
   histo_ratio_good -> SetMarkerSize(0.5);
   histo_ratio_good ->SetTitle("");
   histo_ratio_good -> GetXaxis() -> SetLabelSize(0.1);
   histo_ratio_good -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio_good -> GetXaxis() -> SetTitleSize(0.1);
   histo_ratio_good -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio_good -> GetYaxis() -> SetTitleOffset(0.7);
   histo_ratio_good -> Draw("e");

   //histo_ratio_fake -> GetXaxis() -> SetTitle(x_label.c_str());
   //histo_ratio_fake -> GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
   //histo_ratio_fake -> SetMaximum(1.2);
   //histo_ratio_fake -> SetMinimum(0.8);
   histo_ratio_fake -> SetMarkerColor(kBlue);
   histo_ratio_fake -> SetLineColor(kBlue);
   histo_ratio_fake -> SetLineStyle(1);
   histo_ratio_fake -> SetMarkerSize(0.5);
   /*histo_ratio_fake ->SetTitle("");
   histo_ratio_fake -> GetXaxis() -> SetLabelSize(0.07);
   histo_ratio_fake -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio_fake -> GetXaxis() -> SetTitleSize(0.07);
   histo_ratio_fake -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio_fake -> GetYaxis() -> SetTitleOffset(0.7);*/
   histo_ratio_fake -> Draw("e,sames");

   TF1* f_const = new TF1("f_1", "[0]",histo_ratio_good->GetBinCenter(1)-histo_ratio_good->GetBinWidth(1)/2, histo_ratio_good->GetBinCenter(histo_ratio_good->GetNbinsX())+histo_ratio_good->GetBinWidth(histo_ratio_good->GetNbinsX())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kBlack);
   f_const -> SetLineWidth(2);
   f_const -> Draw("same");

   ratio_legend -> AddEntry(histo_ratio_good,"Photons","L");
   ratio_legend -> AddEntry(histo_ratio_fake,"Jets","L");
   ratio_legend->Draw("same");
    
   gPad -> Update();
   /*st_ratio= (TPaveStats*)(histo_ratio_good->GetListOfFunctions()->FindObject("stats"));
   st_ratio->SetX1NDC(0.); //new x start position
   st_ratio->SetX2NDC(0.); //new x end position
   st_ratio->SetY1NDC(0.); //new y start position
   st_ratio->SetY2NDC(0.); //new y end position
   st_ratio->SetTextColor(kBlack);
   st_ratio->Draw("sames");*/

   if(!log){
      c->SaveAs(std::string(outputDir_+Name+".png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+".pdf").c_str(),"pdf");	
   }else{
      c->SaveAs(std::string(outputDir_+Name+"_logY.png").c_str(),"png");
      c->SaveAs(std::string(outputDir_+Name+"_logY.pdf").c_str(),"pdf");
   }

   //delete histo_ratio;
}

void drawROC(TH1F* h_good_old, TH1F* h_good_new, TH1F* h_fake_old, TH1F* h_fake_new, std::string Name, string SuperClusterRef="Single Set", string SuperClusterVal="Local Sets",bool right = 0){

   int xBins = h_good_old -> GetNbinsX(),
       points = 10,
       pointInc = xBins / points,
       pointIdx = 0;

   TGraph *ROC_old = new TGraph(xBins), 
          *ROC_new = new TGraph(xBins),
          *ROC_points_old = new TGraph(points),
          *ROC_points_new = new TGraph(points);

   TLatex *latex[10];

   TLegend* legend = new TLegend(0.60, 0.75, 0.95, 0.87);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(40);  
   legend -> SetTextSize(0.04);

   for(int i = 1; i < xBins+1; i++){
      float effSig_new = h_good_new -> Integral(1,i) / h_good_new -> Integral(),
            effSig_old = h_good_old -> Integral(1,i) / h_good_old -> Integral(),
            effBkg_new = h_fake_new -> Integral(1,i) / h_fake_new -> Integral(),
            effBkg_old = h_fake_old -> Integral(1,i) / h_fake_old -> Integral();

      if(right){
         effSig_new = 1.0 - effSig_new;
         effSig_old = 1.0 - effSig_old;
         effBkg_new = 1.0 - effBkg_new;
         effBkg_old = 1.0 - effBkg_old;
      }

      ROC_old -> SetPoint(i-1, effSig_old, (1.0 - effBkg_old));
      ROC_new -> SetPoint(i-1, effSig_new, (1.0 - effBkg_new));

      float point_dif = 999;
      if(pointIdx == 0) point_dif = 999;
      else point_dif = effSig_old - ROC_points_old ->GetX()[pointIdx - 1];

      if(pointIdx == 0 || point_dif > 0.1){
         //cout<<"point set!\t"<<pointIdx<<endl;
         ROC_points_old -> SetPoint(pointIdx, effSig_old, (1.0 - effBkg_old));
         ROC_points_new -> SetPoint(pointIdx, effSig_new, (1.0 - effBkg_new));
         latex[pointIdx] = new TLatex(ROC_points_old ->GetX()[pointIdx], ROC_points_old ->GetY()[pointIdx],(to_string(h_good_new -> GetXaxis() -> GetBinCenter(i))).c_str());
         pointIdx++;
         
      }
   }

    
   TCanvas* c = new TCanvas();

   ROC_old -> SetLineColor(kRed);
   ROC_new -> SetLineColor(kBlue);

   ROC_points_old -> SetMarkerStyle(47);
   ROC_points_new -> SetMarkerStyle(34);

   ROC_points_old -> SetMarkerColor(kRed);
   ROC_points_new -> SetMarkerColor(kBlue);

   ROC_old -> SetLineWidth(2);
   ROC_new -> SetLineWidth(2);

   ROC_old -> SetTitle((Name).c_str());
   ROC_old -> GetYaxis() -> SetTitle("Background Rejection (1 - Background Efficiency)");
   ROC_old -> GetXaxis() -> SetTitle("Signal Efficiency");


   ROC_old -> SetMarkerStyle(47);
   ROC_old -> SetMarkerColor(kGreen);

   for(int i=0; i<pointIdx; i++){
      latex[i] -> SetTextSize(0.025);
      //ROC_points_old -> GetListOfFunctions()->Add(latex[i]);
   }

   legend -> AddEntry(ROC_old,SuperClusterRef.c_str(),"L");
   legend -> AddEntry(ROC_new,SuperClusterVal.c_str(),"L");

   ROC_old -> Draw("AC");
   ROC_new -> Draw("C SAME");
   //ROC_points_old -> Draw("P SAME");
   //ROC_points_new -> Draw("P SAME");
   legend ->Draw("SAME");


      c->SaveAs(std::string(Name+".png").c_str(),"png");
      c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");


   for(int i=0; i<pointIdx; i++){
      latex[i] -> Delete();
   }

   return;
}

void DrawHistograms(string SuperClusterRef, string SuperClusterVal){
   /*drawHisto(h_goodSC_R9, h_fakeSC_R9, string("R9"), string("hist"), string("h_SC_R9"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_SigmaIetaIeta, h_fakeSC_SigmaIetaIeta, string("SigmaIetaIeta"), string("hist"), string("h_SC_SigmaIetaIeta"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_SigmaIphiIphi, h_fakeSC_SigmaIphiIphi, string("SigmaIphiIphi"), string("hist"), string("h_SC_SigmaIphiIphi"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_fullR9, h_fakeSC_fullR9, string("R9"), string("hist"), string("h_SC_full5x5_R9"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_fullSigmaIetaIeta, h_fakeSC_fullSigmaIetaIeta, string("SigmaIetaIeta"), string("hist"), string("h_SC_full5x5_SigmaIetaIeta"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_fullSigmaIphiIphi, h_fakeSC_fullSigmaIphiIphi, string("SigmaIphiIphi"), string("hist"), string("h_SC_full5x5_SigmaIphiIphi"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_etaWidth, h_fakeSC_etaWidth, string("etaWidth"), string("hist"), string("h_SC_etaWidth"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_phiWidth, h_fakeSC_phiWidth, string("phiWidth"), string("hist"), string("h_SC_phiWidth"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_eta, h_fakeSC_eta, string("eta"), string("hist"), string("h_SC_eta"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_phi, h_fakeSC_phi, string("phi"), string("hist"), string("h_SC_phi"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_et, h_fakeSC_et, string("et"), string("hist"), string("h_SC_et"), 0, false, string(""), string("Photon SC"), string("Jet SC"));
   drawHisto(h_goodSC_energy, h_fakeSC_energy, string("energy"), string("hist"), string("h_SC_energy"), 0, false, string(""), string("Photon SC"), string("Jet SC"));*/
   
   
   for(int i =0; i<3; i++){
      
      drawHisto_ratio(h_goodSC_R9_old[i], h_goodSC_R9_new[i], h_fakeSC_R9_old[i], h_fakeSC_R9_new[i], 
                      string("R9"), string("hist"), string("h_SC_R9"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_SigmaIetaIeta_old[i],h_goodSC_SigmaIetaIeta_new[i], h_fakeSC_SigmaIetaIeta_old[i], h_fakeSC_SigmaIetaIeta_new[i], 
                      string("SigmaIetaIeta"), string("hist"), string("h_SC_SigmaIetaIeta"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_SigmaIphiIphi_old[i],h_goodSC_SigmaIphiIphi_new[i], h_fakeSC_SigmaIphiIphi_old[i], h_fakeSC_SigmaIphiIphi_new[i], 
                      string("SigmaIphiIphi"), string("hist"), string("h_SC_SigmaIphiIphi"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_fullR9_old[i], h_goodSC_fullR9_new[i], h_fakeSC_fullR9_old[i], h_fakeSC_fullR9_new[i], 
                      string("R9"), string("hist"), string("h_SC_full5x5_R9"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_fullSigmaIetaIeta_old[i],h_goodSC_fullSigmaIetaIeta_new[i], h_fakeSC_fullSigmaIetaIeta_old[i], h_fakeSC_fullSigmaIetaIeta_new[i], 
                      string("SigmaIetaIeta"), string("hist"), string("h_SC_full5x5_SigmaIetaIeta"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_fullSigmaIphiIphi_old[i],h_goodSC_fullSigmaIphiIphi_new[i], h_fakeSC_fullSigmaIphiIphi_old[i], h_fakeSC_fullSigmaIphiIphi_new[i], 
                      string("SigmaIphiIphi"), string("hist"), string("h_SC_full5x5_SigmaIphiIphi"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_etaWidth_old[i], h_goodSC_etaWidth_new[i], h_fakeSC_etaWidth_old[i], h_fakeSC_etaWidth_new[i], 
                      string("etaWidth"), string("hist"), string("h_SC_etaWidth"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_phiWidth_old[i], h_goodSC_phiWidth_new[i], h_fakeSC_phiWidth_old[i], h_fakeSC_phiWidth_new[i], 
                      string("phiWidth"), string("hist"), string("h_SC_phiWidth"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.6,1.4,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_eta_old[i], h_goodSC_eta_new[i], h_fakeSC_eta_old[i], h_fakeSC_eta_new[i], 
                      string("eta"), string("hist"), string("h_SC_eta"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.7,1.3,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_phi_old[i], h_goodSC_phi_new[i], h_fakeSC_phi_old[i], h_fakeSC_phi_new[i], 
                      string("phi"), string("hist"), string("h_SC_phi"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.7,1.3,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_et_old[i], h_goodSC_et_new[i], h_fakeSC_et_old[i], h_fakeSC_et_new[i], 
                      string("et"), string("hist"), string("h_SC_et"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.7,1.3,SuperClusterRef,SuperClusterVal);
      drawHisto_ratio(h_goodSC_energy_old[i], h_goodSC_energy_new[i], h_fakeSC_energy_old[i], h_fakeSC_energy_new[i], 
                      string("energy"), string("hist"), string("h_SC_energy"+binning[i]), 0, string(("Photon SC - " + SuperClusterRef).c_str()), string(("Photon SC - " + SuperClusterVal).c_str()),
                      string(("Jet SC - " + SuperClusterRef).c_str()), string(("Jet SC - " + SuperClusterVal).c_str()),0.7,1.3,SuperClusterRef,SuperClusterVal);
      

      drawROC(h_goodSC_R9_old[i], h_goodSC_R9_new[i], h_fakeSC_R9_old[i], h_fakeSC_R9_new[i],string("ROC_R9"+binning[i]),SuperClusterRef,SuperClusterVal, 1);
      drawROC(h_goodSC_SigmaIetaIeta_old[i],h_goodSC_SigmaIetaIeta_new[i], h_fakeSC_SigmaIetaIeta_old[i], h_fakeSC_SigmaIetaIeta_new[i],string("ROC_SigmaIetaIeta"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_SigmaIphiIphi_old[i],h_goodSC_SigmaIphiIphi_new[i], h_fakeSC_SigmaIphiIphi_old[i], h_fakeSC_SigmaIphiIphi_new[i],string("ROC_SigmaIphiIphi"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_fullR9_old[i], h_goodSC_fullR9_new[i], h_fakeSC_fullR9_old[i], h_fakeSC_fullR9_new[i],string("ROC_fullR9"+binning[i]),SuperClusterRef,SuperClusterVal, 1);
      drawROC(h_goodSC_fullSigmaIetaIeta_old[i],h_goodSC_fullSigmaIetaIeta_new[i], h_fakeSC_fullSigmaIetaIeta_old[i], h_fakeSC_fullSigmaIetaIeta_new[i],string("ROC_fullSigmaIetaIeta"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_fullSigmaIphiIphi_old[i],h_goodSC_fullSigmaIphiIphi_new[i], h_fakeSC_fullSigmaIphiIphi_old[i], h_fakeSC_fullSigmaIphiIphi_new[i],string("ROC_fullSigmaIphiIphi"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_etaWidth_old[i], h_goodSC_etaWidth_new[i], h_fakeSC_etaWidth_old[i], h_fakeSC_etaWidth_new[i],string("ROC_etaWidth"+binning[i])),SuperClusterRef,SuperClusterVal;
      drawROC(h_goodSC_phiWidth_old[i], h_goodSC_phiWidth_new[i], h_fakeSC_phiWidth_old[i], h_fakeSC_phiWidth_new[i],string("ROC_phiWidth"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_eta_old[i], h_goodSC_eta_new[i], h_fakeSC_eta_old[i], h_fakeSC_eta_new[i],string("ROC_eta"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_phi_old[i], h_goodSC_phi_new[i], h_fakeSC_phi_old[i], h_fakeSC_phi_new[i],string("ROC_phi"+binning[i]),SuperClusterRef,SuperClusterVal);
      drawROC(h_goodSC_et_old[i], h_goodSC_et_new[i], h_fakeSC_et_old[i], h_fakeSC_et_new[i],string("ROC_et"+binning[i]),SuperClusterRef,SuperClusterVal,1);
      drawROC(h_goodSC_energy_old[i], h_goodSC_energy_new[i], h_fakeSC_energy_old[i], h_fakeSC_energy_new[i],string("ROC_energy"+binning[i]),SuperClusterRef,SuperClusterVal,1);
   }


}

void GJet_Validation_LocalVSingle(){

    //Single vs Mustache
    /*string superClusterRef_label = "Mustache";
    string superClusterVal_label = "retunedMustache";
    string superClusterRef_collection = "Mustache";
    string superClusterVal_collection = "SingleSet";
    */

    //Single vs Local
    string superClusterRef_label = "SingleSet";
    string superClusterVal_label = "LocalSets";
    string superClusterRef_collection = "SingleSet";
    string superClusterVal_collection = "LocalSets";
    

    //Local vs Mustache
    /*string superClusterRef_label = "Mustache";
    string superClusterVal_label = "retunedMustache";
    string superClusterRef_collection = "Mustache";
    string superClusterVal_collection = "LocalSets";
    */
    ReadInfiles_Integrated(superClusterRef_collection, superClusterVal_collection);
    DrawHistograms(superClusterRef_collection, superClusterVal_collection);
}
