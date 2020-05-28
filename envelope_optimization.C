
#include<TApplication.h>
#include<TFile.h>
#include<TMath.h>
#include<TMinuit.h>
#include<TROOT.h>
#include<TSystem.h>
#include<TTree.h>
#include<TVector2.h>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TLegend.h>
#include<TH2F.h>
#include<TStyle.h>

#include<algorithm>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iostream>
#include<mutex>
#include<string>
#include<thread>
#include<vector>


#ifdef __MAKECINT__ 
#pragma link C++ class vector<vector<double> >+;  
#pragma link C++ class vector<vector<float> >+; 
#pragma link C++ class vector<vector<int> >+;  
#pragma link C++ class vector<vector<bool> >+; 
#pragma link C++ class vector<vector<map<int,float>> >+; 
#endif


/* ---------- Variables ---------- */
    //ofstream outfile;
    ofstream data_outfile;
    ofstream fit_outfile;

    TH2F* caloClusters_width_vs_logET[15];
    TH2F* caloClusters_shape_etWeight[15][15];
    TH2F* caloClusters_shape[15][15];
    TH2F* cluster_dPhi_vs_loget[4];

    TH2F* caloClusters_shape_xAxisWeight_Numerator[15][15];
    TH2F* caloClusters_shape_xAxisWeight_Denom[15][15];
    TH2F* caloClusters_shape_yAxisWeight_Numerator[15][15];
    TH2F* caloClusters_shape_yAxisWeight_Denom[15][15];

    float pi=3.1415927;
    float twopi= 2*pi;

    //current values
    float original_w00 = -0.00571429;
    float original_w01 = -0.002;
    float original_w10 = 0.0135714;
    float original_w11 = 0.001;

    //original dPhi values
    vector<float> yoffset_orig = {0.07151, 0.05058, -0.0913, -0.6246},
                  scale_orig   = {0.5656,  0.7131,  44.04,   13.17},
                  xoffset_orig = {0.2931,  0.01668, -5.326,  -7.037},
                  width_orig   = {0.2976,  0.4114,  1.184,   2.836};

    
    vector<float> yoffset_opt,
                  scale_opt,
                  xoffset_opt,
                  width_opt;

    //vectors and variables
    vector<float> fit_w00_upper;
    vector<float> fit_w01_upper;
    vector<float> fit_w10_upper;
    vector<float> fit_w11_upper;

    vector<float> fit_w00_lower;
    vector<float> fit_w01_lower;
    vector<float> fit_w10_lower;
    vector<float> fit_w11_lower;


    string titles_etas[15] = { "0 < |#eta| #leq 0.2",   "0.2 < |#eta| #leq 0.4", "0.4 < |#eta| #leq 0.6", 
                               "0.6 < |#eta| #leq 0.8", "0.8 < |#eta| #leq 1",   "1 < |#eta| #leq 1.2", 
                               "1.2 < |#eta| #leq 1.4", "1.4 < |#eta| #leq 1.6", "1.6 < |#eta| #leq 1.8", 
                               "1.8 < |#eta| #leq 2",   "2 < |#eta| #leq 2.2",   "2.2 < |#eta| #leq 2.4", 
                               "2.4 < |#eta| #leq 2.6", "2.6 < |#eta| #leq 2.8", "2.8 < |#eta| #leq 3" };

    string titles_loget[15] = { "-1 < log_{10}(E_{T}) < -0.8",      "-0.8 #leq log_{10}(E_{T}) < -0.6", 
                                "-0.6 #leq log_{10}(E_{T}) < -0.4", "-0.4 #leq log_{10}(E_{T}) < -0.2", 
                                "-0.2 #leq log_{10}(E_{T}) < 0",    "0 < log_{10}(E_{T}) < 0.2", 
                                "0.2 #leq log_{10}(E_{T}) < 0.4",   "0.4 #leq log_{10}(E_{T}) < 0.6", 
                                "0.6 #leq log_{10}(E_{T}) < 0.8",   "0.8 #leq log_{10}(E_{T}) < 1",
                                "1 < log_{10}(E_{T}) < 1.2",        "1.2 #leq log_{10}(E_{T}) < 1.4", 
                                "1.4 #leq log_{10}(E_{T}) < 1.6",   "1.6 #leq log_{10}(E_{T}) < 1.8", 
                                "1.8 #leq log_{10}(E_{T}) < 2" };

    string filename_etas[15] = { "00", "01", "02", "03", "04", "05", "06", "07", "08", 
                                 "09", "10", "11", "12", "13", "14" };

    string filename_loget[15] = { "-1 < log_{10}(E_{T}) < -0.8",      "-0.8 #leq log_{10}(E_{T}) < -0.6", 
                                "-0.6 #leq log_{10}(E_{T}) < -0.4", "-0.4 #leq log_{10}(E_{T}) < -0.2", 
                                "-0.2 #leq log_{10}(E_{T}) < 0",    "0 < log_{10}(E_{T}) < 0.2", 
                                "0.2 #leq log_{10}(E_{T}) < 0.4",   "0.4 #leq log_{10}(E_{T}) < 0.6", 
                                "0.6 #leq log_{10}(E_{T}) < 0.8",   "0.8 #leq log_{10}(E_{T}) < 1",
                                "1 < log_{10}(E_{T}) < 1.2",        "1.2 #leq log_{10}(E_{T}) < 1.4", 
                                "1.4 #leq log_{10}(E_{T}) < 1.6",   "1.6 #leq log_{10}(E_{T}) < 1.8", 
                                "1.8 #leq log_{10}(E_{T}) < 2" };


/* Tree Values */
    vector<int>     *genParticle_id;
    vector<float>   *genParticle_energy;
    vector<float>   *genParticle_pt;
    vector<float>   *genParticle_eta;
    vector<float>   *genParticle_phi;
    vector<vector<int> > *genParticle_pfCluster_dR_genScore_MatchedIndex;
    vector<vector<int> > *genParticle_superCluster_dR_genScore_MatchedIndex;
    vector<vector<int> > *genParticle_retunedSuperCluster_dR_genScore_MatchedIndex;
    vector<int>     *caloParticle_id;
    vector<float>   *caloParticle_genEnergy;
    vector<float>   *caloParticle_simEnergy;
    vector<float>   *caloParticle_genPt;
    vector<float>   *caloParticle_simPt;
    vector<float>   *caloParticle_genEta;
    vector<float>   *caloParticle_simEta;
    vector<float>   *caloParticle_genPhi;
    vector<float>   *caloParticle_simPhi;
    vector<int>     *caloParticle_simIeta;
    vector<int>     *caloParticle_simIphi;
    vector<int>     *caloParticle_simIz;
    vector<vector<int> > *caloParticle_pfCluster_dR_simScore_MatchedIndex;
    vector<vector<int> > *caloParticle_pfCluster_sim_fraction_old_MatchedIndex;
    vector<vector<int> > *caloParticle_pfCluster_simScore_MatchedIndex;
    vector<vector<int> > *caloParticle_superCluster_dR_simScore_MatchedIndex;
    vector<vector<int> > *caloParticle_superCluster_sim_fraction_old_MatchedIndex;
    vector<vector<int> > *caloParticle_superCluster_simScore_MatchedIndex;
    vector<vector<int> > *caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex;
    vector<vector<int> > *caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex;
    vector<vector<int> > *caloParticle_retunedSuperCluster_simScore_MatchedIndex;
    vector<vector<float> > *simHit_energy;
    vector<vector<float> > *simHit_eta;
    vector<vector<float> > *simHit_phi;
    vector<vector<int> > *simHit_ieta;
    vector<vector<int> > *simHit_iphi;
    vector<vector<int> > *simHit_iz;
    vector<float>   *recHit_noPF_energy;
    vector<float>   *recHit_noPF_eta;
    vector<float>   *recHit_noPF_phi;
    vector<int>     *recHit_noPF_ieta;
    vector<int>     *recHit_noPF_iphi;
    vector<int>     *recHit_noPF_iz;
    vector<float>   *pfRecHit_unClustered_energy;
    vector<float>   *pfRecHit_unClustered_eta;
    vector<float>   *pfRecHit_unClustered_phi;
    vector<int>     *pfRecHit_unClustered_ieta;
    vector<int>     *pfRecHit_unClustered_iphi;
    vector<int>     *pfRecHit_unClustered_iz;
    vector<float>   *pfCluster_energy;
    vector<float>   *pfCluster_eta;
    vector<float>   *pfCluster_phi;
    vector<int>     *pfCluster_ieta;
    vector<int>     *pfCluster_iphi;
    vector<int>     *pfCluster_iz;
    vector<int>     *pfCluster_nXtals;
    vector<vector<int> > *pfCluster_superClustersIndex;
    vector<vector<int> > *pfCluster_retunedSuperClustersIndex;
    vector<int>     *pfCluster_dR_genScore_MatchedIndex;
    vector<int>     *pfCluster_dR_simScore_MatchedIndex;
    vector<int>     *pfCluster_sim_fraction_old_MatchedIndex;
    vector<int>     *pfCluster_simScore_MatchedIndex;
    vector<vector<double> > *pfCluster_dR_genScore;
    vector<vector<double> > *pfCluster_dR_simScore;
    vector<vector<double> > *pfCluster_sim_fraction_old;
    vector<vector<double> > *pfCluster_simScore;
    vector<vector<float> > *pfClusterHit_energy;
    vector<vector<float> > *pfClusterHit_rechitEnergy;
    vector<vector<float> > *pfClusterHit_eta;
    vector<vector<float> > *pfClusterHit_phi;
    vector<vector<int> > *pfClusterHit_ieta;
    vector<vector<int> > *pfClusterHit_iphi;
    vector<vector<int> > *pfClusterHit_iz;
    vector<float>   *superCluster_energy;
    vector<float>   *superCluster_eta;
    vector<float>   *superCluster_phi;
    vector<float>   *superCluster_etaWidth;
    vector<float>   *superCluster_phiWidth;
    vector<float>   *superCluster_R;
    vector<int>     *superCluster_nPFClusters;
    vector<int>     *superCluster_ieta;
    vector<int>     *superCluster_iphi;
    vector<int>     *superCluster_iz;
    vector<int>     *superCluster_seedIndex;
    vector<vector<int> > *superCluster_pfClustersIndex;
    vector<vector<float> > *superCluster_psCluster_energy;
    vector<vector<float> > *superCluster_psCluster_eta;
    vector<vector<float> > *superCluster_psCluster_phi;

    vector<float>   *retunedSuperCluster_energy;
    vector<float>   *retunedSuperCluster_eta;
    vector<float>   *retunedSuperCluster_phi;
    vector<float>   *retunedSuperCluster_etaWidth;
    vector<float>   *retunedSuperCluster_phiWidth;
    vector<float>   *retunedSuperCluster_R;
    vector<int>     *retunedSuperCluster_nPFClusters;
    vector<int>     *retunedSuperCluster_ieta;
    vector<int>     *retunedSuperCluster_iphi;
    vector<int>     *retunedSuperCluster_iz;
    vector<int>     *retunedSuperCluster_seedIndex;
    vector<vector<int> > *retunedSuperCluster_pfClustersIndex;
    vector<vector<float> > *retunedSuperCluster_psCluster_energy;
    vector<vector<float> > *retunedSuperCluster_psCluster_eta;
    vector<vector<float> > *retunedSuperCluster_psCluster_phi;

    vector<int>     *superCluster_dR_genScore_MatchedIndex;
    vector<int>     *superCluster_dR_simScore_MatchedIndex;
    vector<int>     *superCluster_sim_fraction_old_MatchedIndex;
    vector<int>     *superCluster_simScore_MatchedIndex;
    vector<int>     *retunedSuperCluster_dR_genScore_MatchedIndex;
    vector<int>     *retunedSuperCluster_dR_simScore_MatchedIndex;
    vector<int>     *retunedSuperCluster_sim_fraction_old_MatchedIndex;
    vector<int>     *retunedSuperCluster_simScore_MatchedIndex;
    vector<vector<double> > *superCluster_dR_genScore;
    vector<vector<double> > *superCluster_dR_simScore;
    vector<vector<double> > *superCluster_sim_fraction_old;
    vector<vector<double> > *superCluster_simScore;
    vector<vector<double> > *retunedSuperCluster_dR_genScore;
    vector<vector<double> > *retunedSuperCluster_dR_simScore;
    vector<vector<double> > *retunedSuperCluster_sim_fraction_old;
    vector<vector<double> > *retunedSuperCluster_simScore;
    
    vector<float>   *pfCluster_etaWidth;
    vector<float>   *pfCluster_phiWidth;
    vector<float>   *pfCluster_e5x5;
    vector<float>   *pfCluster_e2x2Ratio;
    vector<float>   *pfCluster_e3x3Ratio;
    vector<float>   *pfCluster_eMaxRatio;
    vector<float>   *pfCluster_e2ndRatio;
    vector<float>   *pfCluster_eTopRatio;
    vector<float>   *pfCluster_eRightRatio;
    vector<float>   *pfCluster_eBottomRatio;
    vector<float>   *pfCluster_eLeftRatio;
    vector<float>   *pfCluster_e2x5MaxRatio;
    vector<float>   *pfCluster_e2x5TopRatio;
    vector<float>   *pfCluster_e2x5RightRatio;
    vector<float>   *pfCluster_e2x5BottomRatio;
    vector<float>   *pfCluster_e2x5LeftRatio;
    vector<float>   *pfCluster_swissCross;
    vector<float>   *pfCluster_r9;
    vector<float>   *pfCluster_sigmaIetaIeta;
    vector<float>   *pfCluster_sigmaIetaIphi;
    vector<float>   *pfCluster_sigmaIphiIphi;
    vector<float>   *pfCluster_full5x5_e5x5;
    vector<float>   *pfCluster_full5x5_e2x2Ratio;
    vector<float>   *pfCluster_full5x5_e3x3Ratio;
    vector<float>   *pfCluster_full5x5_eMaxRatio;
    vector<float>   *pfCluster_full5x5_e2ndRatio;
    vector<float>   *pfCluster_full5x5_eTopRatio;
    vector<float>   *pfCluster_full5x5_eRightRatio;
    vector<float>   *pfCluster_full5x5_eBottomRatio;
    vector<float>   *pfCluster_full5x5_eLeftRatio;
    vector<float>   *pfCluster_full5x5_e2x5MaxRatio;
    vector<float>   *pfCluster_full5x5_e2x5TopRatio;
    vector<float>   *pfCluster_full5x5_e2x5RightRatio;
    vector<float>   *pfCluster_full5x5_e2x5BottomRatio;
    vector<float>   *pfCluster_full5x5_e2x5LeftRatio;
    vector<float>   *pfCluster_full5x5_swissCross;
    vector<float>   *pfCluster_full5x5_r9;
    vector<float>   *pfCluster_full5x5_sigmaIetaIeta;
    vector<float>   *pfCluster_full5x5_sigmaIetaIphi;
    vector<float>   *pfCluster_full5x5_sigmaIphiIphi;

    vector<float>   *superCluster_e5x5;
    vector<float>   *superCluster_e2x2Ratio;
    vector<float>   *superCluster_e3x3Ratio;
    vector<float>   *superCluster_eMaxRatio;
    vector<float>   *superCluster_e2ndRatio;
    vector<float>   *superCluster_eTopRatio;
    vector<float>   *superCluster_eRightRatio;
    vector<float>   *superCluster_eBottomRatio;
    vector<float>   *superCluster_eLeftRatio;
    vector<float>   *superCluster_e2x5MaxRatio;
    vector<float>   *superCluster_e2x5TopRatio;
    vector<float>   *superCluster_e2x5RightRatio;
    vector<float>   *superCluster_e2x5BottomRatio;
    vector<float>   *superCluster_e2x5LeftRatio;
    vector<float>   *superCluster_swissCross;
    vector<float>   *superCluster_r9;
    vector<float>   *superCluster_sigmaIetaIeta;
    vector<float>   *superCluster_sigmaIetaIphi;
    vector<float>   *superCluster_sigmaIphiIphi;
    vector<float>   *superCluster_full5x5_e5x5;
    vector<float>   *superCluster_full5x5_e2x2Ratio;
    vector<float>   *superCluster_full5x5_e3x3Ratio;
    vector<float>   *superCluster_full5x5_eMaxRatio;
    vector<float>   *superCluster_full5x5_e2ndRatio;
    vector<float>   *superCluster_full5x5_eTopRatio;
    vector<float>   *superCluster_full5x5_eRightRatio;
    vector<float>   *superCluster_full5x5_eBottomRatio;
    vector<float>   *superCluster_full5x5_eLeftRatio;
    vector<float>   *superCluster_full5x5_e2x5MaxRatio;
    vector<float>   *superCluster_full5x5_e2x5TopRatio;
    vector<float>   *superCluster_full5x5_e2x5RightRatio;
    vector<float>   *superCluster_full5x5_e2x5BottomRatio;
    vector<float>   *superCluster_full5x5_e2x5LeftRatio;
    vector<float>   *superCluster_full5x5_swissCross;
    vector<float>   *superCluster_full5x5_r9;
    vector<float>   *superCluster_full5x5_sigmaIetaIeta;
    vector<float>   *superCluster_full5x5_sigmaIetaIphi;
    vector<float>   *superCluster_full5x5_sigmaIphiIphi;
    vector<float>   *superCluster_HoEraw;
    vector<float>   *superCluster_HoErawBC;

    vector<float>   *retunedSuperCluster_e5x5;
    vector<float>   *retunedSuperCluster_e2x2Ratio;
    vector<float>   *retunedSuperCluster_e3x3Ratio;
    vector<float>   *retunedSuperCluster_eMaxRatio;
    vector<float>   *retunedSuperCluster_e2ndRatio;
    vector<float>   *retunedSuperCluster_eTopRatio;
    vector<float>   *retunedSuperCluster_eRightRatio;
    vector<float>   *retunedSuperCluster_eBottomRatio;
    vector<float>   *retunedSuperCluster_eLeftRatio;
    vector<float>   *retunedSuperCluster_e2x5MaxRatio;
    vector<float>   *retunedSuperCluster_e2x5TopRatio;
    vector<float>   *retunedSuperCluster_e2x5RightRatio;
    vector<float>   *retunedSuperCluster_e2x5BottomRatio;
    vector<float>   *retunedSuperCluster_e2x5LeftRatio;
    vector<float>   *retunedSuperCluster_swissCross;
    vector<float>   *retunedSuperCluster_r9;
    vector<float>   *retunedSuperCluster_sigmaIetaIeta;
    vector<float>   *retunedSuperCluster_sigmaIetaIphi;
    vector<float>   *retunedSuperCluster_sigmaIphiIphi;
    vector<float>   *retunedSuperCluster_full5x5_e5x5;
    vector<float>   *retunedSuperCluster_full5x5_e2x2Ratio;
    vector<float>   *retunedSuperCluster_full5x5_e3x3Ratio;
    vector<float>   *retunedSuperCluster_full5x5_eMaxRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2ndRatio;
    vector<float>   *retunedSuperCluster_full5x5_eTopRatio;
    vector<float>   *retunedSuperCluster_full5x5_eRightRatio;
    vector<float>   *retunedSuperCluster_full5x5_eBottomRatio;
    vector<float>   *retunedSuperCluster_full5x5_eLeftRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2x5MaxRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2x5TopRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2x5RightRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2x5BottomRatio;
    vector<float>   *retunedSuperCluster_full5x5_e2x5LeftRatio;
    vector<float>   *retunedSuperCluster_full5x5_swissCross;
    vector<float>   *retunedSuperCluster_full5x5_r9;
    vector<float>   *retunedSuperCluster_full5x5_sigmaIetaIeta;
    vector<float>   *retunedSuperCluster_full5x5_sigmaIetaIphi;
    vector<float>   *retunedSuperCluster_full5x5_sigmaIphiIphi;
    vector<float>   *retunedSuperCluster_HoEraw;
    vector<float>   *retunedSuperCluster_HoErawBC;

    TTree* EventTree;
    Int_t EvMax;
/* */

string cur_time(){
//returns the current time, for logging
    std::time_t tt = std::time(NULL);
    std::string s = std::ctime(&tt);
    return s.substr(0, s.size()-1);
}

string trimString(string str){
//trim whitespace at end of string
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

bool comparePFpt(const std::pair<int, double>&i, const std::pair<int, double>&j)
{
   return i.second > j.second;
}

void EventLoop(){
//Loops over all events and fills histograms
    int eventLoopMax = EvMax;

    for(int iev=0; iev<eventLoopMax;++iev){
        EventTree->GetEvent(iev);
        if(iev % 50 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;

        //get all pfClusters matched to caloParticles
        std::vector< std::vector< std::pair<int, double> > > pfMatched_idx_pt(caloParticle_id->size());
        for(unsigned int iPF=0; iPF<pfCluster_energy->size(); iPF++){
            int caloMatch = pfCluster_simScore_MatchedIndex->at(iPF);
            if(caloMatch >= 0)
                pfMatched_idx_pt.at(caloMatch).push_back(std::make_pair(iPF, pfCluster_energy->at(iPF)/TMath::CosH(pfCluster_eta->at(iPF))));
        }

        //loop over caloParticles
        for(unsigned int idxCalo=0; idxCalo<caloParticle_id->size(); idxCalo++){

            //order matched pfClusters in pt
            sort(pfMatched_idx_pt.at(idxCalo).begin(),pfMatched_idx_pt.at(idxCalo).end(),comparePFpt);

            if(pfMatched_idx_pt.at(idxCalo).size() == 0) continue;

            int seedCluster = pfMatched_idx_pt.at(idxCalo).at(0).first;
            float seedEta = pfCluster_eta->at(seedCluster);
            int dPhi_bin = -1, seedEta_bin = -1;

            if     (fabs(seedEta) >= 0.0 && fabs(seedEta) < 0.2){ seedEta_bin = 0; }
            else if(fabs(seedEta) >= 0.2 && fabs(seedEta) < 0.4){ seedEta_bin = 1; }
            else if(fabs(seedEta) >= 0.4 && fabs(seedEta) < 0.6){ seedEta_bin = 2; } 
            else if(fabs(seedEta) >= 0.6 && fabs(seedEta) < 0.8){ seedEta_bin = 3; }
            else if(fabs(seedEta) >= 0.8 && fabs(seedEta) < 1.0){ seedEta_bin = 4; }
            else if(fabs(seedEta) >= 1.0 && fabs(seedEta) < 1.2){ seedEta_bin = 5; } 
            else if(fabs(seedEta) >= 1.2 && fabs(seedEta) < 1.4){ seedEta_bin = 6; }
            else if(fabs(seedEta) >= 1.4 && fabs(seedEta) < 1.6){ seedEta_bin = 7; }
            else if(fabs(seedEta) >= 1.6 && fabs(seedEta) < 1.8){ seedEta_bin = 8; }
            else if(fabs(seedEta) >= 1.8 && fabs(seedEta) < 2.0){ seedEta_bin = 9; }
            else if(fabs(seedEta) >= 2.0 && fabs(seedEta) < 2.2){ seedEta_bin = 10; } 
            else if(fabs(seedEta) >= 2.2 && fabs(seedEta) < 2.4){ seedEta_bin = 11; } 
            else if(fabs(seedEta) >= 2.4 && fabs(seedEta) < 2.6){ seedEta_bin = 12; } 
            else if(fabs(seedEta) >= 2.6 && fabs(seedEta) < 2.8){ seedEta_bin = 13; } 
            else if(fabs(seedEta) >= 2.8 && fabs(seedEta) < 3.0){ seedEta_bin = 14; } 

            if     (fabs(seedEta) >= 0.0   && fabs(seedEta) < 1.479){ dPhi_bin = 0; }
            else if(fabs(seedEta) >= 1.479 && fabs(seedEta) < 1.75) { dPhi_bin = 1; }
            else if(fabs(seedEta) >= 1.75  && fabs(seedEta) < 2.0)  { dPhi_bin = 2; } 
            else if(fabs(seedEta) >= 2.0   && fabs(seedEta) < 3.0)  { dPhi_bin = 3; }

            //loop over all matched pfClusters
            for(auto pfMatch_idx_pt_el : pfMatched_idx_pt.at(idxCalo)){
                int idxPF = pfMatch_idx_pt_el.first;
                double etPF = pfMatch_idx_pt_el.second;

                //skip the seeds
                if(idxPF == seedCluster) continue;

                float dEta = (1-2*(pfCluster_eta->at(seedCluster) < 0)) * (pfCluster_eta->at(idxPF) - pfCluster_eta->at(seedCluster));
                float dPhi = pfCluster_phi->at(idxPF) - pfCluster_phi->at(seedCluster);
                if(dPhi > pi) dPhi -=twopi;
                if(dPhi < -pi) dPhi += twopi;
                float logET = log10(etPF);

                int logET_bin = -1;

                if     (logET >= -1.0 && logET < -0.8) logET_bin = 0;
                else if(logET >= -0.8 && logET < -0.6) logET_bin = 1;
                else if(logET >= -0.6 && logET < -0.4) logET_bin = 2;
                else if(logET >= -0.4 && logET < -0.2) logET_bin = 3;
                else if(logET >= -0.2 && logET < 0.0)  logET_bin = 4;
                else if(logET >= 0.0  && logET < 0.2)  logET_bin = 5;
                else if(logET >= 0.2  && logET < 0.4)  logET_bin = 6;
                else if(logET >= 0.4  && logET < 0.6)  logET_bin = 7;
                else if(logET >= 0.6  && logET < 0.8)  logET_bin = 8;
                else if(logET >= 0.8  && logET < 1.0)  logET_bin = 9;
                else if(logET >= 1.0  && logET < 1.2)  logET_bin = 10;
                else if(logET >= 1.2  && logET < 1.4)  logET_bin = 11;
                else if(logET >= 1.4  && logET < 1.6)  logET_bin = 12;
                else if(logET >= 1.6  && logET < 1.8)  logET_bin = 13;
                else if(logET >= 1.8  && logET < 2.0)  logET_bin = 14;

                if(logET_bin < 0 || seedEta_bin < 0) continue;

                //fill histograms
                caloClusters_shape_etWeight[seedEta_bin][logET_bin] -> Fill(dPhi, dEta, etPF);
                caloClusters_shape[seedEta_bin][logET_bin] -> Fill(dPhi, dEta);
                cluster_dPhi_vs_loget[dPhi_bin] -> Fill(logET, fabs(dPhi), etPF);

                caloClusters_shape_xAxisWeight_Numerator[seedEta_bin][logET_bin] -> Fill(dPhi, dEta, (etPF * pfCluster_phiWidth->at(idxPF)));
                caloClusters_shape_xAxisWeight_Denom[seedEta_bin][logET_bin] -> Fill(dPhi, dEta, etPF);
                caloClusters_shape_yAxisWeight_Numerator[seedEta_bin][logET_bin] -> Fill(dPhi, dEta, (etPF * pfCluster_etaWidth->at(idxPF)));
                caloClusters_shape_yAxisWeight_Denom[seedEta_bin][logET_bin] -> Fill(dPhi, dEta, etPF);
            }

        }// loop over calo particels


    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;
}

void InitHistograms(){

    fit_outfile.open("envelope/fit_status.txt", std::ofstream::out | std::ofstream::trunc);

    double min_dEta = -0.15,
           max_dEta = 0.2,
           min_dPhi = -0.6,
           max_dPhi = 0.6;
    
    int bins_dEta = 105,
        bins_dPhi = 120;

    for(int i=0;i<15;i++){ //seedEta bins
        caloClusters_width_vs_logET[i] = new TH2F(("caloClusters_width_vs_logET_etaBin_"+to_string(i)).c_str(),("caloClusters_width_vs_logET_etaBin_"+to_string(i)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
        for(int j=0; j<15;j++){ //loget bins
            caloClusters_shape[i][j] = new TH2F(("caloClusters_shape_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_etWeight[i][j] = new TH2F(("caloClusters_shape_etWeight_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_etWeight_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_xAxisWeight_Numerator[i][j] = new TH2F(("caloClusters_shape_xAxisWeight_Numerator_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_xAxisWeight_Numerator_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_yAxisWeight_Numerator[i][j] = new TH2F(("caloClusters_shape_yAxisWeight_Numerator_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_yAxisWeight_Numerator_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_xAxisWeight_Denom[i][j] = new TH2F(("caloClusters_shape_xAxisWeight_Denom_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_xAxisWeight_Denom_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_yAxisWeight_Denom[i][j] = new TH2F(("caloClusters_shape_yAxisWeight_Denom_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_yAxisWeight_Denom_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
        }
    }
    
    for(int ii=0; ii<4; ii++){
        cluster_dPhi_vs_loget[ii] = new TH2F(("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),200, -1.0, 3.0, 200, 0.0, 1.0);
    }

}

void InitTree(TString infileName,float weight=1){

    cout<<cur_time()<<"\tProcessing Tree...\n";
    TFile* infile = new TFile(infileName);
    infile->cd("recosimdumper");
    EventTree=(TTree*)gDirectory->Get("caloTree");


    //EventTree->SetBranchAddress("eventId", &eventId);
    //EventTree->SetBranchAddress("lumiId", &lumiId);
    //EventTree->SetBranchAddress("runId", &runId);
    //EventTree->SetBranchAddress("nVtx", &nVtx);
    //EventTree->SetBranchAddress("rho", &rho);
    EventTree->SetBranchAddress("genParticle_id", &genParticle_id);
    EventTree->SetBranchAddress("genParticle_energy", &genParticle_energy);
    EventTree->SetBranchAddress("genParticle_pt", &genParticle_pt);
    EventTree->SetBranchAddress("genParticle_eta", &genParticle_eta);
    EventTree->SetBranchAddress("genParticle_phi", &genParticle_phi);
    EventTree->SetBranchAddress("genParticle_pfCluster_dR_genScore_MatchedIndex", &genParticle_pfCluster_dR_genScore_MatchedIndex);
    EventTree->SetBranchAddress("genParticle_superCluster_dR_genScore_MatchedIndex", &genParticle_superCluster_dR_genScore_MatchedIndex);
    //EventTree->SetBranchAddress("genParticle_retunedSuperCluster_dR_genScore_MatchedIndex", &genParticle_retunedSuperCluster_dR_genScore_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_id", &caloParticle_id);//_caloParticle_id);
    EventTree->SetBranchAddress("caloParticle_genEnergy", &caloParticle_genEnergy);
    EventTree->SetBranchAddress("caloParticle_simEnergy", &caloParticle_simEnergy);
    EventTree->SetBranchAddress("caloParticle_genPt", &caloParticle_genPt);
    EventTree->SetBranchAddress("caloParticle_simPt", &caloParticle_simPt);
    EventTree->SetBranchAddress("caloParticle_genEta", &caloParticle_genEta);
    EventTree->SetBranchAddress("caloParticle_simEta", &caloParticle_simEta);
    EventTree->SetBranchAddress("caloParticle_genPhi", &caloParticle_genPhi);
    EventTree->SetBranchAddress("caloParticle_simPhi", &caloParticle_simPhi);
    EventTree->SetBranchAddress("caloParticle_simIeta", &caloParticle_simIeta);
    EventTree->SetBranchAddress("caloParticle_simIphi", &caloParticle_simIphi);
    EventTree->SetBranchAddress("caloParticle_simIz", &caloParticle_simIz);
    EventTree->SetBranchAddress("caloParticle_pfCluster_dR_simScore_MatchedIndex", &caloParticle_pfCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_pfCluster_sim_fraction_old_MatchedIndex", &caloParticle_pfCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_pfCluster_simScore_MatchedIndex", &caloParticle_pfCluster_simScore_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_superCluster_dR_simScore_MatchedIndex", &caloParticle_superCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_superCluster_sim_fraction_old_MatchedIndex", &caloParticle_superCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_superCluster_simScore_MatchedIndex", &caloParticle_superCluster_simScore_MatchedIndex);
    /*EventTree->SetBranchAddress("caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex", &caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex", &caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_retunedSuperCluster_simScore_MatchedIndex", &caloParticle_retunedSuperCluster_simScore_MatchedIndex);*/
    EventTree->SetBranchAddress("simHit_energy", &simHit_energy);
    EventTree->SetBranchAddress("simHit_eta", &simHit_eta);
    EventTree->SetBranchAddress("simHit_phi", &simHit_phi);
    EventTree->SetBranchAddress("simHit_ieta", &simHit_ieta);
    EventTree->SetBranchAddress("simHit_iphi", &simHit_iphi);
    EventTree->SetBranchAddress("simHit_iz", &simHit_iz);
    EventTree->SetBranchAddress("recHit_noPF_energy", &recHit_noPF_energy);
    EventTree->SetBranchAddress("recHit_noPF_eta", &recHit_noPF_eta);
    EventTree->SetBranchAddress("recHit_noPF_phi", &recHit_noPF_phi);
    EventTree->SetBranchAddress("recHit_noPF_ieta", &recHit_noPF_ieta);
    EventTree->SetBranchAddress("recHit_noPF_iphi", &recHit_noPF_iphi);//_recHit_noPF_iphi);
    EventTree->SetBranchAddress("recHit_noPF_iz", &recHit_noPF_iz);//_recHit_noPF_iz);
    EventTree->SetBranchAddress("pfRecHit_unClustered_energy", &pfRecHit_unClustered_energy);//_pfRecHit_unClustered_energy);
    EventTree->SetBranchAddress("pfRecHit_unClustered_eta", &pfRecHit_unClustered_eta);//_pfRecHit_unClustered_eta);
    EventTree->SetBranchAddress("pfRecHit_unClustered_phi", &pfRecHit_unClustered_phi);//_pfRecHit_unClustered_phi);
    EventTree->SetBranchAddress("pfRecHit_unClustered_ieta", &pfRecHit_unClustered_ieta);//_pfRecHit_unClustered_ieta);
    EventTree->SetBranchAddress("pfRecHit_unClustered_iphi", &pfRecHit_unClustered_iphi);//_pfRecHit_unClustered_iphi);
    EventTree->SetBranchAddress("pfRecHit_unClustered_iz", &pfRecHit_unClustered_iz);//_pfRecHit_unClustered_iz);
    EventTree->SetBranchAddress("pfCluster_energy", &pfCluster_energy);//_pfCluster_energy);
    EventTree->SetBranchAddress("pfCluster_eta", &pfCluster_eta);//_pfCluster_eta);
    EventTree->SetBranchAddress("pfCluster_phi", &pfCluster_phi);//_pfCluster_phi);
    EventTree->SetBranchAddress("pfCluster_ieta", &pfCluster_ieta);//_pfCluster_ieta);
    EventTree->SetBranchAddress("pfCluster_iphi", &pfCluster_iphi);//_pfCluster_iphi);
    EventTree->SetBranchAddress("pfCluster_iz", &pfCluster_iz);//_pfCluster_iz);
    EventTree->SetBranchAddress("pfCluster_nXtals", &pfCluster_nXtals);//_pfCluster_nXtals);
    EventTree->SetBranchAddress("pfCluster_superClustersIndex", &pfCluster_superClustersIndex);//_pfCluster_superClustersIndex);
    //EventTree->SetBranchAddress("pfCluster_retunedSuperClustersIndex", &pfCluster_retunedSuperClustersIndex);//_pfCluster_retunedSuperClustersIndex);
    EventTree->SetBranchAddress("pfCluster_dR_genScore_MatchedIndex", &pfCluster_dR_genScore_MatchedIndex);//_pfCluster_dR_genScore_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_dR_simScore_MatchedIndex", &pfCluster_dR_simScore_MatchedIndex);//_pfCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_sim_fraction_old_MatchedIndex", &pfCluster_sim_fraction_old_MatchedIndex);//_pfCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_simScore_MatchedIndex", &pfCluster_simScore_MatchedIndex);//_pfCluster_simScore_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_dR_genScore", &pfCluster_dR_genScore);//_pfCluster_dR_genScore);
    EventTree->SetBranchAddress("pfCluster_dR_simScore", &pfCluster_dR_simScore);//_pfCluster_dR_simScore);
    EventTree->SetBranchAddress("pfCluster_sim_fraction_old", &pfCluster_sim_fraction_old);//_pfCluster_sim_fraction_old);
    EventTree->SetBranchAddress("pfCluster_simScore", &pfCluster_simScore);//_pfCluster_simScore);
    EventTree->SetBranchAddress("pfClusterHit_energy", &pfClusterHit_energy);//_pfClusterHit_energy);
    EventTree->SetBranchAddress("pfClusterHit_rechitEnergy", &pfClusterHit_rechitEnergy);//_pfClusterHit_rechitEnergy);
    EventTree->SetBranchAddress("pfClusterHit_eta", &pfClusterHit_eta);//_pfClusterHit_eta);
    EventTree->SetBranchAddress("pfClusterHit_phi", &pfClusterHit_phi);//_pfClusterHit_phi);
    EventTree->SetBranchAddress("pfClusterHit_ieta", &pfClusterHit_ieta);//_pfClusterHit_ieta);
    EventTree->SetBranchAddress("pfClusterHit_iphi", &pfClusterHit_iphi);//_pfClusterHit_iphi);
    EventTree->SetBranchAddress("pfClusterHit_iz", &pfClusterHit_iz);//_pfClusterHit_iz);
    EventTree->SetBranchAddress("superCluster_energy", &superCluster_energy);//_superCluster_energy);
    EventTree->SetBranchAddress("superCluster_eta", &superCluster_eta);//_superCluster_eta);
    EventTree->SetBranchAddress("superCluster_phi", &superCluster_phi);//_superCluster_phi);
    EventTree->SetBranchAddress("superCluster_etaWidth", &superCluster_etaWidth);//_superCluster_etaWidth);
    EventTree->SetBranchAddress("superCluster_phiWidth", &superCluster_phiWidth);//_superCluster_phiWidth);
    EventTree->SetBranchAddress("superCluster_R", &superCluster_R);//_superCluster_R);
    EventTree->SetBranchAddress("superCluster_nPFClusters", &superCluster_nPFClusters);//_superCluster_nPFClusters);
    EventTree->SetBranchAddress("superCluster_ieta", &superCluster_ieta);//_superCluster_ieta);
    EventTree->SetBranchAddress("superCluster_iphi", &superCluster_iphi);//_superCluster_iphi);
    EventTree->SetBranchAddress("superCluster_iz", &superCluster_iz);//_superCluster_iz);
    EventTree->SetBranchAddress("superCluster_seedIndex", &superCluster_seedIndex);//_superCluster_seedIndex);
    EventTree->SetBranchAddress("superCluster_pfClustersIndex", &superCluster_pfClustersIndex);//_superCluster_pfClustersIndex);
    EventTree->SetBranchAddress("superCluster_psCluster_energy", &superCluster_psCluster_energy);//_superCluster_psCluster_energy);
    EventTree->SetBranchAddress("superCluster_psCluster_eta", &superCluster_psCluster_eta);//_superCluster_psCluster_eta);
    EventTree->SetBranchAddress("superCluster_psCluster_phi", &superCluster_psCluster_phi);//_superCluster_psCluster_phi);
    /*EventTree->SetBranchAddress("retunedSuperCluster_energy", &retunedSuperCluster_energy);//_retunedSuperCluster_energy);
    EventTree->SetBranchAddress("retunedSuperCluster_eta", &retunedSuperCluster_eta);//_retunedSuperCluster_eta);
    EventTree->SetBranchAddress("retunedSuperCluster_phi", &retunedSuperCluster_phi);//_retunedSuperCluster_phi);
    EventTree->SetBranchAddress("retunedSuperCluster_etaWidth", &retunedSuperCluster_etaWidth);//_retunedSuperCluster_etaWidth);
    EventTree->SetBranchAddress("retunedSuperCluster_phiWidth", &retunedSuperCluster_phiWidth);//_retunedSuperCluster_phiWidth);
    EventTree->SetBranchAddress("retunedSuperCluster_R", &retunedSuperCluster_R);//_retunedSuperCluster_R);
    EventTree->SetBranchAddress("retunedSuperCluster_nPFClusters", &retunedSuperCluster_nPFClusters);//_retunedSuperCluster_nPFClusters);
    EventTree->SetBranchAddress("retunedSuperCluster_ieta", &retunedSuperCluster_ieta);//_retunedSuperCluster_ieta);
    EventTree->SetBranchAddress("retunedSuperCluster_iphi", &retunedSuperCluster_iphi);//_retunedSuperCluster_iphi);
    EventTree->SetBranchAddress("retunedSuperCluster_iz", &retunedSuperCluster_iz);//_retunedSuperCluster_iz);
    EventTree->SetBranchAddress("retunedSuperCluster_seedIndex", &retunedSuperCluster_seedIndex);//_retunedSuperCluster_seedIndex);
    EventTree->SetBranchAddress("retunedSuperCluster_pfClustersIndex", &retunedSuperCluster_pfClustersIndex);//_retunedSuperCluster_pfClustersIndex);
    EventTree->SetBranchAddress("retunedSuperCluster_psCluster_energy", &retunedSuperCluster_psCluster_energy);//_retunedSuperCluster_psCluster_energy);
    EventTree->SetBranchAddress("retunedSuperCluster_psCluster_eta", &retunedSuperCluster_psCluster_eta);//_retunedSuperCluster_psCluster_eta);
    EventTree->SetBranchAddress("retunedSuperCluster_psCluster_phi", &retunedSuperCluster_psCluster_phi);//_retunedSuperCluster_psCluster_phi);*/
    EventTree->SetBranchAddress("superCluster_dR_genScore_MatchedIndex", &superCluster_dR_genScore_MatchedIndex);//_superCluster_dR_genScore_MatchedIndex);
    EventTree->SetBranchAddress("superCluster_dR_simScore_MatchedIndex", &superCluster_dR_simScore_MatchedIndex);//_superCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("superCluster_sim_fraction_old_MatchedIndex", &superCluster_sim_fraction_old_MatchedIndex);//_superCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("superCluster_simScore_MatchedIndex", &superCluster_simScore_MatchedIndex);//_superCluster_simScore_MatchedIndex);
    /*EventTree->SetBranchAddress("retunedSuperCluster_dR_genScore_MatchedIndex", &retunedSuperCluster_dR_genScore_MatchedIndex);//_retunedSuperCluster_dR_genScore_MatchedIndex);
    EventTree->SetBranchAddress("retunedSuperCluster_dR_simScore_MatchedIndex", &retunedSuperCluster_dR_simScore_MatchedIndex);//_retunedSuperCluster_dR_simScore_MatchedIndex);
    EventTree->SetBranchAddress("retunedSuperCluster_sim_fraction_old_MatchedIndex", &retunedSuperCluster_sim_fraction_old_MatchedIndex);//_retunedSuperCluster_sim_fraction_old_MatchedIndex);
    EventTree->SetBranchAddress("retunedSuperCluster_simScore_MatchedIndex", &retunedSuperCluster_simScore_MatchedIndex);//_retunedSuperCluster_simScore_MatchedIndex);
 */EventTree->SetBranchAddress("superCluster_dR_genScore", &superCluster_dR_genScore);//_superCluster_dR_genScore);
    EventTree->SetBranchAddress("superCluster_dR_simScore", &superCluster_dR_simScore);//_superCluster_dR_simScore);
    EventTree->SetBranchAddress("superCluster_sim_fraction_old", &superCluster_sim_fraction_old);//_superCluster_sim_fraction_old);
    EventTree->SetBranchAddress("superCluster_simScore", &superCluster_simScore);//_superCluster_simScore);
    /*EventTree->SetBranchAddress("retunedSuperCluster_dR_genScore", &retunedSuperCluster_dR_genScore);//_retunedSuperCluster_dR_genScore);
    EventTree->SetBranchAddress("retunedSuperCluster_dR_simScore", &retunedSuperCluster_dR_simScore);//_retunedSuperCluster_dR_simScore);
    EventTree->SetBranchAddress("retunedSuperCluster_sim_fraction_old", &retunedSuperCluster_sim_fraction_old);//_retunedSuperCluster_sim_fraction_old);
    EventTree->SetBranchAddress("retunedSuperCluster_simScore", &retunedSuperCluster_simScore);//_retunedSuperCluster_simScore);*/

    EventTree->SetBranchAddress("pfCluster_etaWidth", &pfCluster_etaWidth);//_pfCluster_e5x5);
    EventTree->SetBranchAddress("pfCluster_phiWidth", &pfCluster_phiWidth);//_pfCluster_e2x2Ratio);
    EventTree->SetBranchAddress("pfCluster_e5x5", &pfCluster_e5x5);//_pfCluster_e5x5);
    EventTree->SetBranchAddress("pfCluster_e2x2Ratio", &pfCluster_e2x2Ratio);//_pfCluster_e2x2Ratio);
    EventTree->SetBranchAddress("pfCluster_e3x3Ratio", &pfCluster_e3x3Ratio);//_pfCluster_e3x3Ratio);
    EventTree->SetBranchAddress("pfCluster_eMaxRatio", &pfCluster_eMaxRatio);//_pfCluster_eMaxRatio);
    EventTree->SetBranchAddress("pfCluster_e2ndRatio", &pfCluster_e2ndRatio);//_pfCluster_e2ndRatio);
    EventTree->SetBranchAddress("pfCluster_eTopRatio", &pfCluster_eTopRatio);//_pfCluster_eTopRatio);
    EventTree->SetBranchAddress("pfCluster_eRightRatio", &pfCluster_eRightRatio);//_pfCluster_eRightRatio);
    EventTree->SetBranchAddress("pfCluster_eBottomRatio", &pfCluster_eBottomRatio);//_pfCluster_eBottomRatio);
    EventTree->SetBranchAddress("pfCluster_eLeftRatio", &pfCluster_eLeftRatio);//_pfCluster_eLeftRatio);
    EventTree->SetBranchAddress("pfCluster_e2x5MaxRatio", &pfCluster_e2x5MaxRatio);//_pfCluster_e2x5MaxRatio);
    EventTree->SetBranchAddress("pfCluster_e2x5TopRatio", &pfCluster_e2x5TopRatio);//_pfCluster_e2x5TopRatio);
    EventTree->SetBranchAddress("pfCluster_e2x5RightRatio", &pfCluster_e2x5RightRatio);//_pfCluster_e2x5RightRatio);
    EventTree->SetBranchAddress("pfCluster_e2x5BottomRatio", &pfCluster_e2x5BottomRatio);//_pfCluster_e2x5BottomRatio);
    EventTree->SetBranchAddress("pfCluster_e2x5LeftRatio", &pfCluster_e2x5LeftRatio);//_pfCluster_e2x5LeftRatio);
    EventTree->SetBranchAddress("pfCluster_swissCross", &pfCluster_swissCross);//_pfCluster_swissCross);
    EventTree->SetBranchAddress("pfCluster_r9", &pfCluster_r9);//_pfCluster_r9);
    EventTree->SetBranchAddress("pfCluster_sigmaIetaIeta", &pfCluster_sigmaIetaIeta);//_pfCluster_sigmaIetaIeta);
    EventTree->SetBranchAddress("pfCluster_sigmaIetaIphi", &pfCluster_sigmaIetaIphi);//_pfCluster_sigmaIetaIphi);
    EventTree->SetBranchAddress("pfCluster_sigmaIphiIphi", &pfCluster_sigmaIphiIphi);//_pfCluster_sigmaIphiIphi);
    EventTree->SetBranchAddress("pfCluster_full5x5_e5x5", &pfCluster_full5x5_e5x5);//_pfCluster_full5x5_e5x5);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x2Ratio", &pfCluster_full5x5_e2x2Ratio);//_pfCluster_full5x5_e2x2Ratio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e3x3Ratio", &pfCluster_full5x5_e3x3Ratio);//_pfCluster_full5x5_e3x3Ratio);
    EventTree->SetBranchAddress("pfCluster_full5x5_eMaxRatio", &pfCluster_full5x5_eMaxRatio);//_pfCluster_full5x5_eMaxRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2ndRatio", &pfCluster_full5x5_e2ndRatio);//_pfCluster_full5x5_e2ndRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_eTopRatio", &pfCluster_full5x5_eTopRatio);//_pfCluster_full5x5_eTopRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_eRightRatio", &pfCluster_full5x5_eRightRatio);//_pfCluster_full5x5_eRightRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_eBottomRatio", &pfCluster_full5x5_eBottomRatio);//_pfCluster_full5x5_eBottomRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_eLeftRatio", &pfCluster_full5x5_eLeftRatio);//_pfCluster_full5x5_eLeftRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x5MaxRatio", &pfCluster_full5x5_e2x5MaxRatio);//_pfCluster_full5x5_e2x5MaxRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x5TopRatio", &pfCluster_full5x5_e2x5TopRatio);//_pfCluster_full5x5_e2x5TopRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x5RightRatio", &pfCluster_full5x5_e2x5RightRatio);//_pfCluster_full5x5_e2x5RightRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x5BottomRatio", &pfCluster_full5x5_e2x5BottomRatio);//_pfCluster_full5x5_e2x5BottomRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_e2x5LeftRatio", &pfCluster_full5x5_e2x5LeftRatio);//_pfCluster_full5x5_e2x5LeftRatio);
    EventTree->SetBranchAddress("pfCluster_full5x5_swissCross", &pfCluster_full5x5_swissCross);//_pfCluster_full5x5_swissCross);
    EventTree->SetBranchAddress("pfCluster_full5x5_r9", &pfCluster_full5x5_r9);//_pfCluster_full5x5_r9);
    EventTree->SetBranchAddress("pfCluster_full5x5_sigmaIetaIeta", &pfCluster_full5x5_sigmaIetaIeta);//_pfCluster_full5x5_sigmaIetaIeta);
    EventTree->SetBranchAddress("pfCluster_full5x5_sigmaIetaIphi", &pfCluster_full5x5_sigmaIetaIphi);//_pfCluster_full5x5_sigmaIetaIphi);
    EventTree->SetBranchAddress("pfCluster_full5x5_sigmaIphiIphi", &pfCluster_full5x5_sigmaIphiIphi);//_pfCluster_full5x5_sigmaIphiIphi);
    EventTree->SetBranchAddress("superCluster_e5x5", &superCluster_e5x5);//_superCluster_e5x5);
    EventTree->SetBranchAddress("superCluster_e2x2Ratio", &superCluster_e2x2Ratio);//_superCluster_e2x2Ratio);
    EventTree->SetBranchAddress("superCluster_e3x3Ratio", &superCluster_e3x3Ratio);//_superCluster_e3x3Ratio);
    EventTree->SetBranchAddress("superCluster_eMaxRatio", &superCluster_eMaxRatio);//_superCluster_eMaxRatio);
    EventTree->SetBranchAddress("superCluster_e2ndRatio", &superCluster_e2ndRatio);//_superCluster_e2ndRatio);
    EventTree->SetBranchAddress("superCluster_eTopRatio", &superCluster_eTopRatio);//_superCluster_eTopRatio);
    EventTree->SetBranchAddress("superCluster_eRightRatio", &superCluster_eRightRatio);//_superCluster_eRightRatio);
    EventTree->SetBranchAddress("superCluster_eBottomRatio", &superCluster_eBottomRatio);//_superCluster_eBottomRatio);
    EventTree->SetBranchAddress("superCluster_eLeftRatio", &superCluster_eLeftRatio);//_superCluster_eLeftRatio);
    EventTree->SetBranchAddress("superCluster_e2x5MaxRatio", &superCluster_e2x5MaxRatio);//_superCluster_e2x5MaxRatio);
    EventTree->SetBranchAddress("superCluster_e2x5TopRatio", &superCluster_e2x5TopRatio);//_superCluster_e2x5TopRatio);
    EventTree->SetBranchAddress("superCluster_e2x5RightRatio", &superCluster_e2x5RightRatio);//_superCluster_e2x5RightRatio);
    EventTree->SetBranchAddress("superCluster_e2x5BottomRatio", &superCluster_e2x5BottomRatio);//_superCluster_e2x5BottomRatio);
    EventTree->SetBranchAddress("superCluster_e2x5LeftRatio", &superCluster_e2x5LeftRatio);//_superCluster_e2x5LeftRatio);
    EventTree->SetBranchAddress("superCluster_swissCross", &superCluster_swissCross);//_superCluster_swissCross);
    EventTree->SetBranchAddress("superCluster_r9", &superCluster_r9);//_superCluster_r9);
    EventTree->SetBranchAddress("superCluster_sigmaIetaIeta", &superCluster_sigmaIetaIeta);//_superCluster_sigmaIetaIeta);
    EventTree->SetBranchAddress("superCluster_sigmaIetaIphi", &superCluster_sigmaIetaIphi);//_superCluster_sigmaIetaIphi);
    EventTree->SetBranchAddress("superCluster_sigmaIphiIphi", &superCluster_sigmaIphiIphi);//_superCluster_sigmaIphiIphi);
    EventTree->SetBranchAddress("superCluster_full5x5_e5x5", &superCluster_full5x5_e5x5);//_superCluster_full5x5_e5x5);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x2Ratio", &superCluster_full5x5_e2x2Ratio);//_superCluster_full5x5_e2x2Ratio);
    EventTree->SetBranchAddress("superCluster_full5x5_e3x3Ratio", &superCluster_full5x5_e3x3Ratio);//_superCluster_full5x5_e3x3Ratio);
    EventTree->SetBranchAddress("superCluster_full5x5_eMaxRatio", &superCluster_full5x5_eMaxRatio);//_superCluster_full5x5_eMaxRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2ndRatio", &superCluster_full5x5_e2ndRatio);//_superCluster_full5x5_e2ndRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_eTopRatio", &superCluster_full5x5_eTopRatio);//_superCluster_full5x5_eTopRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_eRightRatio", &superCluster_full5x5_eRightRatio);//_superCluster_full5x5_eRightRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_eBottomRatio", &superCluster_full5x5_eBottomRatio);//_superCluster_full5x5_eBottomRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_eLeftRatio", &superCluster_full5x5_eLeftRatio);//_superCluster_full5x5_eLeftRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x5MaxRatio", &superCluster_full5x5_e2x5MaxRatio);//_superCluster_full5x5_e2x5MaxRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x5TopRatio", &superCluster_full5x5_e2x5TopRatio);//_superCluster_full5x5_e2x5TopRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x5RightRatio", &superCluster_full5x5_e2x5RightRatio);//_superCluster_full5x5_e2x5RightRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x5BottomRatio", &superCluster_full5x5_e2x5BottomRatio);//_superCluster_full5x5_e2x5BottomRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_e2x5LeftRatio", &superCluster_full5x5_e2x5LeftRatio);//_superCluster_full5x5_e2x5LeftRatio);
    EventTree->SetBranchAddress("superCluster_full5x5_swissCross", &superCluster_full5x5_swissCross);//_superCluster_full5x5_swissCross);
    EventTree->SetBranchAddress("superCluster_full5x5_r9", &superCluster_full5x5_r9);//_superCluster_full5x5_r9);
    EventTree->SetBranchAddress("superCluster_full5x5_sigmaIetaIeta", &superCluster_full5x5_sigmaIetaIeta);//_superCluster_full5x5_sigmaIetaIeta);
    EventTree->SetBranchAddress("superCluster_full5x5_sigmaIetaIphi", &superCluster_full5x5_sigmaIetaIphi);//_superCluster_full5x5_sigmaIetaIphi);
    EventTree->SetBranchAddress("superCluster_full5x5_sigmaIphiIphi", &superCluster_full5x5_sigmaIphiIphi);//_superCluster_full5x5_sigmaIphiIphi);
    EventTree->SetBranchAddress("superCluster_HoEraw", &superCluster_HoEraw);//_superCluster_HoEraw);
    EventTree->SetBranchAddress("superCluster_HoErawBC", &superCluster_HoErawBC);//_superCluster_HoErawBC);
    /*EventTree->SetBranchAddress("retunedSuperCluster_e5x5", &retunedSuperCluster_e5x5);//_retunedSuperCluster_e5x5);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x2Ratio", &retunedSuperCluster_e2x2Ratio);//_retunedSuperCluster_e2x2Ratio);
    EventTree->SetBranchAddress("retunedSuperCluster_e3x3Ratio", &retunedSuperCluster_e3x3Ratio);//_retunedSuperCluster_e3x3Ratio);
    EventTree->SetBranchAddress("retunedSuperCluster_eMaxRatio", &retunedSuperCluster_eMaxRatio);//_retunedSuperCluster_eMaxRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2ndRatio", &retunedSuperCluster_e2ndRatio);//_retunedSuperCluster_e2ndRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_eTopRatio", &retunedSuperCluster_eTopRatio);//_retunedSuperCluster_eTopRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_eRightRatio", &retunedSuperCluster_eRightRatio);//_retunedSuperCluster_eRightRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_eBottomRatio", &retunedSuperCluster_eBottomRatio);//_retunedSuperCluster_eBottomRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_eLeftRatio", &retunedSuperCluster_eLeftRatio);//_retunedSuperCluster_eLeftRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x5MaxRatio", &retunedSuperCluster_e2x5MaxRatio);//_retunedSuperCluster_e2x5MaxRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x5TopRatio", &retunedSuperCluster_e2x5TopRatio);//_retunedSuperCluster_e2x5TopRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x5RightRatio", &retunedSuperCluster_e2x5RightRatio);//_retunedSuperCluster_e2x5RightRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x5BottomRatio", &retunedSuperCluster_e2x5BottomRatio);//_retunedSuperCluster_e2x5BottomRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_e2x5LeftRatio", &retunedSuperCluster_e2x5LeftRatio);//_retunedSuperCluster_e2x5LeftRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_swissCross", &retunedSuperCluster_swissCross);//_retunedSuperCluster_swissCross);
    EventTree->SetBranchAddress("retunedSuperCluster_r9", &retunedSuperCluster_r9);//_retunedSuperCluster_r9);
    EventTree->SetBranchAddress("retunedSuperCluster_sigmaIetaIeta", &retunedSuperCluster_sigmaIetaIeta);//_retunedSuperCluster_sigmaIetaIeta);
    EventTree->SetBranchAddress("retunedSuperCluster_sigmaIetaIphi", &retunedSuperCluster_sigmaIetaIphi);//_retunedSuperCluster_sigmaIetaIphi);
    EventTree->SetBranchAddress("retunedSuperCluster_sigmaIphiIphi", &retunedSuperCluster_sigmaIphiIphi);//_retunedSuperCluster_sigmaIphiIphi);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e5x5", &retunedSuperCluster_full5x5_e5x5);//_retunedSuperCluster_full5x5_e5x5);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x2Ratio", &retunedSuperCluster_full5x5_e2x2Ratio);//_retunedSuperCluster_full5x5_e2x2Ratio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e3x3Ratio", &retunedSuperCluster_full5x5_e3x3Ratio);//_retunedSuperCluster_full5x5_e3x3Ratio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_eMaxRatio", &retunedSuperCluster_full5x5_eMaxRatio);//_retunedSuperCluster_full5x5_eMaxRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2ndRatio", &retunedSuperCluster_full5x5_e2ndRatio);//_retunedSuperCluster_full5x5_e2ndRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_eTopRatio", &retunedSuperCluster_full5x5_eTopRatio);//_retunedSuperCluster_full5x5_eTopRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_eRightRatio", &retunedSuperCluster_full5x5_eRightRatio);//_retunedSuperCluster_full5x5_eRightRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_eBottomRatio", &retunedSuperCluster_full5x5_eBottomRatio);//_retunedSuperCluster_full5x5_eBottomRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_eLeftRatio", &retunedSuperCluster_full5x5_eLeftRatio);//_retunedSuperCluster_full5x5_eLeftRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x5MaxRatio", &retunedSuperCluster_full5x5_e2x5MaxRatio);//_retunedSuperCluster_full5x5_e2x5MaxRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x5TopRatio", &retunedSuperCluster_full5x5_e2x5TopRatio);//_retunedSuperCluster_full5x5_e2x5TopRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x5RightRatio", &retunedSuperCluster_full5x5_e2x5RightRatio);//_retunedSuperCluster_full5x5_e2x5RightRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x5BottomRatio", &retunedSuperCluster_full5x5_e2x5BottomRatio);//_retunedSuperCluster_full5x5_e2x5BottomRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_e2x5LeftRatio", &retunedSuperCluster_full5x5_e2x5LeftRatio);//_retunedSuperCluster_full5x5_e2x5LeftRatio);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_swissCross", &retunedSuperCluster_full5x5_swissCross);//_retunedSuperCluster_full5x5_swissCross);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_r9", &retunedSuperCluster_full5x5_r9);//_retunedSuperCluster_full5x5_r9);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_sigmaIetaIeta", &retunedSuperCluster_full5x5_sigmaIetaIeta);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_sigmaIetaIphi", &retunedSuperCluster_full5x5_sigmaIetaIphi);
    EventTree->SetBranchAddress("retunedSuperCluster_full5x5_sigmaIphiIphi", &retunedSuperCluster_full5x5_sigmaIphiIphi);
    EventTree->SetBranchAddress("retunedSuperCluster_HoEraw", &retunedSuperCluster_HoEraw);
    EventTree->SetBranchAddress("retunedSuperCluster_HoErawBC", &retunedSuperCluster_HoErawBC);*/



    EvMax=EventTree->GetEntries();

    cout<<cur_time()<<"\tTree successfully processed!\n";

    EventLoop();
    infile->Close();

    return;
}

Double_t fitFunc_upper(Double_t *x, Double_t *par){
//upper parabola fit function definition
    float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = par[4],
          seedeta = par[5];

    float p00 = -0.107537;
	float p01 = 0.590969;
	float p02 = -0.076494;
	float p10 = -0.0268843;
	float p11 = 0.147742;
	float p12 = -0.0191235;

    float c_upper = (p00 * pow(seedeta*sin(seedeta), 2)) + (p01 * seedeta*sin(seedeta)) + p02;

    float d_upper =( w10*seedeta*sin(seedeta) ) + (w11 / sqrt(1.1 + loget));
    float d_lower =( w00*seedeta*sin(seedeta) ) + (w01 / sqrt(1.1 + loget));
    float b_upper = d_upper - 0.5*(d_lower + d_upper);

    float a_upper = ((1 / (4 * c_upper))) - fabs(b_upper);

    Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f))*(x[0]*x[0]) + b_upper + 0.0087;

    return upper_curve;
}

Double_t fitFunc_lower(Double_t *x, Double_t *par){
//upper parabola fit function definition
    float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = par[4],
          seedeta = par[5];

    float p00 = -0.107537;
	float p01 = 0.590969;
	float p02 = -0.076494;
	float p10 = -0.0268843;
	float p11 = 0.147742;
	float p12 = -0.0191235;

    float c_lower = (p10 * pow(seedeta*sin(seedeta), 2)) + (p11 * seedeta*sin(seedeta)) + p12;

    float d_upper =( w10*seedeta*sin(seedeta) ) + (w11 / sqrt(1.1 + loget));
    float d_lower =( w00*seedeta*sin(seedeta) ) + (w01 / sqrt(1.1 + loget));
    float b_lower = d_lower - 0.5*(d_lower + d_upper);

    float a_lower = ((1 / (4 * c_lower))) - fabs(b_lower);

    Double_t lower_curve = (std::max((1 / (4 * a_lower)),0.0f))*(x[0]*x[0]) + b_lower;

    return lower_curve;
}

float getXError(int xBin, int yBin, int etaBin, int logetBin){
    float sum = caloClusters_shape_xAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    float total_entries = caloClusters_shape_xAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    //float sum = accumulate(SIPIP_seedEtaBin_logetBin_xBin.at(etaBin).at(logetBin).at(xBin).begin(),SIPIP_seedEtaBin_logetBin_xBin.at(etaBin).at(logetBin).at(xBin).end(),0.0);
    //float total_entries = SIPIP_seedEtaBin_logetBin_xBin.at(etaBin).at(logetBin).at(xBin).size();
    float weight = sum / total_entries;

    //cout<<"SUM\t" <<sum<<"\t"<<total_entries<<"\t"<<weight<<endl;

    return weight;
}

float getYError(int xBin, int yBin, int etaBin, int logetBin){
    float sum = caloClusters_shape_yAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    float total_entries = caloClusters_shape_yAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    //float sum = accumulate(SIEIE_seedEtaBin_logetBin_yBin.at(etaBin).at(logetBin).at(yBin).begin(),SIEIE_seedEtaBin_logetBin_yBin.at(etaBin).at(logetBin).at(yBin).end(),0.0);
    //float total_entries = SIEIE_seedEtaBin_logetBin_yBin.at(etaBin).at(logetBin).at(yBin).size();
    float weight = sum / total_entries;

    return weight;
}

void separation_opt(bool fit){
//optimization of parabola separation parameters
    vector<float> mid_etas = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9};
    vector<float> mid_logets = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9};

    TGraphErrors *envelope_plot_upper[15][15];
    TF1 *curve_fits_upper[15][15];

    TGraphErrors *envelope_plot_lower[15][15];
    TF1 *curve_fits_lower[15][15];

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    //loop over seed eta bins - 15    
    for(int seed_eta_idx = 5; seed_eta_idx  <6; seed_eta_idx++){ 
        //loop over loget bins - 15
        for(int loget_idx = 0; loget_idx < 15; loget_idx++){  
            int envelope_percent = 98,
                x_bins_total = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetNbinsX(),
                y_bins_total = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetNbinsY();

            int point_counter = 0;

            vector<float> x_points_upper, x_points_lower, y_points_upper, y_points_lower, 
                          x_errors_upper, y_errors_upper, x_errors_lower, y_errors_lower;

            c1->Clear();

            //loop over x bins
            for(int xbin_idx = 1; xbin_idx < x_bins_total+1; xbin_idx++){
                float xbin_integral = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);
                int n_entries = caloClusters_shape[seed_eta_idx][loget_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);

                if(n_entries < 5) continue;

                int ybin_idx_upper = -1, ybin_idx_lower = -1;
                float sum_temp_upper = 0, sum_temp_lower = 0;

                float xbin_percentage = ((float)envelope_percent * xbin_integral) / 100.0;

                if(xbin_integral <= 0.0) continue;

                float xBin_val = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetXaxis() -> GetBinCenter(xbin_idx);
                if(xBin_val > -0.05 && xBin_val < 0.05) continue; //skip the ring
                //if(xBin_val < -0.4 || xBin_val > 0.4) continue;

                //loop over y bins - upper
                for(int ybin_idx=1; ybin_idx<y_bins_total+1; ybin_idx++){
                    
                    float curr_bin_content = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_upper += curr_bin_content;

                    if(sum_temp_upper >= xbin_percentage){
                        ybin_idx_upper = ybin_idx;
                        break;
                    }
                } //y bin loop
                //loop over y bins - lower
                for(int ybin_idx = y_bins_total + 1; ybin_idx > 0; ybin_idx--){
                    //float ybin_integral 
                    float curr_bin_content = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_lower += curr_bin_content;

                    if(sum_temp_lower >= xbin_percentage){
                        ybin_idx_lower = ybin_idx;
                        break;
                    }
                } //y bin loop

                if (ybin_idx_upper != -1){
                    float envelope_x_upper = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                          envelope_y_upper = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_upper);

                    x_points_upper.push_back(envelope_x_upper);
                    x_errors_upper.push_back(getXError(xbin_idx, ybin_idx_upper, seed_eta_idx, loget_idx));
                    
                    y_points_upper.push_back(envelope_y_upper);
                    y_errors_upper.push_back(getYError(xbin_idx, ybin_idx_upper, seed_eta_idx, loget_idx));

                }
                if (ybin_idx_lower != -1){
                    float envelope_x_lower = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                          envelope_y_lower = caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_lower);

                    x_points_lower.push_back(envelope_x_lower);
                    x_errors_lower.push_back(getXError(xbin_idx, ybin_idx_lower, seed_eta_idx, loget_idx));

                    y_points_lower.push_back(envelope_y_lower);
                    y_errors_lower.push_back(getYError(xbin_idx, ybin_idx_lower, seed_eta_idx, loget_idx));

                }
            } //x bin loop

            //fill graphs for fitting
            envelope_plot_upper[seed_eta_idx][loget_idx] = new TGraphErrors(x_points_upper.size(), &x_points_upper[0], &y_points_upper[0], &x_errors_upper[0], &y_errors_upper[0]);//x_bins_total);
            envelope_plot_lower[seed_eta_idx][loget_idx] = new TGraphErrors(x_points_lower.size(), &x_points_lower[0], &y_points_lower[0], &x_errors_lower[0], &y_errors_lower[0]);//x_bins_total);

            if(envelope_plot_upper[seed_eta_idx][loget_idx]->GetN() <= 5) continue;

            //fitting
            if(fit){
                curve_fits_upper[seed_eta_idx][loget_idx] = new TF1(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str(),fitFunc_upper,-0.6,0.6,6);
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParameters(original_w00,original_w01,original_w10,original_w11,0.0,0.0);
                curve_fits_upper[seed_eta_idx][loget_idx] -> FixParameter(4, mid_logets.at(loget_idx));
                curve_fits_upper[seed_eta_idx][loget_idx] -> FixParameter(5, mid_etas.at(seed_eta_idx));

                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(0, "w00");
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(1, "w01");
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(2, "w10");
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(3, "w11");
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(4, "loget");
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(5, "seedeta");


                bool fitSuccess = false;
                double NDF = 0, 
                    chiSqperNDF = 999,
                    fitIteration = 0,
                    chiSq = 999;
                string fitStatus = "";

                fit_outfile<<"Eta: "<<seed_eta_idx<<" LogET: "<<loget_idx<<endl;

                //fit minimization
                while((!fitSuccess || chiSqperNDF > 1) && fitIteration < 20){
                    fitStatus.clear();
                    cout<<"UPPER FIT:\tEta: "<<seed_eta_idx<<"\tloget: "<<loget_idx<<"\titeration: "<<fitIteration<<endl;

                    envelope_plot_upper[seed_eta_idx][loget_idx] -> Fit(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str());

                    fitStatus = gMinuit->fCstatu;
                    fitStatus = trimString(fitStatus);
                    chiSq = curve_fits_upper[seed_eta_idx][loget_idx]->GetChisquare();
                    NDF = curve_fits_upper[seed_eta_idx][loget_idx]->GetNDF();

                    if(NDF != 0) chiSqperNDF = chiSq / NDF;

                    if((fitStatus.compare("CONVERGED") == 0) || (fitStatus.compare("OK") == 0) || (fitStatus.compare("SUCCESSFUL") == 0)) 
                    { fitSuccess = true; }

                    fitIteration++;
                }
                fit_outfile<<"UPPER FIT:\tIterations: "<<fitIteration<<" Status: "<<fitStatus<<endl;
                

                curve_fits_lower[seed_eta_idx][loget_idx] = new TF1(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str(),fitFunc_lower,-0.6,0.6,6);
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParameters(original_w00,original_w01,original_w10,original_w11,0.0,0.0);
                curve_fits_lower[seed_eta_idx][loget_idx] -> FixParameter(4,mid_logets.at(loget_idx));
                curve_fits_lower[seed_eta_idx][loget_idx] -> FixParameter(5,mid_etas.at(seed_eta_idx));

                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(0, "w00");
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(1, "w01");
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(2, "w10");
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(3, "w11");
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(4, "loget");
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(5, "seedeta");



                fitSuccess = false;
                NDF = 0;
                chiSqperNDF = 999;
                fitIteration = 0;
                chiSq = 999;
                fitStatus = "";

                //fit minimization
                while((!fitSuccess || chiSqperNDF > 1) && fitIteration < 20){
                    fitStatus.clear();
                    cout<<"LOWER FIT:\tEta: "<<seed_eta_idx<<"\tloget: "<<loget_idx<<"\titeration: "<<fitIteration<<endl;
                    envelope_plot_lower[seed_eta_idx][loget_idx]->Fit(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str());
            
                    fitStatus = gMinuit->fCstatu;
                    fitStatus = trimString(fitStatus);
                    chiSq = curve_fits_lower[seed_eta_idx][loget_idx]->GetChisquare();
                    NDF = curve_fits_lower[seed_eta_idx][loget_idx]->GetNDF();

                    if(NDF != 0) chiSqperNDF = chiSq / NDF;

                    if((fitStatus.compare("CONVERGED") == 0) || (fitStatus.compare("OK") == 0) || (fitStatus.compare("SUCCESSFUL") == 0)) 
                    { fitSuccess = true; }

                    fitIteration++;
                }
                fit_outfile<<"LOWER FIT:\tIterations: "<<fitIteration<<" Status: "<<fitStatus<<endl<<endl;
            }

            //plotting
            gStyle->SetOptStat(0);

            caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> SetTitle((titles_etas[seed_eta_idx] + "     " + titles_loget[loget_idx]).c_str());
            caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetXaxis() -> SetTitle("dPhi");
            caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> GetYaxis() -> SetTitle("dEta");

            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerStyle(34);  
            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerColor(2);  //red
            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerSize(3);

            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerStyle(47);
            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerColor(6);  //pink
            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerSize(3);

            caloClusters_shape_etWeight[seed_eta_idx][loget_idx] -> Draw("COLZ");
            envelope_plot_upper[seed_eta_idx][loget_idx] -> Draw("P, sames");
            envelope_plot_lower[seed_eta_idx][loget_idx] -> Draw("P, sames");

            if(fit){
                curve_fits_upper[seed_eta_idx][loget_idx] -> SetLineWidth(4);
                curve_fits_lower[seed_eta_idx][loget_idx] -> SetLineWidth(4);

                curve_fits_lower[seed_eta_idx][loget_idx] -> SetLineColor(6);  //pink
                curve_fits_upper[seed_eta_idx][loget_idx] -> Draw("SAME");
                curve_fits_lower[seed_eta_idx][loget_idx] -> Draw("SAME");
            }

            c1->SaveAs(("envelope/envelope_fit_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_etas[loget_idx]+".png").c_str());
            c1->SaveAs(("envelope/envelope_fit_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_etas[loget_idx]+".pdf").c_str());

            //aggregate parameter values
            if(fit){
                fit_w00_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(0));
                fit_w01_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(1));
                fit_w10_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(2));
                fit_w11_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(3));

                fit_w00_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(0));
                fit_w01_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(1));
                fit_w10_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(2));
                fit_w11_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(3));
            }

        } //loget index
    } //seed eta index
}

Double_t fitFunc_dPhi(Double_t *x, Double_t *par){
//dynamix delta phi fit function definition
    float yoffset = par[0],
          scale   = par[1],
          xoffset = par[2],
          width   = par[3];

    //float width = 0;

    //float cutoff = 0.14,
    //      saturation = 0.60;

    width = 1.0 / width;

    Double_t maxdphi = yoffset + scale / (1 + std::exp((x[0] - xoffset) * width));

    return maxdphi;
}

void dPhi_opt(){
//optimization of dynamic dPhi window parameters
    TGraph *dphi_envelope_plot[4];
    TF1 *dphi_fits[4];
    TF1 *dphi_original[4];

    TLegend* dphi_legends[4]; 

    TCanvas *c2 = new TCanvas("c2","c2",3600,2400);

    float xrange_min = -0.8,
          xrange_max = 1.5;

    //loop over dphi function areas
    for(int dphi_idx = 0; dphi_idx<4; dphi_idx++){

        int envelope_percent = 90,
            x_bins_total     = cluster_dPhi_vs_loget[dphi_idx] -> GetNbinsX(),
            y_bins_total     = cluster_dPhi_vs_loget[dphi_idx] -> GetNbinsY();

        int point_counter = 0;

        c2->Clear();
        dphi_envelope_plot[dphi_idx] = new TGraph();//x_bins_total);

        //loop over x bins
        for(int xbin_idx=0; xbin_idx<x_bins_total+1; xbin_idx++){
            float xbin_integral = cluster_dPhi_vs_loget[dphi_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);

            float envelope_biny = 0, 
                  sum_temp = 0;

            float xbin_percentage = ((float)envelope_percent * xbin_integral) / 100.0;

            //skip if the bin has no contents
            if(xbin_integral <= 0.0) continue;

            //loop over y bins
            for(int ybin_idx=0; ybin_idx<y_bins_total+1; ybin_idx++){
                //float ybin_integral 
                float curr_bin_content = cluster_dPhi_vs_loget[dphi_idx] -> GetBinContent(xbin_idx, ybin_idx);
                sum_temp += curr_bin_content;

                if(sum_temp >= xbin_percentage){
                    envelope_biny = ybin_idx;
                    break;
                }
            } //y bin loop
            float envelope_x = cluster_dPhi_vs_loget[dphi_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                  envelope_y = cluster_dPhi_vs_loget[dphi_idx] -> GetYaxis() -> GetBinCenter(envelope_biny);

            //set the point
            dphi_envelope_plot[dphi_idx] -> SetPoint(point_counter++, envelope_x, envelope_y);
        } //x bin loop

        //fitting
        dphi_fits[dphi_idx] = new TF1(("dphiFit"+to_string(dphi_idx)).c_str(),fitFunc_dPhi,xrange_min,xrange_max,4);
        //dphi_fits[dphi_idx] -> SetParameters(yoffset_orig[dphi_idx], scale_orig[dphi_idx], xoffset_orig[dphi_idx], width_orig[dphi_idx]);

        dphi_fits[dphi_idx] -> SetParName(0, "yoffset");
        dphi_fits[dphi_idx] -> SetParName(1, "scale");
        dphi_fits[dphi_idx] -> SetParName(2, "xoffset");
        dphi_fits[dphi_idx] -> SetParName(3, "width");

        bool fitSuccess = false;
        double NDF = 0, 
               chiSqperNDF = 999,
               fitIteration = 0,
               chiSq = 999;
        string fitStatus = "";

        //fit minimization
        while(fitIteration < 5 || chiSqperNDF > 1){
            dphi_envelope_plot[dphi_idx] -> Fit(("dphiFit"+to_string(dphi_idx)).c_str(), "RBSEM", "", xrange_min, xrange_max);

            fitStatus = gMinuit->fCstatu;
            chiSq = dphi_fits[dphi_idx]->GetChisquare();
            NDF = dphi_fits[dphi_idx]->GetNDF();

            if(NDF != 0) chiSqperNDF = chiSq / NDF;

            //cout << fitStatus <<endl;
            //if(fitStatus.compare("CONVERGED") || fitStatus.compare("OK") || fitStatus.compare("SUCCESSFUL")) fitSuccess = true;

            fitIteration++;
        }


        dphi_original[dphi_idx] = new TF1(("dphiOriginal"+to_string(dphi_idx)).c_str(),fitFunc_dPhi,xrange_min,xrange_max,4);
        dphi_original[dphi_idx] -> SetParameters(yoffset_orig[dphi_idx], scale_orig[dphi_idx], xoffset_orig[dphi_idx], width_orig[dphi_idx]);

        //plotting
        gStyle -> SetOptStat(0);

        dphi_legends[dphi_idx] = new TLegend(1.5,0.7,2.5,0.9,"","");

        cluster_dPhi_vs_loget[dphi_idx] -> GetXaxis() -> SetTitle("loget");
        cluster_dPhi_vs_loget[dphi_idx] -> GetYaxis() -> SetTitle("|dPhi|");

        dphi_envelope_plot[dphi_idx] -> GetXaxis() -> SetTitle("loget");
        dphi_envelope_plot[dphi_idx] -> GetYaxis() -> SetTitle("|dPhi|");

        dphi_envelope_plot[dphi_idx] -> SetMarkerStyle(34);
        dphi_envelope_plot[dphi_idx] -> SetMarkerColor(2);
        dphi_envelope_plot[dphi_idx] -> SetMarkerSize(3);
        dphi_envelope_plot[dphi_idx] -> SetLineWidth(6);

        dphi_fits[dphi_idx] -> SetLineWidth(6);

        dphi_original[dphi_idx] -> SetLineColor(7);
        dphi_original[dphi_idx] -> SetLineWidth(6);

        dphi_legends[dphi_idx] -> AddEntry(dphi_fits[dphi_idx], "Fit DPhi Function","l");
        dphi_legends[dphi_idx] -> AddEntry(dphi_original[dphi_idx], "Original DPhi Function","l");

        cluster_dPhi_vs_loget[dphi_idx] -> Scale(1.0 / cluster_dPhi_vs_loget[dphi_idx] -> Integral());

        cluster_dPhi_vs_loget[dphi_idx] -> Draw("COLZ");
        dphi_envelope_plot[dphi_idx] -> Draw("PSAME");
        dphi_fits[dphi_idx] -> Draw("SAME");
        dphi_original[dphi_idx] -> Draw("SAME");
        dphi_legends[dphi_idx]->Draw("SAME");

        c2 -> SaveAs(("dphi/dphi_fit_"+to_string(dphi_idx)+".png").c_str());

        yoffset_opt.push_back(dphi_fits[dphi_idx]->GetParameter(0));
        scale_opt.push_back(dphi_fits[dphi_idx]->GetParameter(1));
        xoffset_opt.push_back(dphi_fits[dphi_idx]->GetParameter(2));
        width_opt.push_back(dphi_fits[dphi_idx]->GetParameter(3));

    } //dphi index
}

void final_params(){
//aggregate the values of the parabola separation fits
    float w00_upper = (accumulate(fit_w00_upper.begin(),fit_w00_upper.end(),0.0)) / (float)fit_w00_upper.size(),
          w01_upper = (accumulate(fit_w01_upper.begin(),fit_w01_upper.end(),0.0)) / (float)fit_w01_upper.size(),
          w10_upper = (accumulate(fit_w10_upper.begin(),fit_w10_upper.end(),0.0)) / (float)fit_w10_upper.size(),
          w11_upper = (accumulate(fit_w11_upper.begin(),fit_w11_upper.end(),0.0)) / (float)fit_w11_upper.size(),
          w00_lower = (accumulate(fit_w00_lower.begin(),fit_w00_lower.end(),0.0)) / (float)fit_w00_lower.size(),
          w01_lower = (accumulate(fit_w01_lower.begin(),fit_w01_lower.end(),0.0)) / (float)fit_w01_lower.size(),
          w10_lower = (accumulate(fit_w10_lower.begin(),fit_w10_lower.end(),0.0)) / (float)fit_w10_lower.size(),
          w11_lower = (accumulate(fit_w11_lower.begin(),fit_w11_lower.end(),0.0)) / (float)fit_w11_lower.size();

    cout<<cur_time()<<"\tNumber of fits: "<<fit_w00_upper.size()<<endl;
    cout<<cur_time()<<"\n\tUpper Fit Results:\n\t\tw00: " << w00_upper << "\n\t\tw01: "<< w01_upper<< "\n\t\tw10: "<< w10_upper<< "\n\t\tw11: "<< w11_upper <<endl;
    cout<<cur_time()<<"\n\tLower Fit Results:\n\t\tw00: " << w00_lower << "\n\t\tw01: "<< w01_lower<< "\n\t\tw10: "<< w10_lower<< "\n\t\tw11: "<< w11_lower <<endl;


    fit_outfile<<"\n\n\nNumber of fits: "<<fit_w00_upper.size()<<endl;
    fit_outfile<<"\nUpper Fit Results:\n\t\tw00: " << w00_upper << "\n\t\tw01: "<< w01_upper<< "\n\t\tw10: "<< w10_upper<< "\n\t\tw11: "<< w11_upper <<endl;
    fit_outfile<<"\nLower Fit Results:\n\t\tw00: " << w00_lower << "\n\t\tw01: "<< w01_lower<< "\n\t\tw10: "<< w10_lower<< "\n\t\tw11: "<< w11_lower <<endl;
}

void final_params_dPhi(){
//aggregate the values of the dphi fits

    vector<string> yoffset_string = {"constexpr double yoffsetEB = ", "constexpr double yoffsetEE_0 = ","constexpr double yoffsetEE_1 =","constexpr double yoffsetEE_2 ="};
    vector<string> xoffset_string = {"constexpr double xoffsetEB = ", "constexpr double xoffsetEE_0 = ","constexpr double xoffsetEE_1 =","constexpr double xoffsetEE_2 ="};
    vector<string> scale_string = {"constexpr double scaleEB = ", "constexpr double scaleEE_0 = ","constexpr double scaleEE_1 =","constexpr double scaleEE_2 ="};
    vector<string> width_string = {"constexpr double widthEB = ", "constexpr double widthEE_0 = ","constexpr double widthEE_1 =","constexpr double widthEE_2 ="};

    cout<< "\n\n\n";
    for(int i = 0; i<4; i++){
        cout<< yoffset_string[i] << yoffset_opt[i] <<";"<<endl;
        cout<< scale_string[i] << scale_opt[i] <<";"<<endl;
        cout<< xoffset_string[i] << xoffset_opt[i] <<";"<<endl;
        cout<< width_string[i] << width_opt[i] <<";"<<endl<<endl;
    }

}

void envelope_optimization(){
//main program
    //InitTree("root://eoscms.cern.ch///eos/cms/store/user/lzygala/par_con_files/finalDumperTree_4G_Run3_PU.root");

    InitHistograms();
    //InitTree("root://eoscms.cern.ch///eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourGammasGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_Total.root");
    //InitTree("root://eoscms.cern.ch///eos/cms/store/user/lzygala/4Gamma_sample/4G_dumper_partial.root");


    const char* inDir;
    char* dir_xroot;
    char* dir;
    void* dirp;

    inDir = "/eos/cms/store/user/lzygala/4Gamma_sample/output_local";
    dir_xroot = "root://eoscms.cern.ch///eos/cms/store/user/lzygala/4Gamma_sample/output_local";

    dir = gSystem->ExpandPathName(inDir);
    dirp = gSystem->OpenDirectory(dir);

    const char* ext = ".root";

    const char* entry;
    const char* filename[10000];
    Int_t n = 0;
    TString str;

    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        str = entry;
        if(str.EndsWith(ext))
        filename[n++] = gSystem->ConcatFileName(dir_xroot, entry);
    }

    for (Int_t i = 0; i < n; i++){ //n
        //Printf("\tfile -> %s", filename[i]);
        Printf("\n%s\tfile -> %i / %i\t%s", (cur_time()).c_str(), i, n, filename[i]);
        InitTree(filename[i]);
        //EventLoop();
    }

    separation_opt(1);
    //dPhi_opt();
    //final_params_dPhi();
    final_params();

    fit_outfile.close();
}