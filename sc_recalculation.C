
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
#include<TLegend.h>
#include<TLine.h>
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

/* ---------- Variables ---------- */
    float pi=3.1415927;
    float twopi= 2*pi;

    float w00_opt = -0.00016325,
          w01_opt = 0.00330476,
          w10_opt = 0.00663276,
          w11_opt = -0.00461881;

    //current values
    float w00_orig = -0.00571429,
          w01_orig = -0.002,
          w10_orig = 0.0135714,
          w11_orig = 0.001;

    
    //original dPhi values
    vector<float> yoffset_orig = {0.07151, 0.05058, -0.0913, -0.6246},
                  scale_orig   = {0.5656,  0.7131,  44.04,   13.17},
                  xoffset_orig = {0.2931,  0.01668, -5.326,  -7.037},
                  width_orig   = {0.2976,  0.4114,  1.184,   2.836};

    
    //original dPhi values
    vector<float> yoffset_opt = {-0.0374196, 0.083572, -0.101296, -0.42682},
                  scale_opt  =  {2.53626,    0.656847, 32.2549,   10.9146},
                  xoffset_opt = {-1.18163,   0.0207699, -5.5264,  -7.44059},
                  width_opt   = {0.782166,   0.298995,  1.29536,   2.85026};
/*
    constexpr double yoffsetEB = -0.0374196;
    constexpr double scaleEB = 2.53626;
    constexpr double xoffsetEB = -1.18163;
    constexpr double widthEB = 0.782166];

    constexpr double yoffsetEE_0 = 0.083572;
    constexpr double scaleEE_0 = 0.656847;
    constexpr double xoffsetEE_0 = 0.0207699;
    constexpr double widthEE_0 = 0.298995;

    constexpr double yoffsetEE_1 =-0.101296;
    constexpr double scaleEE_1 = 32.2549;
    constexpr double xoffsetEE_1 = -5.5264;
    constexpr double widthEE_1 = 1.29536;

    constexpr double yoffsetEE_2 =-0.42682;
    constexpr double scaleEE_2 = 10.9146;
    constexpr double xoffsetEE_2 = -7.44059;
    constexpr double widthEE_2 = 2.85026; */

    vector<float> list_w00, list_w01, list_w10, list_w11, list_eff;

    //plotting
    TH2F* reformed_mustaches[7][12];
    TH1F* effInt_caloInSC_reformed;
    TH1F* effInt_totalCalo;
    TH1F* effInt_totalSC_reformed;
    TH1F* effInt_caloInSC_original;
    TH1F* effInt_totalSC_original;

    string titles_etas[12] = {"0 < |#eta| #leq 0.25", "0.25 < |#eta| #leq 0.45", "0.45 < |#eta| #leq 0.65", "0.65 < |#eta| #leq 0.85", "0.85 < |#eta| #leq 1.05", "1.05 < |#eta| #leq 1.25", "1.25 < |#eta| #leq 1.45", "1.45 < |#eta| #leq 1.65", "1.65 < |#eta| #leq 1.85", "1.85 < |#eta| #leq 2.05", "2.05 < |#eta| #leq 2.25", "2.25 < |#eta| #leq 2.50"};
    string titles_etas_dPhi[12] = {"0 #leq |#eta| < 0.2", "0.2 #leq |#eta| < 0.4", "0.4 #leq |#eta| < 0.6", "0.6 #leq |#eta| < 0.8", "0.8 #leq |#eta| < 1.0", "1.0 #leq |#eta| < 1.2", "1.2 #leq |#eta| < 1.4", "1.4 #leq |#eta| < 1.479","1.479 #leq |#eta| < 1.75", "1.75 #leq |#eta| < 2.0", "2.0 #leq |#eta| < 2.25", "2.25 #geq |#eta|"};
    string et_bin_titles[7] = {"All E_{T}","0 < E_{T} < 0.5", "0.5 #leq E_{T} < 1", "1 #leq E_{T} < 2", "2 #leq E_{T} < 3", "3 #leq E_{T} < 6", "E_{T} #geq 6"};

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);
    float sc_reformed_maxET[7][12],
          sc_reformed_minET[7][12],
          sc_original_maxET[7][12],
          sc_original_minET[7][12],
          sc_reformed_maxET_seedEta[7][12],
          sc_reformed_minET_seedEta[7][12];

/* Tree Values */
    Int_t N_ECALClusters;
    Int_t nVtx;

    std::vector<int> *genParticle_id;
    std::vector<float> *genParticle_energy;
    std::vector<float> *genParticle_pt;
    std::vector<float> *genParticle_eta;
    std::vector<float> *genParticle_phi;

    std::vector<std::vector<int>> *genParticle_pfCluster_dR_genScore_MatchedIndex;
    std::vector<std::vector<int>> *genParticle_superCluster_dR_genScore_MatchedIndex;


    std::vector<int> *caloParticle_id;
    std::vector<float> *caloParticle_genEnergy;
    std::vector<float> *caloParticle_simEnergy;
    std::vector<float> *caloParticle_genPt;
    std::vector<float> *caloParticle_simPt;
    std::vector<float> *caloParticle_genEta;
    std::vector<float> *caloParticle_simEta;
    std::vector<float> *caloParticle_genPhi;
    std::vector<float> *caloParticle_simPhi;
    std::vector<int> *caloParticle_simIeta;
    std::vector<int> *caloParticle_simIphi;
    std::vector<int> *caloParticle_simIz;

    std::vector<std::vector<int>> *caloParticle_pfCluster_dR_simScore_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_n_shared_xtals_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_sim_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_sim_fraction_min1_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_sim_fraction_min3_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_sim_rechit_diff_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex;

    std::vector<std::vector<int>> *caloParticle_superCluster_dR_simScore_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_n_shared_xtals_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_sim_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_sim_fraction_min1_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_sim_fraction_min3_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_sim_rechit_diff_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_sim_rechit_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex;
    std::vector<std::vector<int>> *caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex;

    std::vector<std::vector<float>> *simHit_energy;
    std::vector<std::vector<float>> *simHit_eta;
    std::vector<std::vector<float>> *simHit_phi;
    std::vector<std::vector<float>> *simHit_ieta;
    std::vector<std::vector<float>> *simHit_iphi;
    std::vector<std::vector<float>> *simHit_iz;
    
    std::vector<float> *recHit_noPF_energy;
    std::vector<float> *recHit_noPF_eta;
    std::vector<float> *recHit_noPF_phi;
    std::vector<int> *recHit_noPF_ieta;
    std::vector<int> *recHit_noPF_iphi;
    std::vector<int> *recHit_noPF_iz;

    std::vector<float> *pfRecHit_unClustered_energy;
    std::vector<float> *pfRecHit_unClustered_eta;
    std::vector<float> *pfRecHit_unClustered_phi;
    std::vector<int> *pfRecHit_unClustered_ieta;
    std::vector<int> *pfRecHit_unClustered_iphi;
    std::vector<int> *pfRecHit_unClustered_iz;

    std::vector<float> *pfCluster_energy;
    std::vector<float> *pfCluster_eta;
    std::vector<float> *pfCluster_phi;
    std::vector<int> *pfCluster_ieta;
    std::vector<int> *pfCluster_iphi;
    std::vector<int> *pfCluster_iz;

    std::vector<std::vector<int>> *pfCluster_superClustersIndex;

    std::vector<int> *pfCluster_dR_genScore_MatchedIndex;
    std::vector<int> *pfCluster_dR_simScore_MatchedIndex;
    std::vector<int> *pfCluster_n_shared_xtals_MatchedIndex;
    std::vector<int> *pfCluster_sim_fraction_MatchedIndex;
    std::vector<int> *pfCluster_sim_fraction_min1_MatchedIndex;
    std::vector<int> *pfCluster_sim_fraction_min3_MatchedIndex;
    std::vector<int> *pfCluster_sim_rechit_diff_MatchedIndex;
    std::vector<int> *pfCluster_sim_rechit_fraction_MatchedIndex;
    std::vector<int> *pfCluster_global_sim_rechit_fraction_MatchedIndex;
    std::vector<int> *pfCluster_hgcal_caloToCluster_MatchedIndex;
    std::vector<int> *pfCluster_hgcal_clusterToCalo_MatchedIndex;

    std::vector<std::vector<double>> *pfCluster_dR_genScore;
    std::vector<std::vector<double>> *pfCluster_dR_simScore;
    std::vector<std::vector<int>> *pfCluster_n_shared_xtals;
    std::vector<std::vector<double>> *pfCluster_sim_fraction;
    std::vector<std::vector<double>> *pfCluster_sim_fraction_min1;
    std::vector<std::vector<double>> *pfCluster_sim_fraction_min3;
    std::vector<std::vector<double>> *pfCluster_sim_rechit_diff;
    std::vector<std::vector<double>> *pfCluster_sim_rechit_fraction;
    std::vector<std::vector<double>> *pfCluster_global_sim_rechit_fraction;
    std::vector<std::vector<double>> *pfCluster_hgcal_caloToCluster;
    std::vector<std::vector<double>> *pfCluster_hgcal_clusterToCalo;

    std::vector<std::vector<float>> *pfClusterHit_energy;
    std::vector<std::vector<float>> *pfClusterHit_rechitEnergy;
    std::vector<std::vector<float>> *pfClusterHit_eta;
    std::vector<std::vector<float>> *pfClusterHit_phi;
    std::vector<std::vector<int>> *pfClusterHit_ieta;
    std::vector<std::vector<int>> *pfClusterHit_iphi;
    std::vector<std::vector<int>> *pfClusterHit_iz;

    std::vector<float> *superCluster_energy;
    std::vector<float> *superCluster_eta;
    std::vector<float> *superCluster_phi;
    std::vector<float> *superCluster_etaWidth;
    std::vector<float> *superCluster_phiWidth;
    std::vector<float> *superCluster_R;
    std::vector<int> *superCluster_ieta;
    std::vector<int> *superCluster_iphi;
    std::vector<int> *superCluster_iz;

    std::vector<int> *superCluster_seedIndex;
    std::vector<std::vector<int>> *superCluster_pfClustersIndex;

    std::vector<std::vector<float>> *psCluster_energy;
    std::vector<std::vector<float>> *psCluster_eta;
    std::vector<std::vector<float>> *psCluster_phi;

    std::vector<int> *superCluster_dR_genScore_MatchedIndex;
    std::vector<int> *superCluster_dR_simScore_MatchedIndex;
    std::vector<int> *superCluster_n_shared_xtals_MatchedIndex;
    std::vector<int> *superCluster_sim_fraction_MatchedIndex;
    std::vector<int> *superCluster_sim_fraction_min1_MatchedIndex;
    std::vector<int> *superCluster_sim_fraction_min3_MatchedIndex;
    std::vector<int> *superCluster_sim_rechit_diff_MatchedIndex;
    std::vector<int> *superCluster_sim_rechit_fraction_MatchedIndex;
    std::vector<int> *superCluster_global_sim_rechit_fraction_MatchedIndex;
    std::vector<int> *superCluster_hgcal_caloToCluster_MatchedIndex;
    std::vector<int> *superCluster_hgcal_clusterToCalo_MatchedIndex;

    std::vector<std::vector<double>> *superCluster_dR_genScore;
    std::vector<std::vector<double>> *superCluster_dR_simScore;
    std::vector<std::vector<int>> *superCluster_n_shared_xtals;
    std::vector<std::vector<double>> *superCluster_sim_fraction;
    std::vector<std::vector<double>> *superCluster_sim_fraction_min1;
    std::vector<std::vector<double>> *superCluster_sim_fraction_min3;
    std::vector<std::vector<double>> *superCluster_sim_rechit_diff;
    std::vector<std::vector<double>> *superCluster_sim_rechit_fraction;
    std::vector<std::vector<double>> *superCluster_global_sim_rechit_fraction;
    std::vector<std::vector<double>> *superCluster_hgcal_caloToCluster;
    std::vector<std::vector<double>> *superCluster_hgcal_clusterToCalo;

    TTree* EventTree;
    Int_t EvMax;
/* */

string cur_time(){
//returns the current time, for logging
    std::time_t tt = std::time(NULL);
    std::string s = std::ctime(&tt);
    return s.substr(0, s.size()-1);
}

bool DynamicDPhi_orig(float clust_ET, float absSeedEta, float dPhiCurrent){
//return true if cluster is within dPhi window, false if not
    constexpr double yoffsetEB = 7.151e-02;
    constexpr double scaleEB = 5.656e-01;
    constexpr double xoffsetEB = 2.931e-01;
    constexpr double widthEB = 2.976e-01;

    constexpr double yoffsetEE_0 = 5.058e-02;
    constexpr double scaleEE_0 = 7.131e-01;
    constexpr double xoffsetEE_0 = 1.668e-02;
    constexpr double widthEE_0 = 4.114e-01;

    constexpr double yoffsetEE_1 = -9.913e-02;
    constexpr double scaleEE_1 = 4.404e+01;
    constexpr double xoffsetEE_1 = -5.326e+00;
    constexpr double widthEE_1 = 1.184e+00;

    constexpr double yoffsetEE_2 = -6.346e-01;
    constexpr double scaleEE_2 = 1.317e+01;
    constexpr double xoffsetEE_2 = -7.037e+00;
    constexpr double widthEE_2 = 2.836e+00;

    const int etaBin = ((int)(absSeedEta >= 1.479) + (int)(absSeedEta >= 1.75) + (int)(absSeedEta >= 2.0));
    const double logClustEt = log10(clust_ET);
    double yoffset, scale, xoffset, width, saturation, cutoff, maxdphi;

    switch (etaBin) {
            case 0:  // EB
            yoffset = yoffsetEB;
            scale = scaleEB;
            xoffset = xoffsetEB;
            width = 1.0 / widthEB;
            saturation = 0.14;
            cutoff = 0.60;
            break;
            case 1:  // 1.479 -> 1.75
            yoffset = yoffsetEE_0;
            scale = scaleEE_0;
            xoffset = xoffsetEE_0;
            width = 1.0 / widthEE_0;
            saturation = 0.14;
            cutoff = 0.55;
            break;
            case 2:  // 1.75 -> 2.0
            yoffset = yoffsetEE_1;
            scale = scaleEE_1;
            xoffset = xoffsetEE_1;
            width = 1.0 / widthEE_1;
            saturation = 0.12;
            cutoff = 0.45;
            break;
            case 3:  // 2.0 and up
            yoffset = yoffsetEE_2;
            scale = scaleEE_2;
            xoffset = xoffsetEE_2;
            width = 1.0 / widthEE_2;
            saturation = 0.12;
            cutoff = 0.30;
            break;
    }

    maxdphi = yoffset + scale / (1 + std::exp((logClustEt - xoffset) * width));
    maxdphi = std::min(maxdphi, cutoff);
    maxdphi = std::max(maxdphi, saturation);

    return fabs(dPhiCurrent) < maxdphi;
}

bool DynamicDPhi_opt(float clust_ET, float absSeedEta, float dPhiCurrent){
//return true if cluster is within dPhi window, false if not
    /*constexpr double yoffsetEB = yoffset_opt[0];
    constexpr double scaleEB = scale_opt[0];
    constexpr double xoffsetEB = xoffset_opt[0];
    constexpr double widthEB = width_opt[0];

    constexpr double yoffsetEE_0 = yoffset_opt[1];
    constexpr double scaleEE_0 = scale_opt[1];
    constexpr double xoffsetEE_0 = xoffset_opt[1];
    constexpr double widthEE_0 = width_opt[1];

    constexpr double yoffsetEE_1 =yoffset_opt[2];
    constexpr double scaleEE_1 = scale_opt[2];
    constexpr double xoffsetEE_1 = xoffset_opt[2];
    constexpr double widthEE_1 = width_opt[2];

    constexpr double yoffsetEE_2 =yoffset_opt[3];
    constexpr double scaleEE_2 = scale_opt[3];
    constexpr double xoffsetEE_2 = xoffset_opt[3];
    constexpr double widthEE_2 = width_opt[3];*/
        constexpr double yoffsetEB = -0.0374196;
    constexpr double scaleEB = 2.53626;
    constexpr double xoffsetEB = -1.18163;
    constexpr double widthEB = 0.782166;

    constexpr double yoffsetEE_0 = 0.083572;
    constexpr double scaleEE_0 = 0.656847;
    constexpr double xoffsetEE_0 = 0.0207699;
    constexpr double widthEE_0 = 0.298995;

    constexpr double yoffsetEE_1 =-0.101296;
    constexpr double scaleEE_1 = 32.2549;
    constexpr double xoffsetEE_1 = -5.5264;
    constexpr double widthEE_1 = 1.29536;

    constexpr double yoffsetEE_2 =-0.42682;
    constexpr double scaleEE_2 = 10.9146;
    constexpr double xoffsetEE_2 = -7.44059;
    constexpr double widthEE_2 = 2.85026;

    const int etaBin = ((int)(absSeedEta >= 1.479) + (int)(absSeedEta >= 1.75) + (int)(absSeedEta >= 2.0));
    const double logClustEt = log10(clust_ET);
    double yoffset, scale, xoffset, width, saturation, cutoff, maxdphi;

    switch (etaBin) {
            case 0:  // EB
            yoffset = yoffsetEB;
            scale = scaleEB;
            xoffset = xoffsetEB;
            width = 1.0 / widthEB;
            saturation = 0.14;
            cutoff = 0.60;
            break;
            case 1:  // 1.479 -> 1.75
            yoffset = yoffsetEE_0;
            scale = scaleEE_0;
            xoffset = xoffsetEE_0;
            width = 1.0 / widthEE_0;
            saturation = 0.14;
            cutoff = 0.55;
            break;
            case 2:  // 1.75 -> 2.0
            yoffset = yoffsetEE_1;
            scale = scaleEE_1;
            xoffset = xoffsetEE_1;
            width = 1.0 / widthEE_1;
            saturation = 0.12;
            cutoff = 0.45;
            break;
            case 3:  // 2.0 and up
            yoffset = yoffsetEE_2;
            scale = scaleEE_2;
            xoffset = xoffsetEE_2;
            width = 1.0 / widthEE_2;
            saturation = 0.12;
            cutoff = 0.30;
            break;
    }

    maxdphi = yoffset + scale / (1 + std::exp((logClustEt - xoffset) * width));
    maxdphi = std::min(maxdphi, cutoff);
    maxdphi = std::max(maxdphi, saturation);

    return fabs(dPhiCurrent) < maxdphi;
}

float DynamicDPhi_orig_val(float clust_ET, float absSeedEta){
//return true if cluster is within dPhi window, false if not
    constexpr double yoffsetEB = 7.151e-02;
    constexpr double scaleEB = 5.656e-01;
    constexpr double xoffsetEB = 2.931e-01;
    constexpr double widthEB = 2.976e-01;

    constexpr double yoffsetEE_0 = 5.058e-02;
    constexpr double scaleEE_0 = 7.131e-01;
    constexpr double xoffsetEE_0 = 1.668e-02;
    constexpr double widthEE_0 = 4.114e-01;

    constexpr double yoffsetEE_1 = -9.913e-02;
    constexpr double scaleEE_1 = 4.404e+01;
    constexpr double xoffsetEE_1 = -5.326e+00;
    constexpr double widthEE_1 = 1.184e+00;

    constexpr double yoffsetEE_2 = -6.346e-01;
    constexpr double scaleEE_2 = 1.317e+01;
    constexpr double xoffsetEE_2 = -7.037e+00;
    constexpr double widthEE_2 = 2.836e+00;

    const int etaBin = ((int)(absSeedEta >= 1.479) + (int)(absSeedEta >= 1.75) + (int)(absSeedEta >= 2.0));
    const double logClustEt = log10(clust_ET);
    double yoffset, scale, xoffset, width, saturation, cutoff, maxdphi;

    switch (etaBin) {
            case 0:  // EB
            yoffset = yoffsetEB;
            scale = scaleEB;
            xoffset = xoffsetEB;
            width = 1.0 / widthEB;
            saturation = 0.14;
            cutoff = 0.60;
            break;
            case 1:  // 1.479 -> 1.75
            yoffset = yoffsetEE_0;
            scale = scaleEE_0;
            xoffset = xoffsetEE_0;
            width = 1.0 / widthEE_0;
            saturation = 0.14;
            cutoff = 0.55;
            break;
            case 2:  // 1.75 -> 2.0
            yoffset = yoffsetEE_1;
            scale = scaleEE_1;
            xoffset = xoffsetEE_1;
            width = 1.0 / widthEE_1;
            saturation = 0.12;
            cutoff = 0.45;
            break;
            case 3:  // 2.0 and up
            yoffset = yoffsetEE_2;
            scale = scaleEE_2;
            xoffset = xoffsetEE_2;
            width = 1.0 / widthEE_2;
            saturation = 0.12;
            cutoff = 0.30;
            break;
    }

    maxdphi = yoffset + scale / (1 + std::exp((logClustEt - xoffset) * width));
    maxdphi = std::min(maxdphi, cutoff);
    maxdphi = std::max(maxdphi, saturation);

    return maxdphi;
}

float DynamicDPhi_opt_val(float clust_ET, float absSeedEta){
//return true if cluster is within dPhi window, false if not
    /*constexpr double yoffsetEB = yoffset_opt[0];
    constexpr double scaleEB = scale_opt[0];
    constexpr double xoffsetEB = xoffset_opt[0];
    constexpr double widthEB = width_opt[0];

    constexpr double yoffsetEE_0 = yoffset_opt[1];
    constexpr double scaleEE_0 = scale_opt[1];
    constexpr double xoffsetEE_0 = xoffset_opt[1];
    constexpr double widthEE_0 = width_opt[1];

    constexpr double yoffsetEE_1 =yoffset_opt[2];
    constexpr double scaleEE_1 = scale_opt[2];
    constexpr double xoffsetEE_1 = xoffset_opt[2];
    constexpr double widthEE_1 = width_opt[2];

    constexpr double yoffsetEE_2 =yoffset_opt[3];
    constexpr double scaleEE_2 = scale_opt[3];
    constexpr double xoffsetEE_2 = xoffset_opt[3];
    constexpr double widthEE_2 = width_opt[3];*/
        constexpr double yoffsetEB = -0.0374196;
    constexpr double scaleEB = 2.53626;
    constexpr double xoffsetEB = -1.18163;
    constexpr double widthEB = 0.782166;

    constexpr double yoffsetEE_0 = 0.083572;
    constexpr double scaleEE_0 = 0.656847;
    constexpr double xoffsetEE_0 = 0.0207699;
    constexpr double widthEE_0 = 0.298995;

    constexpr double yoffsetEE_1 =-0.101296;
    constexpr double scaleEE_1 = 32.2549;
    constexpr double xoffsetEE_1 = -5.5264;
    constexpr double widthEE_1 = 1.29536;

    constexpr double yoffsetEE_2 =-0.42682;
    constexpr double scaleEE_2 = 10.9146;
    constexpr double xoffsetEE_2 = -7.44059;
    constexpr double widthEE_2 = 2.85026;

    const int etaBin = ((int)(absSeedEta >= 1.479) + (int)(absSeedEta >= 1.75) + (int)(absSeedEta >= 2.0));
    const double logClustEt = log10(clust_ET);
    double yoffset, scale, xoffset, width, saturation, cutoff, maxdphi;

    switch (etaBin) {
            case 0:  // EB
            yoffset = yoffsetEB;
            scale = scaleEB;
            xoffset = xoffsetEB;
            width = 1.0 / widthEB;
            saturation = 0.14;
            cutoff = 0.60;
            break;
            case 1:  // 1.479 -> 1.75
            yoffset = yoffsetEE_0;
            scale = scaleEE_0;
            xoffset = xoffsetEE_0;
            width = 1.0 / widthEE_0;
            saturation = 0.14;
            cutoff = 0.55;
            break;
            case 2:  // 1.75 -> 2.0
            yoffset = yoffsetEE_1;
            scale = scaleEE_1;
            xoffset = xoffsetEE_1;
            width = 1.0 / widthEE_1;
            saturation = 0.12;
            cutoff = 0.45;
            break;
            case 3:  // 2.0 and up
            yoffset = yoffsetEE_2;
            scale = scaleEE_2;
            xoffset = xoffsetEE_2;
            width = 1.0 / widthEE_2;
            saturation = 0.12;
            cutoff = 0.30;
            break;
    }

    maxdphi = yoffset + scale / (1 + std::exp((logClustEt - xoffset) * width));
    maxdphi = std::min(maxdphi, cutoff);
    maxdphi = std::max(maxdphi, saturation);

    return maxdphi;
}

bool upper_parabola(float eta_max, float clusET, float dEta, float dPhi, float w00, float w01, float w10, float w11){
//returns true if cluster is located under parabola
	float p00 = -0.107537;
	float p01 = 0.590969;
	float p02 = -0.076494;
	float p10 = -0.0268843;
	float p11 = 0.147742;
	float p12 = -0.0191235;

    float c_upper = (p00 * pow(eta_max*sin(eta_max), 2)) + (p01 * eta_max*sin(eta_max)) + p02;

    float d_upper =( w10*eta_max*sin(eta_max) ) + (w11 / sqrt(1.1 + log10(clusET)));
    float d_lower =( w00*eta_max*sin(eta_max) ) + (w01 / sqrt(1.1 + log10(clusET)));
    float b_upper = d_upper - 0.5*(d_lower + d_upper);

    float a_upper = ((1 / (4 * c_upper))) - fabs(b_upper);
    float upper_curve = (std::max((1 / (4 * a_upper)),0.0f)*pow(dPhi,2)) + std::max(b_upper, 0.0087f) + 0.0087;

    return dEta < upper_curve;
}

bool lower_parabola(float eta_max, float clusET, float dEta, float dPhi, float w00, float w01, float w10, float w11){
//returns true if cluster is located above parabola
    float p00 = -0.107537;
    float p01 = 0.590969;
    float p02 = -0.076494;
    float p10 = -0.0268843;
    float p11 = 0.147742;
    float p12 = -0.0191235;

    float c_lower = (p10 * pow(eta_max*sin(eta_max), 2)) + (p11 * eta_max*sin(eta_max)) + p12;

    float d_upper =( w10*eta_max*sin(eta_max) ) + (w11 / sqrt(1.1 + log10(clusET)));
    float d_lower =( w00*eta_max*sin(eta_max) ) + (w01 / sqrt(1.1 + log10(clusET)));
    float b_lower = d_lower - 0.5*(d_lower + d_upper);

    float a_lower = ((1 / (4 * c_lower))) - fabs(b_lower);
    float lower_curve = (std::max((1 / (4 * a_lower)),0.0f)*pow(dPhi,2)) + std::min(b_lower, -0.0087f);

    return dEta > lower_curve;
}

void InitTree(TString infileName){
//Gets all the information from the tree and fills the cluster vectors
    //outfile_efficiences = new TFile("efficiencies_2.root","UPDATE");

    //caloParticle_widths.resize

    /*double min_dEta = -0.15,
           max_dEta = 0.2,
           min_dPhi = -0.6,
           max_dPhi = 0.6;
    
    int bins_dEta = 105,
          bins_dPhi = 120;

    effInt_caloInSC = new TH1F("effInt_caloInSC","effInt_caloInSC",60,0,3);
    effInt_totalCalo = new TH1F("effInt_totalCalo","effInt_totalCalo",60,0,3);
    effInt_totalSC = new TH1F("effInt_totalSC","effInt_totalSC",60,0,3);

    for(int i=0;i<7;i++){
        for(int j=0;j<12;j++){
            mustaches_pileup_vs_seed[i][j] = new TH2F(("dEtaVsDphi_pileup_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("dEtaVsDphi_pileup_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            mustaches_caloClusters_vs_seed[i][j]  = new TH2F(("dEtaVsDphi_calo_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("dEtaVsDphi_calo_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            mustaches_OTHERcaloClusters_vs_seed[i][j]  = new TH2F(("dEtaVsDphi_OTHERcalo_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("dEtaVsDphi_OTHERcalo_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_inSC[i][j]  = new TH2F(("caloClusters_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("caloClusters_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_other_inSC[i][j]  = new TH2F(("caloClusters_other_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("caloClusters_other_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_total[i][j]  = new TH2F(("caloClusters_total_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("caloClusters_total_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            
            
            pileup_inSC[i][j]  = new TH2F(("pileup_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("pileup_inSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_NOTinSC[i][j]  = new TH2F(("caloClusters_NOTinSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("caloClusters_NOTinSC_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);

            totalCaloClusters_inSC[i][j] = 0;
            totalCaloClusters[i][j] = 0;
            totalClusters_inSC[i][j] = 0;
        }
    }

    ring_inSC = new TH2F("ring_inSC","ring_inSC",50,-0.25,0.25,50,-0.25,0.25);
    ring_NOTinSC = new TH2F("ring_inSC","ring_inSC",50,-0.25,0.25,50,-0.25,0.25);
    */
   double min_dEta = -0.15,
           max_dEta = 0.2,
           min_dPhi = -0.6,
           max_dPhi = 0.6;
    
    int bins_dEta = 105,
          bins_dPhi = 120;
    for(int i=0;i<7;i++){
        for(int j=0;j<12;j++){
            reformed_mustaches[i][j] = new TH2F(("reformed_mustaches_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),("reformed_mustaches_etBin_"+to_string(i)+"_etaBin_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            sc_reformed_maxET[i][j] = -999;
            sc_reformed_minET[i][j] = 999;
        }
    }
    effInt_caloInSC_reformed = new TH1F("effInt_caloInSC_reformed","effInt_caloInSC_reformed",60,-3,3);
    effInt_totalCalo = new TH1F("effInt_totalCalo","effInt_totalCalo",60,-3,3);
    effInt_totalSC_reformed = new TH1F("effInt_totalSC_reformed","effInt_totalSC_reformed",60,-3,3);

    effInt_caloInSC_original = new TH1F("effInt_caloInSC_original","effInt_caloInSC_original",60,-3,3);
    effInt_totalSC_original = new TH1F("effInt_totalSC_original","effInt_totalSC_original",60,-3,3);


    cout<<cur_time()<<"\tProcessing Tree...\n";
   //outfile<<cur_time()<<"\tProcessing Tree...\n";
    TFile* infile = new TFile(infileName);
    infile->cd("recosimdumper");
    EventTree=(TTree*)gDirectory->Get("caloTree");

	EventTree->SetBranchAddress("genParticle_id", &genParticle_id);
	EventTree->SetBranchAddress("genParticle_energy", &genParticle_energy);
	EventTree->SetBranchAddress("genParticle_pt", &genParticle_pt);
	EventTree->SetBranchAddress("genParticle_eta", &genParticle_eta);
	EventTree->SetBranchAddress("genParticle_phi", &genParticle_phi);

	EventTree->SetBranchAddress("genParticle_pfCluster_dR_genScore_MatchedIndex", &genParticle_pfCluster_dR_genScore_MatchedIndex);
	EventTree->SetBranchAddress("genParticle_superCluster_dR_genScore_MatchedIndex", &genParticle_superCluster_dR_genScore_MatchedIndex);

	EventTree->SetBranchAddress("caloParticle_id", &caloParticle_id);
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
	EventTree->SetBranchAddress("caloParticle_pfCluster_n_shared_xtals_MatchedIndex", &caloParticle_pfCluster_n_shared_xtals_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_sim_fraction_MatchedIndex", &caloParticle_pfCluster_sim_fraction_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_sim_fraction_min1_MatchedIndex", &caloParticle_pfCluster_sim_fraction_min1_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_sim_fraction_min3_MatchedIndex", &caloParticle_pfCluster_sim_fraction_min3_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_sim_rechit_diff_MatchedIndex", &caloParticle_pfCluster_sim_rechit_diff_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex", &caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex", &caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex", &caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex", &caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex);

	EventTree->SetBranchAddress("caloParticle_superCluster_dR_simScore_MatchedIndex", &caloParticle_superCluster_dR_simScore_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_n_shared_xtals_MatchedIndex", &caloParticle_superCluster_n_shared_xtals_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_sim_fraction_MatchedIndex", &caloParticle_superCluster_sim_fraction_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_sim_fraction_min1_MatchedIndex", &caloParticle_superCluster_sim_fraction_min1_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_sim_fraction_min3_MatchedIndex", &caloParticle_superCluster_sim_fraction_min3_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_sim_rechit_diff_MatchedIndex", &caloParticle_superCluster_sim_rechit_diff_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_sim_rechit_fraction_MatchedIndex", &caloParticle_superCluster_sim_rechit_fraction_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex", &caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex);
    EventTree->SetBranchAddress("caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex", &caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex);
	EventTree->SetBranchAddress("caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex", &caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex);

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
	EventTree->SetBranchAddress("recHit_noPF_iphi", &recHit_noPF_iphi);
	EventTree->SetBranchAddress("recHit_noPF_iz", &recHit_noPF_iz);

	EventTree->SetBranchAddress("pfRecHit_unClustered_energy", &pfRecHit_unClustered_energy);
	EventTree->SetBranchAddress("pfRecHit_unClustered_eta", &pfRecHit_unClustered_eta);
	EventTree->SetBranchAddress("pfRecHit_unClustered_phi", &pfRecHit_unClustered_phi);
	EventTree->SetBranchAddress("pfRecHit_unClustered_ieta", &pfRecHit_unClustered_ieta);
	EventTree->SetBranchAddress("pfRecHit_unClustered_iphi", &pfRecHit_unClustered_iphi);
	EventTree->SetBranchAddress("pfRecHit_unClustered_iz", &pfRecHit_unClustered_iz);

	EventTree->SetBranchAddress("pfCluster_energy", &pfCluster_energy);
	EventTree->SetBranchAddress("pfCluster_eta", &pfCluster_eta);
	EventTree->SetBranchAddress("pfCluster_phi", &pfCluster_phi);
	EventTree->SetBranchAddress("pfCluster_ieta", &pfCluster_ieta);
	EventTree->SetBranchAddress("pfCluster_iphi", &pfCluster_iphi);
	EventTree->SetBranchAddress("pfCluster_iz", &pfCluster_iz);

	EventTree->SetBranchAddress("pfCluster_superClustersIndex", &pfCluster_superClustersIndex);

	EventTree->SetBranchAddress("pfCluster_dR_genScore_MatchedIndex", &pfCluster_dR_genScore_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_dR_simScore_MatchedIndex", &pfCluster_dR_simScore_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_n_shared_xtals_MatchedIndex", &pfCluster_n_shared_xtals_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_sim_fraction_MatchedIndex", &pfCluster_sim_fraction_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_sim_fraction_min1_MatchedIndex", &pfCluster_sim_fraction_min1_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_sim_fraction_min3_MatchedIndex", &pfCluster_sim_fraction_min3_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_sim_rechit_diff_MatchedIndex", &pfCluster_sim_rechit_diff_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_sim_rechit_fraction_MatchedIndex", &pfCluster_sim_rechit_fraction_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_global_sim_rechit_fraction_MatchedIndex", &pfCluster_global_sim_rechit_fraction_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_hgcal_caloToCluster_MatchedIndex", &pfCluster_hgcal_caloToCluster_MatchedIndex);
	EventTree->SetBranchAddress("pfCluster_hgcal_clusterToCalo_MatchedIndex", &pfCluster_hgcal_clusterToCalo_MatchedIndex);

	EventTree->SetBranchAddress("pfCluster_dR_genScore", &pfCluster_dR_genScore);
	EventTree->SetBranchAddress("pfCluster_dR_simScore", &pfCluster_dR_simScore);
	EventTree->SetBranchAddress("pfCluster_n_shared_xtals", &pfCluster_n_shared_xtals);
	EventTree->SetBranchAddress("pfCluster_sim_fraction", &pfCluster_sim_fraction);
	EventTree->SetBranchAddress("pfCluster_sim_fraction_min1", &pfCluster_sim_fraction_min1);
	EventTree->SetBranchAddress("pfCluster_sim_fraction_min3", &pfCluster_sim_fraction_min3);
	EventTree->SetBranchAddress("pfCluster_sim_rechit_diff", &pfCluster_sim_rechit_diff);
	EventTree->SetBranchAddress("pfCluster_sim_rechit_fraction", &pfCluster_sim_rechit_fraction);
	EventTree->SetBranchAddress("pfCluster_global_sim_rechit_fraction", &pfCluster_global_sim_rechit_fraction);
    EventTree->SetBranchAddress("pfCluster_hgcal_caloToCluster", &pfCluster_hgcal_caloToCluster);
	EventTree->SetBranchAddress("pfCluster_hgcal_clusterToCalo", &pfCluster_hgcal_clusterToCalo);

	EventTree->SetBranchAddress("pfClusterHit_energy", &pfClusterHit_energy);
	EventTree->SetBranchAddress("pfClusterHit_rechitEnergy", &pfClusterHit_rechitEnergy);
	EventTree->SetBranchAddress("pfClusterHit_eta", &pfClusterHit_eta);
	EventTree->SetBranchAddress("pfClusterHit_phi", &pfClusterHit_phi);
	EventTree->SetBranchAddress("pfClusterHit_ieta", &pfClusterHit_ieta);
	EventTree->SetBranchAddress("pfClusterHit_iphi", &pfClusterHit_iphi);
	EventTree->SetBranchAddress("pfClusterHit_iz", &pfClusterHit_iz);

	EventTree->SetBranchAddress("superCluster_energy", &superCluster_energy);
	EventTree->SetBranchAddress("superCluster_eta", &superCluster_eta);
	EventTree->SetBranchAddress("superCluster_phi", &superCluster_phi);
	EventTree->SetBranchAddress("superCluster_R", &superCluster_R);
	EventTree->SetBranchAddress("superCluster_etaWidth", &superCluster_etaWidth);
	EventTree->SetBranchAddress("superCluster_phiWidth", &superCluster_phiWidth);
	EventTree->SetBranchAddress("superCluster_ieta", &superCluster_ieta);
	EventTree->SetBranchAddress("superCluster_iphi", &superCluster_iphi);
	EventTree->SetBranchAddress("superCluster_iz", &superCluster_iz);

	EventTree->SetBranchAddress("superCluster_seedIndex", &superCluster_seedIndex);
	EventTree->SetBranchAddress("superCluster_pfClustersIndex", &superCluster_pfClustersIndex);

	EventTree->SetBranchAddress("psCluster_energy", &psCluster_energy);
	EventTree->SetBranchAddress("psCluster_eta", &psCluster_eta);
	EventTree->SetBranchAddress("psCluster_phi", &psCluster_phi);

	EventTree->SetBranchAddress("superCluster_dR_genScore_MatchedIndex", &superCluster_dR_genScore_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_dR_simScore_MatchedIndex", &superCluster_dR_simScore_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_n_shared_xtals_MatchedIndex", &superCluster_n_shared_xtals_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_sim_fraction_MatchedIndex", &superCluster_sim_fraction_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_sim_fraction_min1_MatchedIndex", &superCluster_sim_fraction_min1_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_sim_fraction_min3_MatchedIndex", &superCluster_sim_fraction_min3_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_sim_rechit_diff_MatchedIndex", &superCluster_sim_rechit_diff_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_sim_rechit_fraction_MatchedIndex", &superCluster_sim_rechit_fraction_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_global_sim_rechit_fraction_MatchedIndex", &superCluster_global_sim_rechit_fraction_MatchedIndex);
    EventTree->SetBranchAddress("superCluster_hgcal_caloToCluster_MatchedIndex", &superCluster_hgcal_caloToCluster_MatchedIndex);
	EventTree->SetBranchAddress("superCluster_hgcal_clusterToCalo_MatchedIndex", &superCluster_hgcal_clusterToCalo_MatchedIndex);

	EventTree->SetBranchAddress("superCluster_dR_genScore", &superCluster_dR_genScore);
	EventTree->SetBranchAddress("superCluster_dR_simScore", &superCluster_dR_simScore);
	EventTree->SetBranchAddress("superCluster_n_shared_xtals", &superCluster_n_shared_xtals);
	EventTree->SetBranchAddress("superCluster_sim_fraction", &superCluster_sim_fraction);
	EventTree->SetBranchAddress("superCluster_sim_fraction_min1", &superCluster_sim_fraction_min1);
	EventTree->SetBranchAddress("superCluster_sim_fraction_min3", &superCluster_sim_fraction_min3);
	EventTree->SetBranchAddress("superCluster_sim_rechit_diff", &superCluster_sim_rechit_diff);
	EventTree->SetBranchAddress("superCluster_sim_rechit_fraction", &superCluster_sim_rechit_fraction);
	EventTree->SetBranchAddress("superCluster_global_sim_rechit_fraction", &superCluster_global_sim_rechit_fraction);
    EventTree->SetBranchAddress("superCluster_hgcal_caloToCluster", &superCluster_hgcal_caloToCluster);
	EventTree->SetBranchAddress("superCluster_hgcal_clusterToCalo", &superCluster_hgcal_clusterToCalo);

    EvMax=EventTree->GetEntries();

    cout<<cur_time()<<"\tTree successfully processed!\n";
    cout<<cur_time()<<"\tCreating cluster vectors... \n";

   //outfile<<cur_time()<<"\tTree successfully processed!\n";
   //outfile<<cur_time()<<"\tCreating cluster vectors... \n";


    return;
}

void EventLoop(){
    int eventLoopMax = EvMax;

    bool separation_opt = true, 
         dphi_opt = false;

    //loop over events
    for(int iev=0; iev<eventLoopMax; ++iev){
        EventTree->GetEvent(iev);
        if(iev % 1000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev<<endl;

        //loop over caloParticles -  for total clusters in caloparticles plots
        for(int idxCalo=0; idxCalo<caloParticle_id->size(); idxCalo++){
            //loop over pfClusters
            for(int pf_idx = 0; pf_idx < pfCluster_energy->size(); pf_idx++){

                float pf_ET = pfCluster_energy->at(pf_idx) / cosh(pfCluster_eta->at(pf_idx));
                
                //if pfCluster is matched to caloParticle, plot in efficiency plot
                if(pfCluster_sim_fraction_min1_MatchedIndex->at(pf_idx) == idxCalo) effInt_totalCalo->Fill(pfCluster_eta->at(pf_idx), pf_ET);
            } //loop over pfClusters
        }//loop over caloParticles

        std::vector<int> seedCluster_etaBins(superCluster_seedIndex->size());
        //loop over superclusters - for original parameters
        for(int idx=0; idx<superCluster_seedIndex->size(); idx++){
            int currSeed = superCluster_seedIndex->at(idx);  //seed index
            float seedEta = pfCluster_eta->at(currSeed);
            int dPhi_bin = -1;
            float seed_ET = pfCluster_energy->at(currSeed) / cosh(pfCluster_eta->at(currSeed));

            if     (fabs(seedEta) >= 0.0   && fabs(seedEta) < 0.2)   { dPhi_bin = 0; }
            else if(fabs(seedEta) >= 0.2   && fabs(seedEta) < 0.4)   { dPhi_bin = 1; }
            else if(fabs(seedEta) >= 0.4   && fabs(seedEta) < 0.6)   { dPhi_bin = 2; } 
            else if(fabs(seedEta) >= 0.6   && fabs(seedEta) < 0.8)   { dPhi_bin = 3; }
            else if(fabs(seedEta) >= 0.8   && fabs(seedEta) < 1.0)   { dPhi_bin = 4; }
            else if(fabs(seedEta) >= 1.0   && fabs(seedEta) < 1.2)   { dPhi_bin = 5; } 
            else if(fabs(seedEta) >= 1.2   && fabs(seedEta) < 1.4)   { dPhi_bin = 6; }
            else if(fabs(seedEta) >= 1.4   && fabs(seedEta) < 1.479) { dPhi_bin = 7; }
            else if(fabs(seedEta) >= 1.479 && fabs(seedEta) < 1.75)  { dPhi_bin = 8; }
            else if(fabs(seedEta) >= 1.75  && fabs(seedEta) < 2.0)   { dPhi_bin = 9; }
            else if(fabs(seedEta) >= 2.0   && fabs(seedEta) < 2.25)  { dPhi_bin = 10; } 
            else if(fabs(seedEta) >= 2.25)                           { dPhi_bin = 11; } 
            
            seedCluster_etaBins[idx] = dPhi_bin;

            int seed_calo_idx = pfCluster_sim_fraction_min1_MatchedIndex->at(currSeed);
            if(seed_calo_idx < 0) continue;
            if(seed_calo_idx != superCluster_sim_fraction_min1_MatchedIndex->at(idx)) continue;

            //fill efficiency plots with seed info
            effInt_totalSC_original  -> Fill(seedEta, seed_ET);
            effInt_caloInSC_original -> Fill(seedEta, seed_ET);

            //loop over all pfClusters matched to the SuperCluster
            for(int pf_in_sc_idx = 0; pf_in_sc_idx < superCluster_pfClustersIndex->at(idx).size(); pf_in_sc_idx++){

                int pf_idx = superCluster_pfClustersIndex->at(idx).at(pf_in_sc_idx);
                if(pf_idx == currSeed) continue;  //skip the seeds

                int pf_calo_idx = pfCluster_sim_fraction_min1_MatchedIndex->at(pf_idx);

                //dEta and dPhi calculations
                float pf_dEta = (1-2*(pfCluster_eta->at(currSeed) < 0)) * (pfCluster_eta->at(pf_idx) - pfCluster_eta->at(currSeed));
                float pf_dPhi = pfCluster_phi->at(pf_idx) - pfCluster_phi->at(currSeed);
                if(pf_dPhi > pi) pf_dPhi -=twopi;
                if(pf_dPhi < -pi) pf_dPhi += twopi;
                float pf_ET = pfCluster_energy->at(pf_idx) / cosh(pfCluster_eta->at(pf_idx));

                //fill efficiency plots
                effInt_totalSC_original->Fill(pfCluster_eta->at(pf_idx), pf_ET);
                if(seed_calo_idx == pf_calo_idx) effInt_caloInSC_original->Fill(pfCluster_eta->at(pf_idx), pf_ET);
            } //loop over clusters in SC
        }
        
        vector<int> gatheredClusters;
        vector<int> pfClusters_ordered_ET;

        //order the pf Clusters in ET
        for(int pf = 0; pf < pfCluster_energy->size(); pf++){
            float curr_ET = pfCluster_energy->at(pf) / cosh(pfCluster_eta->at(pf));
            if(pfClusters_ordered_ET.size() == 0){ pfClusters_ordered_ET.push_back(pf); continue; } //add first element
            
            for(std::vector<int>::iterator it = pfClusters_ordered_ET.begin(); it != pfClusters_ordered_ET.end(); ++it){
                float lead_ET = pfCluster_energy->at(*it) / cosh(pfCluster_eta->at(*it));
                
                if(curr_ET > lead_ET){ pfClusters_ordered_ET.insert(it, pf); break; }
                if(it+1 == pfClusters_ordered_ET.end()){ pfClusters_ordered_ET.insert(it+1, pf); break; } // looking at the last element

                float trail_ET = pfCluster_energy->at(*(it+1)) / cosh(pfCluster_eta->at(*(it+1)));
                if(curr_ET < lead_ET && curr_ET > trail_ET){ pfClusters_ordered_ET.insert(it+1, pf); break;}

                /*
                if(pfClusters_ordered_ET.size() == 1){
                    if(curr_ET > lead_ET) pfClusters_ordered_ET.insert(it, pf);
                    else pfClusters_ordered_ET.insert(it+1, pf);
                    break;
                }
                if(it+1 == pfClusters_ordered_ET.end()){
                    if(curr_ET > lead_ET) pfClusters_ordered_ET.insert(it, pf);
                    else pfClusters_ordered_ET.insert(it+1, pf);
                    break;
                } 

                float trail_ET = pfCluster_energy->at(*(it+1)) / cosh(pfCluster_eta->at(*(it+1)));
                
                if(curr_ET > lead_ET) {pfClusters_ordered_ET.insert(it, pf); break;}
                else if(curr_ET < lead_ET && curr_ET > trail_ET){ pfClusters_ordered_ET.insert(it+1, pf); break;}
                else {continue;} */
            }
        }

        //loop over all clusters and find those that pass the seeding requirements
        for(int pf_order_idx = 0; pf_order_idx < pfCluster_energy->size(); pf_order_idx++){
            int pf_seed_idx = pfClusters_ordered_ET.at(pf_order_idx);

            //skip if already gathered
            if(find(gatheredClusters.begin(), gatheredClusters.end(), pf_seed_idx) != gatheredClusters.end()) continue;
            int seed_calo_idx = pfCluster_sim_fraction_min1_MatchedIndex->at(pf_seed_idx);
            if(seed_calo_idx < 0) continue;
            float seed_ET = pfCluster_energy->at(pf_seed_idx) / cosh(pfCluster_eta->at(pf_seed_idx));
            if(seed_ET < 1) continue;   //seeding et threshold

            //seed has passed all requirements:
            gatheredClusters.push_back(pf_seed_idx);
            float seedEta = pfCluster_eta->at(pf_seed_idx);
            int dPhi_bin = -1;

            if     (fabs(seedEta) >= 0.0   && fabs(seedEta) < 0.2)   { dPhi_bin = 0; }
            else if(fabs(seedEta) >= 0.2   && fabs(seedEta) < 0.4)   { dPhi_bin = 1; }
            else if(fabs(seedEta) >= 0.4   && fabs(seedEta) < 0.6)   { dPhi_bin = 2; } 
            else if(fabs(seedEta) >= 0.6   && fabs(seedEta) < 0.8)   { dPhi_bin = 3; }
            else if(fabs(seedEta) >= 0.8   && fabs(seedEta) < 1.0)   { dPhi_bin = 4; }
            else if(fabs(seedEta) >= 1.0   && fabs(seedEta) < 1.2)   { dPhi_bin = 5; } 
            else if(fabs(seedEta) >= 1.2   && fabs(seedEta) < 1.4)   { dPhi_bin = 6; }
            else if(fabs(seedEta) >= 1.4   && fabs(seedEta) < 1.479) { dPhi_bin = 7; }
            else if(fabs(seedEta) >= 1.479 && fabs(seedEta) < 1.75)  { dPhi_bin = 8; }
            else if(fabs(seedEta) >= 1.75  && fabs(seedEta) < 2.0)   { dPhi_bin = 9; }
            else if(fabs(seedEta) >= 2.0   && fabs(seedEta) < 2.25)  { dPhi_bin = 10; } 
            else if(fabs(seedEta) >= 2.25)                           { dPhi_bin = 11; } 
            
            //fill efficiency plots with seed info
            effInt_totalSC_reformed  -> Fill(seedEta, seed_ET);
            effInt_caloInSC_reformed -> Fill(seedEta, seed_ET);

            //loop over all clusters to gather them to the seeds
            for(int pf_idx = 0; pf_idx < pfCluster_energy->size(); pf_idx++){
                //skip if already gathered
                if(find(gatheredClusters.begin(), gatheredClusters.end(), pf_idx) != gatheredClusters.end()) continue;

                seed_calo_idx = pfCluster_sim_fraction_min1_MatchedIndex->at(pf_seed_idx);
                int pf_calo_idx = pfCluster_sim_fraction_min1_MatchedIndex->at(pf_idx);

                //dEta and dPhi calculations
                float pf_dEta = (1-2*(pfCluster_eta->at(pf_seed_idx) < 0)) * (pfCluster_eta->at(pf_idx) - pfCluster_eta->at(pf_seed_idx));
                float pf_dPhi = pfCluster_phi->at(pf_idx) - pfCluster_phi->at(pf_seed_idx);
                if(pf_dPhi > pi)  pf_dPhi -= twopi;
                if(pf_dPhi < -pi) pf_dPhi += twopi;
                float pf_ET = pfCluster_energy->at(pf_idx) / cosh(pfCluster_eta->at(pf_idx));

                //if the cluster doesn't fall in the dynamic dPhi window, skip
                if(dphi_opt) { if(!DynamicDPhi_opt (pf_ET, fabs(seedEta), pf_dPhi)) continue; }
                else         { if(!DynamicDPhi_orig(pf_ET, fabs(seedEta), pf_dPhi)) continue; }

                int et_bin = -1;

                if     (pf_ET <  0.5)              et_bin = 1;
                else if(pf_ET >= 0.5 && pf_ET < 1) et_bin = 2;
                else if(pf_ET >= 1   && pf_ET < 2) et_bin = 3;
                else if(pf_ET >= 2   && pf_ET < 3) et_bin = 4;
                else if(pf_ET >= 3   && pf_ET < 6) et_bin = 5;
                else if(pf_ET >= 6)                et_bin = 6;

                if(et_bin == -1 || dPhi_bin == -1) continue; //error checking

                bool par_pass = false;

                if(separation_opt)
                    par_pass =  upper_parabola(seedEta, pf_ET, pf_dEta, pf_dPhi, w00_opt, w01_opt, w10_opt, w11_opt) 
                             && lower_parabola(seedEta, pf_ET, pf_dEta, pf_dPhi, w00_opt, w01_opt, w10_opt, w11_opt);

                else
                    par_pass =  upper_parabola(seedEta, pf_ET, pf_dEta, pf_dPhi, w00_orig, w01_orig, w10_orig, w11_orig) 
                             && lower_parabola(seedEta, pf_ET, pf_dEta, pf_dPhi, w00_orig, w01_orig, w10_orig, w11_orig);

                if(par_pass){
                    gatheredClusters.push_back(pf_idx);
                    reformed_mustaches[et_bin][dPhi_bin]->Fill(pf_dPhi, pf_dEta);
                    reformed_mustaches[0][dPhi_bin]->Fill(pf_dPhi, pf_dEta);

                    effInt_totalSC_reformed -> Fill(pfCluster_eta->at(pf_idx), pf_ET);
                    if(seed_calo_idx == pf_calo_idx) effInt_caloInSC_reformed -> Fill(pfCluster_eta->at(pf_idx), pf_ET);

                    if(pf_ET > sc_reformed_maxET[et_bin][dPhi_bin]){ sc_reformed_maxET[et_bin][dPhi_bin] = pf_ET; sc_reformed_maxET_seedEta[et_bin][dPhi_bin] = seedEta; }
                    if(pf_ET < sc_reformed_minET[et_bin][dPhi_bin]){ sc_reformed_minET[et_bin][dPhi_bin] = pf_ET; sc_reformed_minET_seedEta[et_bin][dPhi_bin] = seedEta; }
                }

            }//loop over clusters to gather them to the seeds
        }//loop over clusters and see the seeds

    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;
}

Double_t upper_parabola_draw(Double_t *x, Double_t *par){
//Draw the upper parabola for given seed_eta, pfCluster ET, and w parameters
    float eta_max = par[0],//1.45;
          max_E_T = par[1],//30.4012;
          w00     = par[2],//w00_opt;
          w01     = par[3],//w01_opt;
          w10     = par[4],//w10_opt;
          w11     = par[5];//w11_opt;

	float p00 = -0.107537;
	float p01 = 0.590969;
	float p02 = -0.076494;
	float p10 = -0.0268843;
	float p11 = 0.147742;
	float p12 = -0.0191235;
    

    float c_upper = (p00 * pow(eta_max*sin(eta_max), 2)) + (p01 * eta_max*sin(eta_max)) + p02;
    float d_upper =( w10*eta_max*sin(eta_max) ) + (w11 / sqrt(1.1 + log10(max_E_T)));
    float d_lower =( w00*eta_max*sin(eta_max) ) + (w01 / sqrt(1.1 + log10(max_E_T)));
    float b_upper = d_upper - 0.5*(d_lower + d_upper);
    float a_upper = ((1 / (4 * c_upper))) - fabs(b_upper);
    Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f)*pow(x[0],2)) + std::max(b_upper, 0.0087f) + 0.0087;

    return upper_curve;
}

Double_t lower_parabola_draw(Double_t *x, Double_t *par){
//Draw the lower parabola for given seed_eta, pfCluster ET, and w parameters
    float eta_max = par[0],//1.45;
          min_E_T = par[1],//0.123621;
          w00     = par[2],//w00_opt;
          w01     = par[3],//w01_opt;
          w10     = par[4],//w10_opt;
          w11     = par[5];//w11_opt;

    float p00 = -0.107537;
    float p01 = 0.590969;
    float p02 = -0.076494;
    float p10 = -0.0268843;
    float p11 = 0.147742;
    float p12 = -0.0191235;


    float c_lower = (p10 * pow(eta_max*sin(eta_max), 2)) + (p11 * eta_max*sin(eta_max)) + p12;
    float d_upper =( w10*eta_max*sin(eta_max) ) + (w11 / sqrt(1.1 + log10(min_E_T)));
    float d_lower =( w00*eta_max*sin(eta_max) ) + (w01 / sqrt(1.1 + log10(min_E_T)));
    float b_lower = d_lower - 0.5*(d_lower + d_upper);
    float a_lower = ((1 / (4 * c_lower))) - fabs(b_lower);
    Double_t lower_curve = (std::max((1 / (4 * a_lower)),0.0f)*pow(x[0],2)) + std::min(b_lower, -0.0087f);

    return lower_curve;
}

void draw_sc(){
//Draw the reformed mustache superclusters
    TFile *recalc_SC= new TFile("recalc_SC.root","RECREATE");
    c1->Clear(); c1->Divide(4,3);
    //ET BINS
    TLegend* mustaches_dphi_etbins_legends[7][12]; 

    TF1* upper_parabolas_dPhi_ETBins_max[7][12];
    TF1* lower_parabolas_dPhi_ETBins_max[7][12];
    TF1* upper_parabolas_dPhi_ETBins_min[7][12];
    TF1* lower_parabolas_dPhi_ETBins_min[7][12];

    TF1* upper_parabolas_dPhi_ETBins_max_original[7][12];
    TF1* lower_parabolas_dPhi_ETBins_max_original[7][12];
    TF1* upper_parabolas_dPhi_ETBins_min_original[7][12];
    TF1* lower_parabolas_dPhi_ETBins_min_original[7][12];

    TLine* max_dphi_lines_etBins_original[7][12];
    TLine* min_dphi_lines_etBins_original[7][12];
    TLine* max_dphi_lines_etBins_left_original[7][12];
    TLine* min_dphi_lines_etBins_left_original[7][12];

    TLine* max_dphi_lines_etBins_opt[7][12];
    TLine* min_dphi_lines_etBins_opt[7][12];
    TLine* max_dphi_lines_etBins_left_opt[7][12];
    TLine* min_dphi_lines_etBins_left_opt[7][12];

    float mid_etas[12] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.4395,1.6145,1.875,2.125,2.625};

    //loop over et Bins
    //for(int x=1; x<2;x++){//
    for(int x=0; x<7;x++){
        c1->cd();
        c1->Clear();
        c1->Divide(4,3);

      //loop over eta Bins
        //for(int y=3; y<4;y++){//
        for(int y=0; y<12;y++){

            //Plotting
            c1->cd(y+1);
            gStyle->SetOptStat(0);
            //gStyle->SetStatX(0.87);
            //gStyle->SetStatY(0.37);
            reformed_mustaches[x][y]->SetTitle((et_bin_titles[x] + "    " + titles_etas_dPhi[y]).c_str());
            reformed_mustaches[x][y]->SetMarkerColor(2); //in red
            reformed_mustaches[x][y]->Draw("COLZ"); c1->Update();

            float DPhi_MAX_original = DynamicDPhi_orig_val(sc_reformed_maxET[x][y],fabs(sc_reformed_maxET_seedEta[x][y])),
                  DPhi_MIN_original = DynamicDPhi_orig_val(sc_reformed_minET[x][y],fabs(sc_reformed_minET_seedEta[x][y])),
                  DPhi_MAX_opt = DynamicDPhi_opt_val(sc_reformed_maxET[x][y],fabs(sc_reformed_maxET_seedEta[x][y])),
                  DPhi_MIN_opt = DynamicDPhi_opt_val(sc_reformed_minET[x][y],fabs(sc_reformed_minET_seedEta[x][y]));

            if(x != 0){ // draw the parabolas 
                upper_parabolas_dPhi_ETBins_max[x][y] = new TF1(("upper_parabolas_dphi_etBin_"+to_string(x)+"etaBin_"+to_string(y)+"_max").c_str(), upper_parabola_draw, -0.6, 0.6,6);
                lower_parabolas_dPhi_ETBins_max[x][y] = new TF1(("lower_parabolas_dphi_etBin_"+to_string(x)+"etaBin_"+to_string(y)+"_max").c_str(), lower_parabola_draw, -0.6, 0.6,6);
                upper_parabolas_dPhi_ETBins_max[x][y]->SetParameters(sc_reformed_maxET_seedEta[x][y],sc_reformed_maxET[x][y],w00_opt,w01_opt,w10_opt,w11_opt);
                lower_parabolas_dPhi_ETBins_max[x][y]->SetParameters(sc_reformed_maxET_seedEta[x][y],sc_reformed_maxET[x][y],w00_opt,w01_opt,w10_opt,w11_opt);

                upper_parabolas_dPhi_ETBins_min[x][y] = new TF1(("upper_parabolas_dphi_etBin_"+to_string(x)+"etaBin_"+to_string(y)+"_min").c_str(), upper_parabola_draw, -0.6, 0.6,6);
                lower_parabolas_dPhi_ETBins_min[x][y] = new TF1(("lower_parabolas_dphi_etBin_"+to_string(x)+"etaBin_"+to_string(y)+"_min").c_str(), lower_parabola_draw, -0.6, 0.6,6);
                upper_parabolas_dPhi_ETBins_min[x][y]->SetParameters(sc_reformed_minET_seedEta[x][y],sc_reformed_minET[x][y],w00_opt,w01_opt,w10_opt,w11_opt);
                lower_parabolas_dPhi_ETBins_min[x][y]->SetParameters(sc_reformed_minET_seedEta[x][y],sc_reformed_minET[x][y],w00_opt,w01_opt,w10_opt,w11_opt);
                
                upper_parabolas_dPhi_ETBins_max_original[x][y] = new TF1(("upper_parabolas_dphi_etBin_original_"+to_string(x)+"etaBin_"+to_string(y)+"_max").c_str(), upper_parabola_draw, -0.6, 0.6,6);
                lower_parabolas_dPhi_ETBins_max_original[x][y] = new TF1(("lower_parabolas_dphi_etBin_original_"+to_string(x)+"etaBin_"+to_string(y)+"_max").c_str(), lower_parabola_draw, -0.6, 0.6,6);
                upper_parabolas_dPhi_ETBins_max_original[x][y]->SetParameters(sc_reformed_maxET_seedEta[x][y],sc_reformed_maxET[x][y],w00_orig,w01_orig,w10_orig,w11_orig);
                lower_parabolas_dPhi_ETBins_max_original[x][y]->SetParameters(sc_reformed_maxET_seedEta[x][y],sc_reformed_maxET[x][y],w00_orig,w01_orig,w10_orig,w11_orig);

                upper_parabolas_dPhi_ETBins_min_original[x][y] = new TF1(("upper_parabolas_dphi_etBin_original_"+to_string(x)+"etaBin_"+to_string(y)+"_min").c_str(), upper_parabola_draw, -0.6, 0.6,6);
                lower_parabolas_dPhi_ETBins_min_original[x][y] = new TF1(("lower_parabolas_dphi_etBin_original_"+to_string(x)+"etaBin_"+to_string(y)+"_min").c_str(), lower_parabola_draw, -0.6, 0.6,6);
                upper_parabolas_dPhi_ETBins_min_original[x][y]->SetParameters(sc_reformed_minET_seedEta[x][y],sc_reformed_minET[x][y],w00_orig,w01_orig,w10_orig,w11_orig);
                lower_parabolas_dPhi_ETBins_min_original[x][y]->SetParameters(sc_reformed_minET_seedEta[x][y],sc_reformed_minET[x][y],w00_orig,w01_orig,w10_orig,w11_orig);

                max_dphi_lines_etBins_original[x][y] = new TLine(DPhi_MAX_original, -0.15, DPhi_MAX_original, 0.2);
                min_dphi_lines_etBins_original[x][y] = new TLine(DPhi_MIN_original, -0.15, DPhi_MIN_original, 0.2);
                max_dphi_lines_etBins_left_original[x][y] = new TLine(-(DPhi_MAX_original), -0.15, -(DPhi_MAX_original), 0.2);
                min_dphi_lines_etBins_left_original[x][y] = new TLine(-(DPhi_MIN_original), -0.15, -(DPhi_MIN_original), 0.2);

                max_dphi_lines_etBins_opt[x][y] = new TLine(DPhi_MAX_opt, -0.15, DPhi_MAX_opt, 0.2);
                min_dphi_lines_etBins_opt[x][y] = new TLine(DPhi_MIN_opt, -0.15, DPhi_MIN_opt, 0.2);
                max_dphi_lines_etBins_left_opt[x][y] = new TLine(-(DPhi_MAX_opt), -0.15, -(DPhi_MAX_opt), 0.2);
                min_dphi_lines_etBins_left_opt[x][y] = new TLine(-(DPhi_MIN_opt), -0.15, -(DPhi_MIN_opt), 0.2);

                
                upper_parabolas_dPhi_ETBins_max[x][y]->SetLineColor(5); //yellow
                lower_parabolas_dPhi_ETBins_max[x][y]->SetLineColor(5); //yellow
                upper_parabolas_dPhi_ETBins_min[x][y]->SetLineColor(2); //red
                lower_parabolas_dPhi_ETBins_min[x][y]->SetLineColor(2); //red

                upper_parabolas_dPhi_ETBins_max[x][y]->SetLineStyle(9); c1->Update();
                lower_parabolas_dPhi_ETBins_max[x][y]->SetLineStyle(9); c1->Update();
                upper_parabolas_dPhi_ETBins_min[x][y]->SetLineStyle(9); c1->Update();
                lower_parabolas_dPhi_ETBins_min[x][y]->SetLineStyle(9); c1->Update();
                upper_parabolas_dPhi_ETBins_max[x][y]->SetLineWidth(7);
                lower_parabolas_dPhi_ETBins_max[x][y]->SetLineWidth(7);
                upper_parabolas_dPhi_ETBins_min[x][y]->SetLineWidth(7);
                lower_parabolas_dPhi_ETBins_min[x][y]->SetLineWidth(7);
                upper_parabolas_dPhi_ETBins_max[x][y]->SetLineStyle(9); c1->Update();
                lower_parabolas_dPhi_ETBins_max[x][y]->SetLineStyle(9); c1->Update();
                upper_parabolas_dPhi_ETBins_min[x][y]->SetLineStyle(9); c1->Update();
                lower_parabolas_dPhi_ETBins_min[x][y]->SetLineStyle(9); c1->Update();

                //upper_parabolas_dPhi_ETBins_max[x][y]->Draw("SAME");
                //lower_parabolas_dPhi_ETBins_max[x][y]->Draw("SAME");
                //upper_parabolas_dPhi_ETBins_min[x][y]->Draw("SAME");
                //lower_parabolas_dPhi_ETBins_min[x][y]->Draw("SAME");
                

                upper_parabolas_dPhi_ETBins_max_original[x][y]->SetLineColor(6); 
                lower_parabolas_dPhi_ETBins_max_original[x][y]->SetLineColor(6);
                upper_parabolas_dPhi_ETBins_min_original[x][y]->SetLineColor(3);
                lower_parabolas_dPhi_ETBins_min_original[x][y]->SetLineColor(3);

                upper_parabolas_dPhi_ETBins_max_original[x][y]->SetLineStyle(10);
                lower_parabolas_dPhi_ETBins_max_original[x][y]->SetLineStyle(10);
                upper_parabolas_dPhi_ETBins_min_original[x][y]->SetLineStyle(10);
                lower_parabolas_dPhi_ETBins_min_original[x][y]->SetLineStyle(10);

                upper_parabolas_dPhi_ETBins_max_original[x][y]->SetLineWidth(5);
                lower_parabolas_dPhi_ETBins_max_original[x][y]->SetLineWidth(5);
                upper_parabolas_dPhi_ETBins_min_original[x][y]->SetLineWidth(5);
                lower_parabolas_dPhi_ETBins_min_original[x][y]->SetLineWidth(5);

                upper_parabolas_dPhi_ETBins_max_original[x][y]->SetLineStyle(10); c1->Update();
                lower_parabolas_dPhi_ETBins_max_original[x][y]->SetLineStyle(10); c1->Update();
                upper_parabolas_dPhi_ETBins_min_original[x][y]->SetLineStyle(10); c1->Update();
                lower_parabolas_dPhi_ETBins_min_original[x][y]->SetLineStyle(10); c1->Update();

                upper_parabolas_dPhi_ETBins_max_original[x][y]->Draw("SAME");
                lower_parabolas_dPhi_ETBins_max_original[x][y]->Draw("SAME");
                upper_parabolas_dPhi_ETBins_min_original[x][y]->Draw("SAME");
                lower_parabolas_dPhi_ETBins_min_original[x][y]->Draw("SAME");



                max_dphi_lines_etBins_opt[x][y]->SetLineColor(5);
                min_dphi_lines_etBins_opt[x][y]->SetLineColor(2);
                max_dphi_lines_etBins_left_opt[x][y]->SetLineColor(5);
                min_dphi_lines_etBins_left_opt[x][y]->SetLineColor(2);

                max_dphi_lines_etBins_opt[x][y]->SetLineWidth(7);
                min_dphi_lines_etBins_opt[x][y]->SetLineWidth(7);
                max_dphi_lines_etBins_left_opt[x][y]->SetLineWidth(7);
                min_dphi_lines_etBins_left_opt[x][y]->SetLineWidth(7);

                max_dphi_lines_etBins_opt[x][y]->Draw("SAME");
                min_dphi_lines_etBins_opt[x][y]->Draw("SAME");
                max_dphi_lines_etBins_left_opt[x][y]->Draw("SAME");
                min_dphi_lines_etBins_left_opt[x][y]->Draw("SAME");


                max_dphi_lines_etBins_original[x][y]->SetLineColor(6);
                min_dphi_lines_etBins_original[x][y]->SetLineColor(3);
                max_dphi_lines_etBins_left_original[x][y]->SetLineColor(6);
                min_dphi_lines_etBins_left_original[x][y]->SetLineColor(3);

                max_dphi_lines_etBins_original[x][y]->SetLineWidth(5);
                min_dphi_lines_etBins_original[x][y]->SetLineWidth(5);
                max_dphi_lines_etBins_left_original[x][y]->SetLineWidth(5);
                min_dphi_lines_etBins_left_original[x][y]->SetLineWidth(5);

                max_dphi_lines_etBins_original[x][y]->SetLineStyle(10);
                min_dphi_lines_etBins_original[x][y]->SetLineStyle(10);
                max_dphi_lines_etBins_left_original[x][y]->SetLineStyle(10);
                min_dphi_lines_etBins_left_original[x][y]->SetLineStyle(10);

                max_dphi_lines_etBins_original[x][y]->Draw("SAME");
                min_dphi_lines_etBins_original[x][y]->Draw("SAME");
                max_dphi_lines_etBins_left_original[x][y]->Draw("SAME");
                min_dphi_lines_etBins_left_original[x][y]->Draw("SAME");

                

                mustaches_dphi_etbins_legends[x][y] = new TLegend(0.1,-0.09,0.5,-0.04,"","");
                //mustaches_dphi_etbins_legends[x][y]->AddEntry(upper_parabolas_dPhi_ETBins_max[x][y], "Maximum ET Cluster - REFORMED","l");
                //mustaches_dphi_etbins_legends[x][y]->AddEntry(upper_parabolas_dPhi_ETBins_min[x][y], "Minimum ET Cluster - REFORMED","l");
                mustaches_dphi_etbins_legends[x][y]->AddEntry(max_dphi_lines_etBins_opt[x][y], "Maximum ET Cluster - REFORMED","l");
                mustaches_dphi_etbins_legends[x][y]->AddEntry(min_dphi_lines_etBins_opt[x][y], "Minimum ET Cluster - REFORMED","l");
                mustaches_dphi_etbins_legends[x][y]->AddEntry(upper_parabolas_dPhi_ETBins_max_original[x][y], "Maximum ET Cluster - ORIGINAL","l");
                mustaches_dphi_etbins_legends[x][y]->AddEntry(upper_parabolas_dPhi_ETBins_min_original[x][y], "Minimum ET Cluster - ORIGINAL","l");
                mustaches_dphi_etbins_legends[x][y]->Draw("SAME");
            }
        }

        c1->cd(); c1->Modified(); c1->Update();
        c1->SaveAs(("reform_SC/mustache_etbin_"+to_string(x)+".png").c_str());
    }
    recalc_SC->Write();
    recalc_SC->Close();


    return;
}

void plotCaloEfficiency_Integral(){
    c1 -> Clear();

    TH1F* effInt_caloInSC_Clone = (TH1F*) effInt_caloInSC_reformed -> Clone();
    //TH1F* effInt_totalSC_Clone = (TH1F*) effInt_totalSC->Clone();

    effInt_caloInSC_reformed->Draw();
    c1->SaveAs("reform_SC/calo_in_sc_ref.png");
    effInt_totalCalo->Draw();
    c1->SaveAs("reform_SC/totalcalo.png");
    effInt_caloInSC_original->Draw();
    c1->SaveAs("reform_SC/calo_in_sc_orig.png");


    //efficiency - reformed
    effInt_caloInSC_reformed -> Divide(effInt_totalCalo);
    
    //calculate error bars
    for(int i = 1; i < 61; i++){
        float eff   = effInt_caloInSC_reformed -> GetBinContent(i);
        float total = effInt_caloInSC_reformed -> GetEntries();
        float error = sqrt(eff*(1-eff)/total);  // binomial error bars

        effInt_caloInSC_reformed -> SetBinError(i, error);
    }

    //efficiency - original
    effInt_caloInSC_original -> Divide(effInt_totalCalo);
    
    //calculate error bars
    for(int i = 1; i < 61; i++){
        float eff   = effInt_caloInSC_original -> GetBinContent(i);
        float total = effInt_caloInSC_original -> GetEntries();
        float error = sqrt(eff*(1-eff)/total);  // binomial error bars

        effInt_caloInSC_original -> SetBinError(i, error);
    }

    effInt_caloInSC_reformed -> GetXaxis() -> SetTitle("PF Cluster #eta");
    effInt_caloInSC_reformed -> GetYaxis() -> SetRangeUser(0.98,1.0);
    effInt_caloInSC_reformed -> SetLineWidth(2);
    effInt_caloInSC_reformed -> SetLineColor(2);  // red
    effInt_caloInSC_reformed -> Draw();

    effInt_caloInSC_original -> GetXaxis() -> SetTitle("PF Cluster #eta");
    effInt_caloInSC_original -> GetYaxis() -> SetRangeUser(0.98,1.0);
    effInt_caloInSC_original -> SetLineWidth(2);
    effInt_caloInSC_original -> SetLineColor(4);  // blue
    effInt_caloInSC_original -> Draw("SAME");

    c1 -> SaveAs("reform_SC/SC_efficiency_integral.png");
    c1 -> Clear();
}

void getData(){
    cout<<cur_time()<<"\tAttempting to retrieve data"<<endl;
    ifstream fin;
    string dir = "/eos/cms/store/user/lzygala/par_con_files/output_files/", 
           filepath;
    DIR *dp;
    struct dirent *dirp;
    struct stat filestat;

    float read_w00, read_w01, read_w10, read_w11, read_eff;
    int files_read = 0;

    dp = opendir(dir.c_str());
    if(dp == NULL){
        cout << cur_time() << "\tError reading files" << endl;
        return;
    }

    cout << cur_time() << "\tReading data from: " << dir.c_str() << endl;

    while(dirp = readdir(dp)){
        filepath = dir + dirp->d_name;

        if(stat(filepath.c_str(), &filestat)) continue;
        if(S_ISDIR(filestat.st_mode)) continue;

        fin.open(filepath.c_str());
        while(!fin.eof()){
            if(fin >> read_w00 >> read_w01 >> read_w10 >> read_w11 >> read_eff){
                list_w00.push_back(read_w00);
                list_w01.push_back(read_w01);
                list_w10.push_back(read_w10);
                list_w11.push_back(read_w11);
                list_eff.push_back(read_eff);
            }
        }
        fin.close();
        if(++files_read % 100 == 0)
            cout << cur_time() << "\tFiles Processed: " << files_read << endl;
    }

    closedir(dp);
    cout << cur_time() << "\tAll data retrieved. Total files processed: " << files_read << endl;
    cout << cur_time() << "\tTotal number of entries: " << list_eff.size() << endl;
    return;
}

void findEff(){

    vector<float> highest_eff_idx;
        
    cout    << cur_time() << "\tFinding best efficiency...\n";

   if(list_eff.size() == 0){ cout << cur_time()<<"\tThere are no groups of values!"<<endl; return;}

    float best_efficiency = *std::max_element(list_eff.begin(), list_eff.end());
    cout    << cur_time() << "\tBest Efficiency: " << best_efficiency << endl;
    std::vector<float>::iterator it=std::find(list_eff.begin(), list_eff.end(), best_efficiency);
    highest_eff_idx.clear();
    int starting = 0;
    while(it != list_eff.end()){
        starting = std::distance(list_eff.begin(), it) + 1;
        highest_eff_idx.push_back(std::distance(list_eff.begin(), it));
        it = std::find(list_eff.begin() + starting,list_eff.end(), best_efficiency);
    }

    cout    << cur_time() << "\tNumber of groups with Best Efficiency: " << highest_eff_idx.size() << endl;
   //outfile << cur_time() << "\tNumber of groups with Best Efficiency: " << highest_eff_idx.size() << endl;

    cout    << cur_time() << "\tFinding best separation values...\n";
   //outfile << cur_time() << "\tFinding best separation values...\n";

  /*
    vector<float> sepAves;
    //Want to choose the parameters that give the smallest separation between the parabolas
    //Want values closest to minimum separation. Should i take the average?
    for(std::vector<float>::iterator iHighIDX = highest_eff_idx.begin(); iHighIDX != highest_eff_idx.end(); ++iHighIDX){
        sepAves.push_back(averageSeparationValue_event(list_w00.at(*iHighIDX), list_w01.at(*iHighIDX), list_w10.at(*iHighIDX), list_w11.at(*iHighIDX)));
    }

    std::vector<float>::iterator smallest_separation = std::min_element(sepAves.begin(),sepAves.end());
    cout    << cur_time() << "\tBest Separation: " << *smallest_separation << endl;
   //outfile << cur_time() << "\tBest Separation: " << *smallest_separation << endl;

    vector<float> min_sep_idx;
    min_sep_idx.clear();
    while(smallest_separation != sepAves.end()){
        starting = std::distance(sepAves.begin(), smallest_separation) + 1;
        min_sep_idx.push_back(highest_eff_idx.at(std::distance(sepAves.begin(), smallest_separation)));
        smallest_separation = std::find(sepAves.begin() + starting,sepAves.end(), *smallest_separation);
    }

    cout    << cur_time() << "\tNumber of groups with Best Separation: " << min_sep_idx.size() << endl;
   //outfile << cur_time() << "\tNumber of groups with Best Separation: " << min_sep_idx.size() << endl;

    //Print all parameter group options
    int optionNum = 1;
    for(std::vector<float>::iterator iMinSepIDX = min_sep_idx.begin(); iMinSepIDX != min_sep_idx.end(); ++iMinSepIDX){
        cout<<"Option "<<optionNum<<":\n";
        cout<<"\tEfficiency:\t"<<total_efficiency.at(*iMinSepIDX)
            <<"\n\tw00:\t"<<test_w00.at(*iMinSepIDX)
            <<"\n\tw01:\t"<<test_w01.at(*iMinSepIDX)
            <<"\n\tw10:\t"<<test_w10.at(*iMinSepIDX)
            <<"\n\tw11:\t"<<test_w11.at(*iMinSepIDX)
            <<"\n";

            optionNum++;

       //outfile<<"Option "<<optionNum<<":\n";
       //outfile<<"\tEfficiency:\t"<<total_efficiency.at(*iMinSepIDX)
            <<"\n\tw00:\t"<<test_w00.at(*iMinSepIDX)
            <<"\n\tw01:\t"<<test_w01.at(*iMinSepIDX)
            <<"\n\tw10:\t"<<test_w10.at(*iMinSepIDX)
            <<"\n\tw11:\t"<<test_w11.at(*iMinSepIDX)
            <<"\n";
    }
  /*
    if(min_sep_idx.size()>1){
        cout    << cur_time() << "\tMore than 1 group has the same separation!\n";
        //outfile << cur_time() << "\tMore than 1 group has the same separation!\n";
       // break;
    }

    best_w00 = test_w00.at(min_sep_idx.at(0));
    best_w10 = test_w10.at(min_sep_idx.at(0));
    best_w01 = test_w01.at(min_sep_idx.at(0));
    best_w11 = test_w11.at(min_sep_idx.at(0));
  */
}

void sc_recalculation(){
    InitTree("root://eoscms.cern.ch///eos/cms/store/user/lzygala/par_con_files/finalDumperTree_4G_Run3_PU.root");
    EventLoop();
    draw_sc();
    plotCaloEfficiency_Integral();
    //draw_sc();
    //getData();
    //findEff();
    //plotEfficiency_Integral();

    return;
}