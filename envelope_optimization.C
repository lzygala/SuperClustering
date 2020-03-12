
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

//Add efficiency per each eta bin

/* ---------- Variables ---------- */
    //ofstream outfile;
    ofstream data_outfile;

    TH2F* caloClusters_width_vs_logET[15];
    TH2F* caloClusters_integral;
    TH2F* caloClusters_shape[15][15];
    TH2F* cluster_dPhi_vs_loget[4];

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

    //vectors and variables
    vector<float> fit_w00_upper;
    vector<float> fit_w01_upper;
    vector<float> fit_w10_upper;
    vector<float> fit_w11_upper;

    vector<float> fit_w00_lower;
    vector<float> fit_w01_lower;
    vector<float> fit_w10_lower;
    vector<float> fit_w11_lower;

    float best_w00, best_w01, best_w10, best_w11;

    vector<float> lower_efficiency;
    vector<float> upper_efficiency;
    vector<float> total_efficiency;

    vector<float> highest_eff_idx;
    vector<float> min_sep_idx;

    //cluster vectors
    //They hold information on the clusters that pass the following requirements:
    //  1. Clusters that are in the caloCluster
    //  1. Cluster ET > 0
    //  2. Cluster is within dPhi window
    //  3. Cluster passes matching requirement
    float totalClusterCount;
    vector<float> cluster_ET;
    vector<float> cluster_Eta;
    vector<float> cluster_dPhi;
    vector<float> cluster_dEta;
    vector<float> cluster_scSeed_Eta;

    vector<float> integers;
    vector<float> upperParabolaValues;
    vector<float> lowerParabolaValues;

    //event vectors - for recalculating superclusters
    vector<vector<int>> eventSeeds; //saves indexes
    vector<vector<float>> eventCluster_ET;
    vector<vector<float>> eventCluster_Eta;
    vector<vector<float>> eventCluster_Phi;
    vector<vector<int>> eventCluster_CaloParticle;
    vector<vector<float>> eventCluster_CaloParticle_ET;
    vector<vector<float>> eventCluster_CaloParticle_Energy;

    //vectors for recalculating superclusters - only saving cluster info near seeds
    vector<float> gatheringClusters_caloET;
    vector<float> gatheringClusters_caloET_clusters;
    vector<float> gatheringClusters_seedEta;
    vector<float> gatheringClusters_seedPhi;
    vector<float> gatheringClusters_seedET;
    vector<vector<float>> gatheringClusters_clusterET;
    vector<vector<float>> gatheringClusters_clusterDEta;
    vector<vector<float>> gatheringClusters_clusterDPhi;


    //separation fitting vectors:
    vector<vector<vector<float>>> caloParticle_widths(15, vector<vector<float>>(30));

    //multithreading stuff
    mutex eff_return_val_mutex, upper_mutex, lower_mutex;

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

    for(int i=0;i<15;i++){
       caloClusters_width_vs_logET[i] = new TH2F(("caloClusters_width_vs_logET_etaBin_"+to_string(i)).c_str(),("caloClusters_width_vs_logET_etaBin_"+to_string(i)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
       for(int j=0; j<15;j++){
           caloClusters_shape[i][j] = new TH2F(("caloClusters_shape_"+to_string(i)+"_"+to_string(j)).c_str(),("caloClusters_shape_"+to_string(i)+"_"+to_string(j)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
       }
    }
    caloClusters_integral = new TH2F("caloClusters_integral","caloClusters_integral",bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
    
    for(int ii=0; ii<4; ii++){
        cluster_dPhi_vs_loget[ii] = new TH2F(("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),200, -1.0, 3.0, 200, 0.0, 1.0);
    }


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
//Loops over all events and fills histograms
    int eventLoopMax = EvMax;

    for(int iev=0; iev<eventLoopMax;++iev){
        EventTree->GetEvent(iev);
        if(iev % 1000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev<<endl;

        //cout<<"A\n";
        std::vector<int> seedCluster_etaBins(superCluster_seedIndex->size());
        std::vector<int> seedCluster_dPhiBins(superCluster_seedIndex->size());
        for(int idx=0; idx<superCluster_seedIndex->size(); idx++){
        //loop over superclusters - get seed cluster eta bin
            int currSeed = superCluster_seedIndex->at(idx);  //seed index
            float seedEta = pfCluster_eta->at(currSeed);
            int dPhi_bin = -1, seedEta_bin;

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
            
            seedCluster_etaBins[idx] = seedEta_bin;
            seedCluster_dPhiBins[idx] = dPhi_bin;

        }// loop over superclusters - get seed cluster eta bin

        
        //loop over calo_particles
        for(int idxCalo=0; idxCalo<caloParticle_id->size(); idxCalo++){
            vector<float> cluster_Eta;
            vector<vector<vector<float>>> cluster_dEta(15, vector<vector<float>>(30));
            
            //find associated supercluster
            std::vector<int>::iterator itSC = superCluster_sim_fraction_min1_MatchedIndex->begin();
            int idxSC = -1, counter = 0;
            while(itSC != superCluster_sim_fraction_min1_MatchedIndex -> end()){
                itSC = find(itSC, superCluster_sim_fraction_min1_MatchedIndex -> end(), idxCalo);
                if(itSC != superCluster_sim_fraction_min1_MatchedIndex -> end()){
                    idxSC = std::distance(superCluster_sim_fraction_min1_MatchedIndex -> begin(), itSC); 
                    counter++; itSC++;
                }
            }
            if(idxSC == -1) continue; //no SuperCluster associated
            if(counter > 1) continue; //more than 1 SuperCluster associated

            int currSeed = superCluster_seedIndex->at(idxSC);
            if(seedCluster_etaBins.at(idxSC) == -1) continue;  // skip the seeds

            int seed_eta_bin  = seedCluster_etaBins.at(idxSC),
                seed_dphi_bin = seedCluster_dPhiBins.at(idxSC);
            
            //Sanity check that caloParticle matched to seed and supercluster are the same????
            if(pfCluster_sim_fraction_min1_MatchedIndex->at(currSeed) != idxCalo){
                //cout<<"THE SEED AND SUPERCLUSTER ARE MATCHED TO DIFFERENT CALO"<<endl;
                continue;
            }


            //loop over the pf clusters
            for(int idxPF=0; idxPF<pfCluster_energy->size(); idxPF++){

                //skip the seeds
                if(currSeed == idxPF) continue;

                //pfCluster calculations
                float dEta = (1-2*(pfCluster_eta->at(currSeed) < 0)) * (pfCluster_eta->at(idxPF) - pfCluster_eta->at(currSeed));
                float dPhi = pfCluster_phi->at(idxPF) - pfCluster_phi->at(currSeed);
                if(dPhi > pi) dPhi -=twopi;
                if(dPhi < -pi) dPhi += twopi;
                float ET = pfCluster_energy->at(idxPF) / cosh(pfCluster_eta->at(idxPF));
                float logET = log10(ET);

                //pfCluster does not belong to caloParticle
                if(pfCluster_sim_fraction_min1_MatchedIndex->at(idxPF) != idxCalo) continue;

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

                /*if(logET >= -1.0 && logET < -0.9) logET_bin = 0;
                else if(logET >= -0.9 && logET < -0.8) logET_bin = 1;
                else if(logET >= -0.8 && logET < -0.7) logET_bin = 2;
                else if(logET >= -0.7 && logET < -0.6) logET_bin = 3;
                else if(logET >= -0.6 && logET < -0.5) logET_bin = 4;
                else if(logET >= -0.5 && logET < -0.4) logET_bin = 5;
                else if(logET >= -0.4 && logET < -0.3) logET_bin = 6;
                else if(logET >= -0.3 && logET < -0.2) logET_bin = 7;
                else if(logET >= -0.2 && logET < -0.1) logET_bin = 8;
                else if(logET >= -0.1 && logET < 0.0) logET_bin = 9;
                else if(logET >= 0.0 && logET < 0.1) logET_bin = 10;
                else if(logET >= 0.1 && logET < 0.2) logET_bin = 11;
                else if(logET >= 0.2 && logET < 0.3) logET_bin = 12;
                else if(logET >= 0.3 && logET < 0.4) logET_bin = 13;
                else if(logET >= 0.4 && logET < 0.5) logET_bin = 14;
                else if(logET >= 0.5 && logET < 0.6) logET_bin = 15;
                else if(logET >= 0.6 && logET < 0.7) logET_bin = 16;
                else if(logET >= 0.7 && logET < 0.8) logET_bin = 17;
                else if(logET >= 0.8 && logET < 0.9) logET_bin = 18;
                else if(logET >= 0.9 && logET < 1.0) logET_bin = 19;
                else if(logET >= 1.0 && logET < 1.1) logET_bin = 20;
                else if(logET >= 1.1 && logET < 1.2) logET_bin = 21;
                else if(logET >= 1.2 && logET < 1.3) logET_bin = 22;
                else if(logET >= 1.3 && logET < 1.4) logET_bin = 23;
                else if(logET >= 1.4 && logET < 1.5) logET_bin = 24;
                else if(logET >= 1.5 && logET < 1.6) logET_bin = 25;
                else if(logET >= 1.6 && logET < 1.7) logET_bin = 26;
                else if(logET >= 1.7 && logET < 1.8) logET_bin = 27;
                else if(logET >= 1.8 && logET < 1.9) logET_bin = 28;
                else if(logET >= 1.9 && logET < 2.0) logET_bin = 29;*/

                if(logET_bin < 0 || seed_eta_bin < 0) continue;

                caloClusters_shape[seed_eta_bin][logET_bin] -> Fill(dPhi, dEta, ET);
                cluster_dPhi_vs_loget[seed_dphi_bin] -> Fill(logET, fabs(dPhi), ET);
                
            }

        }

    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop"<<endl;
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
    //Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f)*(x[0]*x[0])) + std::max(b_upper, 0.0087f) + 0.0087;
    Double_t upper_curve = ((1 / (4 * a_upper)))*(x[0]*x[0]) + b_upper + 0.0087;

    //Double_t fitval = ( par[0]*par[2]*sin(par[2]) ) + (par[1] / sqrt(1.1 + x[0]));

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
    //Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f)*(x[0]*x[0])) + std::max(b_upper, 0.0087f) + 0.0087;
    Double_t lower_curve = ((1 / (4 * a_lower)))*(x[0]*x[0]) + b_lower;

    //Double_t fitval = ( par[0]*par[2]*sin(par[2]) ) + (par[1] / sqrt(1.1 + x[0]));

    return lower_curve;
}

void separation_opt(){
//optimization of parabola separation parameters
    vector<float> mid_etas = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9};
    //vector<float> mid_logets = {-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
    vector<float> mid_logets = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9};

    TGraph *envelope_plot_upper[15][30];
    TF1 *curve_fits_upper[15][30];

    TGraph *envelope_plot_lower[15][30];
    TF1 *curve_fits_lower[15][30];

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);


    for(int seed_eta_idx = 0; seed_eta_idx<15; seed_eta_idx++){ //loop over seed eta bins
        for(int loget_idx = 0; loget_idx<15; loget_idx++){  //loop over loget binsf
            int envelope_percent = 98,
                x_bins_total = caloClusters_shape[seed_eta_idx][loget_idx] -> GetNbinsX(),
                y_bins_total = caloClusters_shape[seed_eta_idx][loget_idx] -> GetNbinsY();

            int point_counter = 0;

            c1->Clear();
            envelope_plot_upper[seed_eta_idx][loget_idx] = new TGraph();//x_bins_total);
            envelope_plot_lower[seed_eta_idx][loget_idx] = new TGraph();//x_bins_total);

            //loop over x bins
            for(int xbin_idx=0; xbin_idx<x_bins_total+1; xbin_idx++){
                float xbin_integral = caloClusters_shape[seed_eta_idx][loget_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);

                float envelope_biny_upper = 0, envelope_biny_lower = 0, sum_temp_upper = 0, sum_temp_lower = 0;

                float xbin_percentage = ((float)envelope_percent * xbin_integral) / 100.0;

                if(xbin_integral <= 0.0) continue;

                //loop over y bins - upper
                for(int ybin_idx=0; ybin_idx<y_bins_total+1; ybin_idx++){
                    //float ybin_integral 
                    float curr_bin_content = caloClusters_shape[seed_eta_idx][loget_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_upper += curr_bin_content;

                    if(sum_temp_upper >= xbin_percentage){
                        envelope_biny_upper = ybin_idx;
                        break;
                    }
                } //y bin loop
                //loop over y bins - lower
                for(int ybin_idx=y_bins_total+1; ybin_idx>0; ybin_idx--){
                    //float ybin_integral 
                    float curr_bin_content = caloClusters_shape[seed_eta_idx][loget_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_lower += curr_bin_content;

                    if(sum_temp_lower >= xbin_percentage){
                        envelope_biny_lower = ybin_idx;
                        break;
                    }
                } //y bin loop
                float envelope_x_upper = caloClusters_shape[seed_eta_idx][loget_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                      envelope_y_upper = caloClusters_shape[seed_eta_idx][loget_idx] -> GetYaxis() -> GetBinCenter(envelope_biny_upper),
                      envelope_x_lower = caloClusters_shape[seed_eta_idx][loget_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                      envelope_y_lower = caloClusters_shape[seed_eta_idx][loget_idx] -> GetYaxis() -> GetBinCenter(envelope_biny_lower);

                //set the point
                envelope_plot_upper[seed_eta_idx][loget_idx]->SetPoint(point_counter++, envelope_x_upper, envelope_y_upper);
                envelope_plot_lower[seed_eta_idx][loget_idx]->SetPoint(point_counter++, envelope_x_lower, envelope_y_lower);
            } //x bin loop

            if(envelope_plot_upper[seed_eta_idx][loget_idx]->GetN() <= 5) continue;

            //fitting
            curve_fits_upper[seed_eta_idx][loget_idx] = new TF1(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str(),fitFunc_upper,-0.6,0.6,6);
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParameters(-0.00571429,-0.002,0.0135714,0.001,0.0,0.0);
            curve_fits_upper[seed_eta_idx][loget_idx] -> FixParameter(5, mid_etas.at(seed_eta_idx));
            curve_fits_upper[seed_eta_idx][loget_idx] -> FixParameter(4, mid_logets.at(loget_idx));

            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(0, "w00");
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(1, "w01");
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(2, "w10");
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(3, "w11");
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(5, "seedeta");
            curve_fits_upper[seed_eta_idx][loget_idx] -> SetParName(4, "loget");


            bool fitSuccess = false;
            double NDF = 0, 
                chiSqperNDF = 999,
                fitIteration = 0,
                chiSq = 999;
            string fitStatus = "";

            //fit minimization
            while(fitIteration < 5){
                envelope_plot_upper[seed_eta_idx][loget_idx] -> Fit(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str());

                fitStatus = gMinuit->fCstatu;
                chiSq = curve_fits_upper[seed_eta_idx][loget_idx]->GetChisquare();
                NDF = curve_fits_upper[seed_eta_idx][loget_idx]->GetNDF();

                if(NDF != 0) chiSqperNDF = chiSq / NDF;

                //cout << fitStatus <<endl;
                //if(fitStatus.compare("CONVERGED") || fitStatus.compare("OK") || fitStatus.compare("SUCCESSFUL")) fitSuccess = true;

                fitIteration++;
            }

            curve_fits_lower[seed_eta_idx][loget_idx] = new TF1(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str(),fitFunc_lower,-0.6,0.6,6);
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParameters(-0.00571429,-0.002,0.0135714,0.001,0.0,0.0);
            curve_fits_lower[seed_eta_idx][loget_idx] -> FixParameter(5,mid_etas.at(seed_eta_idx));
            curve_fits_lower[seed_eta_idx][loget_idx] -> FixParameter(4,mid_logets.at(loget_idx));

            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(0, "w00");
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(1, "w01");
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(2, "w10");
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(3, "w11");
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(5, "seedeta");
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetParName(4, "loget");



            fitSuccess = false;
            NDF = 0;
            chiSqperNDF = 999;
            fitIteration = 0;
            chiSq = 999;
            fitStatus = "";

            //fit minimization
            while(!fitSuccess){// || chiSqperNDF > 1) && fitIteration < 5){
                envelope_plot_lower[seed_eta_idx][loget_idx]->Fit(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)).c_str());
           
                fitStatus = gMinuit->fCstatu;
                chiSq = curve_fits_lower[seed_eta_idx][loget_idx]->GetChisquare();
                NDF = curve_fits_lower[seed_eta_idx][loget_idx]->GetNDF();

                if(NDF != 0) chiSqperNDF = chiSq / NDF;

                cout << "fit status: "<<fitStatus <<endl;
                if(fitStatus.compare("CONVERGED") || fitStatus.compare("OK") || fitStatus.compare("SUCCESSFUL")) fitSuccess = true;

                fitIteration++;
            }

            //plotting
            gStyle->SetOptStat(0);

            caloClusters_shape[seed_eta_idx][loget_idx] -> SetTitle((titles_etas[seed_eta_idx] + "     " + titles_loget[loget_idx]).c_str());

            envelope_plot_upper[seed_eta_idx][loget_idx] ->GetXaxis()->SetTitle("dPhi");
            envelope_plot_upper[seed_eta_idx][loget_idx] ->GetYaxis()->SetTitle("dEta");

            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerStyle(34);
            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerColor(2);
            envelope_plot_upper[seed_eta_idx][loget_idx] -> SetMarkerSize(3);

            envelope_plot_lower[seed_eta_idx][loget_idx] -> GetXaxis() -> SetTitle("dPhi");
            envelope_plot_lower[seed_eta_idx][loget_idx] -> GetYaxis() -> SetTitle("dEta");

            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerStyle(34);
            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerColor(6);
            envelope_plot_lower[seed_eta_idx][loget_idx] -> SetMarkerSize(3);

            curve_fits_upper[seed_eta_idx][loget_idx] -> SetLineWidth(4);
            curve_fits_lower[seed_eta_idx][loget_idx] -> SetLineWidth(4);

            curve_fits_lower[seed_eta_idx][loget_idx] -> SetLineColor(6);

            caloClusters_shape[seed_eta_idx][loget_idx] -> Draw("COLZ");
            envelope_plot_upper[seed_eta_idx][loget_idx] -> Draw("PSAME");
            envelope_plot_lower[seed_eta_idx][loget_idx] -> Draw("PSAME");
            curve_fits_upper[seed_eta_idx][loget_idx] -> Draw("SAME");
            curve_fits_lower[seed_eta_idx][loget_idx] -> Draw("SAME");

            c1->SaveAs(("envelope/upper_envelope_fit_"+to_string(seed_eta_idx)+"_"+to_string(loget_idx)+".png").c_str());

            //aggregate parameter values
            fit_w00_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(0));
            fit_w01_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(1));
            fit_w10_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(2));
            fit_w11_upper.push_back(curve_fits_upper[seed_eta_idx][loget_idx] -> GetParameter(3));

            fit_w00_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(0));
            fit_w01_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(1));
            fit_w10_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(2));
            fit_w11_lower.push_back(curve_fits_lower[seed_eta_idx][loget_idx] -> GetParameter(3));


        } //loget index
    } //seed eta index
}

Double_t fitFunc_dPhi(Double_t *x, Double_t *par){
//dynamix delta phi fit function definition
    float yoffset = par[0],
          scale   = par[1],
          xoffset = par[2],
          width   = par[3];

    float cutoff = 0.14,
          saturation = 0.60;

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

        int envelope_percent = 99,
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
                float curr_bin_content = cluster_dPhi_vs_loget[dphi_idx]->GetBinContent(xbin_idx, ybin_idx);
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
        dphi_fits[dphi_idx] -> SetParameters(yoffset_orig[dphi_idx], scale_orig[dphi_idx], xoffset_orig[dphi_idx], width_orig[dphi_idx]);

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
        while(fitIteration < 5){
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

        dphi_legends[dphi_idx]->AddEntry(dphi_fits[dphi_idx], "Fit DPhi Function","l");
        dphi_legends[dphi_idx]->AddEntry(dphi_original[dphi_idx], "Original DPhi Function","l");

        cluster_dPhi_vs_loget[dphi_idx] -> Scale(1.0 / cluster_dPhi_vs_loget[dphi_idx] -> Integral());

        cluster_dPhi_vs_loget[dphi_idx] -> Draw("COLZ");
        dphi_envelope_plot[dphi_idx] -> Draw("PSAME");
        dphi_fits[dphi_idx] -> Draw("SAME");
        dphi_original[dphi_idx] -> Draw("SAME");
        dphi_legends[dphi_idx]->Draw("SAME");

        c2 -> SaveAs(("dphi/dphi_fit_"+to_string(dphi_idx)+".png").c_str());

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
}

void envelope_optimization(){
//main program
    InitTree("root://eoscms.cern.ch///eos/cms/store/user/lzygala/par_con_files/finalDumperTree_4G_Run3_PU.root");
    EventLoop();
    //separation_opt();
    dPhi_opt();
    final_params();
}