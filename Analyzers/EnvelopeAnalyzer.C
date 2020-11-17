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
#include<TLine.h>
#include<TH2F.h>
#include<TPaveText.h>
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

#include<sys/stat.h>
#include<errno.h>

#ifdef __MAKECINT__ 
#pragma link C++ class vector<vector<double> >+;  
#pragma link C++ class vector<vector<float> >+; 
#pragma link C++ class vector<vector<int> >+;  
#pragma link C++ class vector<vector<bool> >+; 
#pragma link C++ class vector<vector<map<int,float>> >+; 
#endif


/* ---------- Variables ---------- */
    int seedEtaBins = 15,
        logEBins = 15,
        caloClusterShapeDEtaDistBins = 120,
        dPhiWindowEtaBins = 4,
        dPhiWindowETDistBins = 200;

    double minSeedEta = 0.0, 
           maxSeedEta = 3.0, 
           minLogE = -1.0, 
           maxLogE = 2.0;

    //ofstream outfile;
    ofstream data_outfile;
    ofstream fit_outfile;
    TFile* score_infile;

    TH2F* h2_Minimum_simScore;

    vector<TH2F*> cluster_dPhi_vs_loget;

    vector<vector<TH2F*>> caloClusters_shape_eBins_etWeight;
    vector<vector<TH2F*>> caloClusters_shape_eBins_eWeight;
    vector<vector<TH2F*>> caloClusters_shape_eBins;
    vector<vector<vector<TH1F*>>> caloClusters_eBins_dEtaDist;

    vector<vector<vector<TH1F*>>> caloClusters_dEtaDist;
    vector<vector<TH1F*>> caloClusters_dPhiDist;

    float pi=3.1415927;
    float twopi= 2*pi;
    float eta_crack_max = 1.566, eta_crack_min = 1.442;

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
    vector<vector<int> > *caloParticle_pfCluster_sim_fraction_MatchedIndex;
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
    vector<int>     *pfCluster_sim_fraction_MatchedIndex;
    vector<int>     *pfCluster_simScore_MatchedIndex;
    vector<vector<double> > *pfCluster_dR_genScore;
    vector<vector<double> > *pfCluster_dR_simScore;
    vector<vector<double> > *pfCluster_sim_fraction;
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

void setOutput(string path){
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1){
        if( errno == EEXIST ) {
        // alredy exists
            return;
        } else {
        // something else
            std::cout << "cannot create folder error:" << strerror(errno) << std::endl;
        }
    }
    return;
}

bool comparePFpt(const std::pair<int, double>&i, const std::pair<int, double>&j)
{
   return i.second > j.second;
}

float getThreshold(float seedEta, float seedET){
    int xBin = h2_Minimum_simScore -> GetXaxis() -> FindBin(seedET),
        yBin = h2_Minimum_simScore -> GetYaxis() -> FindBin(fabs(seedEta));

    float scoreThreshold = h2_Minimum_simScore -> GetBinContent(xBin, yBin);

    return scoreThreshold;
}

void EventLoop(){
//Loops over all events and fills histograms
    int eventLoopMax = 1000;//EvMax;

    for(int iev=0; iev<eventLoopMax;++iev){
        EventTree->GetEvent(iev);
        if(iev % 1000 == 0) cout<<cur_time()<<"\tProcessing event: "<<iev <<" / "<<eventLoopMax<<endl;

        //get all pfClusters matched to caloParticles
        std::vector< std::vector< std::pair<int, double> > > pfMatched_idx_pt(caloParticle_id->size());
        std::vector< std::vector< std::pair<int, double> > > pfMatched_idx_pt_temp(caloParticle_id->size());
        std::vector< std::vector< std::pair<int, double> > > pfMatched_idx_pt_simFrac(caloParticle_id->size());
        for(unsigned int iPF=0; iPF<pfCluster_energy->size(); iPF++){
            int caloMatch = -1;//pfCluster_sim_fraction_MatchedIndex->at(iPF);//pfCluster_simScore_MatchedIndex->at(iPF);
            int maxScore = std::max_element(pfCluster_sim_fraction->at(iPF).begin(), pfCluster_sim_fraction->at(iPF).end()) - pfCluster_sim_fraction->at(iPF).begin();

            int rechitseed = std::max_element(pfClusterHit_rechitEnergy->at(iPF).begin(), pfClusterHit_rechitEnergy->at(iPF).end()) - pfClusterHit_rechitEnergy->at(iPF).begin();
            float seedEta = pfClusterHit_eta->at(iPF).at(rechitseed);
            float seedET = pfClusterHit_rechitEnergy->at(iPF).at(rechitseed)/TMath::CosH(pfClusterHit_eta->at(iPF).at(rechitseed));

            float simScoreThresh = getThreshold(seedEta, seedET);
            if(pfCluster_sim_fraction->at(iPF).at(maxScore) > 0.01) caloMatch = maxScore;
            if(caloMatch >= 0)
                pfMatched_idx_pt_temp.at(caloMatch).push_back(std::make_pair(iPF, pfCluster_energy->at(iPF)/TMath::CosH(pfCluster_eta->at(iPF))));
        }
        //loop over caloParticles
        for(unsigned int idxCalo=0; idxCalo<caloParticle_id->size(); idxCalo++){
            float genEta = caloParticle_genEta->at(idxCalo);
            //if(fabs(genEta) < eta_crack_max && fabs(genEta)>eta_crack_min) continue;

            //order matched pfClusters in pt
            sort(pfMatched_idx_pt_temp.at(idxCalo).begin(),pfMatched_idx_pt_temp.at(idxCalo).end(),comparePFpt);

            if(pfMatched_idx_pt_temp.at(idxCalo).size() == 0) continue;

            int seedCluster = pfMatched_idx_pt_temp.at(idxCalo).at(0).first;
            float seedEta = pfCluster_eta->at(seedCluster),
                  seedET = pfMatched_idx_pt_temp.at(idxCalo).at(0).second;
            int dPhi_bin = -1, seedEta_bin = -1;

            float simScoreThresh = getThreshold(seedEta, seedET);

            for(unsigned int iPF=0; iPF<pfCluster_energy->size(); iPF++){
                int caloMatch = std::max_element(pfCluster_sim_fraction->at(iPF).begin(), pfCluster_sim_fraction->at(iPF).end()) - pfCluster_sim_fraction->at(iPF).begin();
                if(caloMatch == (int)idxCalo && pfCluster_sim_fraction->at(iPF).at(caloMatch) >= simScoreThresh)
                    pfMatched_idx_pt.at(caloMatch).push_back(std::make_pair(iPF, pfCluster_energy->at(iPF)/TMath::CosH(pfCluster_eta->at(iPF))));  
            }

            sort(pfMatched_idx_pt.at(idxCalo).begin(),pfMatched_idx_pt.at(idxCalo).end(),comparePFpt);


            double seedEtaStep = (maxSeedEta - minSeedEta) / (double)seedEtaBins;
            double logEStep = (maxLogE - minLogE) / (double)logEBins;

            double seedEtaVal = minSeedEta;
            for(int iSeedEta = 0; iSeedEta < seedEtaBins; iSeedEta++){
                if(fabs(seedEta) >= seedEtaVal && fabs(seedEta) < seedEtaVal + seedEtaStep){
                    seedEta_bin = iSeedEta;
                    break;
                } 
                seedEtaVal+=seedEtaStep;
            }

            if     (fabs(seedEta) >= 0.0   && fabs(seedEta) < 1.479){ dPhi_bin = 0; }
            else if(fabs(seedEta) >= 1.479 && fabs(seedEta) < 1.75) { dPhi_bin = 1; }
            else if(fabs(seedEta) >= 1.75  && fabs(seedEta) < 2.0)  { dPhi_bin = 2; } 
            else if(fabs(seedEta) >= 2.0   && fabs(seedEta) < 3.0)  { dPhi_bin = 3; }

            //loop over all matched pfClusters
            for(auto pfMatch_idx_pt_el : pfMatched_idx_pt_temp.at(idxCalo)){
                int idxPF = pfMatch_idx_pt_el.first;
                double etPF = pfMatch_idx_pt_el.second;
                double ePF = pfCluster_energy->at(idxPF);

                //skip the seeds
                if(idxPF == seedCluster) continue;

                float dEta = (1-2*(pfCluster_eta->at(seedCluster) < 0)) * (pfCluster_eta->at(idxPF) - pfCluster_eta->at(seedCluster));
                float dPhi = pfCluster_phi->at(idxPF) - pfCluster_phi->at(seedCluster);
                if(dPhi > pi) dPhi -=twopi;
                if(dPhi < -pi) dPhi += twopi;
                float logET = log10(etPF);
                float logE = log10(ePF);

                int logE_bin = -1;

                double logEVal = minLogE;
                for(int iLogE = 0; iLogE < logEBins; iLogE++){
                    if(logE >= logEVal && logE < logEVal + logEStep){
                        logE_bin = iLogE;
                        break;
                    } 
                    logEVal+=logEStep;
                }


                if(logE_bin < 0 || seedEta_bin < 0){cout<<logE_bin<<"\t"<<seedEta_bin<<endl; continue;}

                //fill histograms
                caloClusters_shape_eBins_etWeight[seedEta_bin][logE_bin] -> Fill(dPhi, dEta, etPF);
                caloClusters_shape_eBins_eWeight[seedEta_bin][logE_bin] -> Fill(dPhi, dEta, ePF);
                caloClusters_shape_eBins[seedEta_bin][logE_bin] -> Fill(dPhi, dEta);

                //index 0 : underflow bin
                //index 1->caloClusterShapeDEtaDistBins : valid bins
                //index caloClusterShapeDEtaDistBins + 1 : overflow bin
                int dPhiBin = caloClusters_shape_eBins_etWeight[seedEta_bin][logE_bin] -> GetXaxis() -> FindBin(dPhi);
                if(dPhiBin < caloClusterShapeDEtaDistBins) 
                    caloClusters_eBins_dEtaDist[seedEta_bin][logE_bin][dPhiBin] -> Fill(dEta, etPF);

                int logETBin = cluster_dPhi_vs_loget[dPhi_bin] -> GetXaxis() -> FindBin(logET);
                if(logETBin < dPhiWindowETDistBins)
                    caloClusters_dPhiDist[dPhi_bin][logETBin] -> Fill(fabs(dPhi), etPF);
            }
        }// loop over calo particels
    }//event loop
    cout<<cur_time()<<"\tFinished Event Loop 2"<<endl;
}

void InitHistograms(){

    fit_outfile.open("envelope/fit_status.txt", std::ofstream::out | std::ofstream::trunc);
    score_infile = TFile::Open("Data/simScore_Minima_withHitFraction.root");
    score_infile->GetObject("h2_Minimum_simScore", h2_Minimum_simScore);

    double min_dEta = -0.15,
           max_dEta = 0.2,
           min_dPhi = -0.6,
           max_dPhi = 0.6,
           min_dEta_large = -0.5,
           max_dEta_large = 0.5;
    
    int bins_dEta = 105,
        bins_dPhi = 120,
        bins_dEta_large = 300;

    caloClusters_shape_eBins.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_shape_eBins_etWeight.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_shape_eBins_eWeight.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_eBins_dEtaDist.resize(seedEtaBins, vector<vector<TH1F*>>(logEBins, vector<TH1F*>(caloClusterShapeDEtaDistBins + 2)));

    cluster_dPhi_vs_loget.resize(dPhiWindowEtaBins);
    caloClusters_dPhiDist.resize(dPhiWindowEtaBins, vector<TH1F*>(dPhiWindowETDistBins + 2));

    for(int seedEtaIdx = 0; seedEtaIdx < seedEtaBins; seedEtaIdx++){
        for(int logEIdx = 0; logEIdx < logEBins; logEIdx++){    
            caloClusters_shape_eBins[seedEtaIdx][logEIdx] = new TH2F(("caloClusters_shape_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),("caloClusters_shape_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_eBins_etWeight[seedEtaIdx][logEIdx] = new TH2F(("caloClusters_shape_eBins_etWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),("caloClusters_shape_eBins_etWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            caloClusters_shape_eBins_eWeight[seedEtaIdx][logEIdx] = new TH2F(("caloClusters_shape_eBins_eWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),("caloClusters_shape_eBins_eWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),bins_dPhi,min_dPhi,max_dPhi,bins_dEta,min_dEta,max_dEta);
            
            for(int k = 0; k < caloClusterShapeDEtaDistBins + 2; k++){
                caloClusters_eBins_dEtaDist[seedEtaIdx][logEIdx][k] = new TH1F(("ET_vs_dEta_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)+"_"+to_string(k)).c_str(),("ET_vs_dEta_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)+"_"+to_string(k)).c_str(), 105, -0.15, 0.2);
            }
        }
    }
    for(int ii = 0; ii < dPhiWindowEtaBins; ii++){
        cluster_dPhi_vs_loget[ii] = new TH2F(("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),200, -1.0, 3.0, 200, 0.0, 1.0);
        
        for(int jj = 0; jj < dPhiWindowETDistBins + 2; jj++){
            caloClusters_dPhiDist[ii][jj] = new TH1F(("ET_vs_dPhi_"+to_string(ii)+"_"+to_string(jj)).c_str(),("ET_vs_dPhi_"+to_string(ii)+"_"+to_string(jj)).c_str(), 200, 0.0, 1.0);
        }
    }
}

void SaveHistograms(string outputFile){
    cout<<cur_time()<<"\tCreating Output File"<<endl;
    TFile *plotFile = new TFile(outputFile.c_str(), "RECREATE");
    plotFile->cd();
    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    for(int seedEta_bin = 0; seedEta_bin < seedEtaBins; seedEta_bin++){
        for(int logET_bin = 0; logET_bin < logEBins; logET_bin++){

            caloClusters_shape_eBins[seedEta_bin][logET_bin]->Draw();
            caloClusters_shape_eBins[seedEta_bin][logET_bin]->Write();

            caloClusters_shape_eBins_eWeight[seedEta_bin][logET_bin]->Draw();
            caloClusters_shape_eBins_eWeight[seedEta_bin][logET_bin]->Write();

            caloClusters_shape_eBins_etWeight[seedEta_bin][logET_bin]->Draw();
            caloClusters_shape_eBins_etWeight[seedEta_bin][logET_bin]->Write();

            for(int dPhi = 0; dPhi < caloClusterShapeDEtaDistBins + 2; dPhi++){
                caloClusters_eBins_dEtaDist[seedEta_bin][logET_bin][dPhi]->Draw();
                caloClusters_eBins_dEtaDist[seedEta_bin][logET_bin][dPhi]->Write();
            }
        }
    }

    for(int dPhi = 0; dPhi < dPhiWindowEtaBins; dPhi++){
        cluster_dPhi_vs_loget[dPhi] -> Draw();
        cluster_dPhi_vs_loget[dPhi] -> Write();
        for(int x = 0; x < dPhiWindowETDistBins + 2; x++){
            caloClusters_dPhiDist[dPhi][x] ->Draw();
            caloClusters_dPhiDist[dPhi][x] ->Write();
        }
    }

   plotFile->Close(); 
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
    EventTree->SetBranchAddress("caloParticle_superCluster_dR_simScore_MatchedIndex", &caloParticle_superCluster_dR_simScore_MatchedIndex);

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
    EventTree->SetBranchAddress("pfCluster_sim_fraction_MatchedIndex", &pfCluster_sim_fraction_MatchedIndex);//_pfCluster_sim_fraction_old_MatchedIndex);
    //EventTree->SetBranchAddress("pfCluster_simScore_MatchedIndex", &pfCluster_simScore_MatchedIndex);//_pfCluster_simScore_MatchedIndex);
    EventTree->SetBranchAddress("pfCluster_dR_genScore", &pfCluster_dR_genScore);//_pfCluster_dR_genScore);
    EventTree->SetBranchAddress("pfCluster_dR_simScore", &pfCluster_dR_simScore);//_pfCluster_dR_simScore);
    EventTree->SetBranchAddress("pfCluster_sim_fraction", &pfCluster_sim_fraction);//_pfCluster_sim_fraction_old);
    //EventTree->SetBranchAddress("pfCluster_simScore", &pfCluster_simScore);//_pfCluster_simScore);
    //EventTree->SetBranchAddress("pfClusterHit_energy", &pfClusterHit_energy);//_pfClusterHit_energy);
    EventTree->SetBranchAddress("pfClusterHit_rechitEnergy", &pfClusterHit_rechitEnergy);//_pfClusterHit_rechitEnergy);
    EventTree->SetBranchAddress("pfClusterHit_eta", &pfClusterHit_eta);//_pfClusterHit_eta);
    EventTree->SetBranchAddress("pfClusterHit_phi", &pfClusterHit_phi);//_pfClusterHit_phi);
    EventTree->SetBranchAddress("pfClusterHit_ieta", &pfClusterHit_ieta);//_pfClusterHit_ieta);
    EventTree->SetBranchAddress("pfClusterHit_iphi", &pfClusterHit_iphi);//_pfClusterHit_iphi);
    EventTree->SetBranchAddress("pfClusterHit_iz", &pfClusterHit_iz);//_pfClusterHit_iz);

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

    EvMax=EventTree->GetEntries();

    cout<<cur_time()<<"\tTree successfully processed!\n";

    EventLoop();
    infile->Close();
    //score_infile->Close();

    return;
}

void EnvelopeAnalyzer(string inputFile, string outputFile){
//main program

    InitHistograms();
    InitTree(inputFile.c_str()); 
    SaveHistograms(outputFile.c_str());

    fit_outfile.close();
}