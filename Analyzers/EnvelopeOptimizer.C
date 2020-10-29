
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
#include<sstream>
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
    //binning
    int seedEtaBins = 15,
        logEBins = 15,
        caloClusterShapeDEtaDistBins = 120,
        dPhiWindowEtaBins = 4,
        dPhiWindowETDistBins = 200;

    double minSeedEta = 0.0, 
           maxSeedEta = 3.0, 
           minLogE = -2.0, 
           maxLogE = 1.0;

    //ofstream outfile;
    ofstream data_outfile;
    ofstream fit_outfile;
    ofstream aveParam_outfile;
    ofstream localParam_outfile;
    ofstream dPhiParam_outfile;
    TFile* score_infile;

    vector<TH2F*> cluster_dPhi_vs_loget;            //[dPhiWindowEtaBins]
    vector<vector<TH1F*>> caloClusters_dPhiDist;    //[dPhiWindowEtaBins][dPhiWindowETDistBins]

    vector<vector<TH2F*>> caloClusters_shape_eBins_etWeight;    //[seedEtaBins][logEBins]
    vector<vector<TH2F*>> caloClusters_shape_eBins_eWeight;     //[seedEtaBins][logEBins]
    vector<vector<TH2F*>> caloClusters_shape_eBins;             //[seedEtaBins][logEBins]
    vector<vector<vector<TH1F*>>> caloClusters_eBins_dEtaDist;  //[seedEtaBins][logEBins][caloClusterShapeDEtaDistBins]

    vector<vector<bool>> fitPlot;   //[seedEtaBins][logEBins]

    vector<vector<float>> fit_w00_up,
                          fit_w01_up,
                          fit_w10_up,
                          fit_w11_up,
                          fit_w00_low,
                          fit_w01_low,
                          fit_w10_low,
                          fit_w11_low;

    float curr_logET;
    float curr_seedEta;

    //current values
    float original_w00 = -0.00571429;
    float original_w01 = -0.002;
    float original_w10 = 0.0135714;
    float original_w11 = 0.001;

    float phase1_sepPars[4] = {-0.00571429,-0.002,0.0135714,0.001};

    float original_p00 = -0.107537;
	float original_p01 = 0.590969;
	float original_p02 = -0.076494;
	float original_p10 = -0.0268843;
	float original_p11 = 0.147742;
	float original_p12 = -0.0191235;

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

    int weightNum = 5;
        //0 = Unweighted
        //1 = Total Clusters
        //2 = Cluster ET
        //3 = Cluster Energy
        //4 = Eta

    string weightTitles[] = { "UNWEIGHTED","TOTALCLUSTER", "CLUSTERET", "CLUSTEREN", "ETA" };

    vector<vector<float>> fit_w00_upper_weighted(weightNum);
    vector<vector<float>> fit_w01_upper_weighted(weightNum);
    vector<vector<float>> fit_w10_upper_weighted(weightNum);
    vector<vector<float>> fit_w11_upper_weighted(weightNum);
    vector<vector<float>> fit_w00_lower_weighted(weightNum);
    vector<vector<float>> fit_w01_lower_weighted(weightNum);
    vector<vector<float>> fit_w10_lower_weighted(weightNum);
    vector<vector<float>> fit_w11_lower_weighted(weightNum);

    vector<vector<float>> fit_w00_upper_weighted_denom(weightNum);
    vector<vector<float>> fit_w01_upper_weighted_denom(weightNum);
    vector<vector<float>> fit_w10_upper_weighted_denom(weightNum);
    vector<vector<float>> fit_w11_upper_weighted_denom(weightNum);
    vector<vector<float>> fit_w00_lower_weighted_denom(weightNum);
    vector<vector<float>> fit_w01_lower_weighted_denom(weightNum);
    vector<vector<float>> fit_w10_lower_weighted_denom(weightNum);
    vector<vector<float>> fit_w11_lower_weighted_denom(weightNum);
    
    vector<string> titles_etas, titles_loge, filename_etas, filename_loges;

    vector<float> mid_etas,
                  mid_loges;

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

Double_t fitFunc_upper(Double_t *x, Double_t *par){
//upper parabola fit function definition
    float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = curr_logET,//par[4],
          seedeta = curr_seedEta;//par[5];

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

    //Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f))*(x[0]*x[0]) + std::max(b_upper, 0.0087f) + 0.0087;
    Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f))*(x[0]*x[0]) + b_upper + 0.0087;

    return upper_curve;
}

Double_t fitFunc_lower(Double_t *x, Double_t *par){
//upper parabola fit function definition
    float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = curr_logET,//par[4],
          seedeta = curr_seedEta;//par[5];

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

    //Double_t lower_curve = (std::max((1 / (4 * a_lower)),0.0f))*(x[0]*x[0]) + std::min(b_lower, -0.0087f);
    Double_t lower_curve = (std::max((1 / (4 * a_lower)),0.0f))*(x[0]*x[0]) + b_lower;

    return lower_curve;
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

float upper_parabola_dEtaValue(float w00, float w01, float w10, float w11, float loget, float seedeta, float dPhi){
    /*float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = curr_logET,//par[4],
          seedeta = curr_seedEta;//par[5];
    */
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

    float upper_curve = (std::max((1 / (4 * a_upper)),0.0f))*(dPhi*dPhi) + std::max(b_upper, 0.0087f) + 0.0087;
    //float upper_curve = (std::max((1 / (4 * a_upper)),0.0f))*(dPhi*dPhi) + b_upper + 0.0087;
    //cout<<"\tup"<<upper_curve<<endl;
    return upper_curve;
}

float lower_parabola_dEtaValue(float w00, float w01, float w10, float w11, float loget, float seedeta, float dPhi){
    /*float w00     = par[0],
          w01     = par[1],
          w10     = par[2],
          w11     = par[3],
          loget   = curr_logET,//par[4],
          seedeta = curr_seedEta;//par[5];
    */

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

    float lower_curve = (std::max((1 / (4 * a_lower)),0.0f))*(dPhi*dPhi) + std::min(b_lower, -0.0087f);
    //Double_t lower_curve = (std::max((1 / (4 * a_lower)),0.0f))*(x[0]*x[0]) + b_lower;
    //cout<<"\tlow"<<lower_curve<<endl;

    return lower_curve;
}

vector<float> XErrors(int xBin, int yBin, int etaBin, int logetBin){
    vector<float> errors;
    errors.push_back(0.0);  //no errors

    //width errors
    float sum = 1;//caloClusters_shape_xAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    float total_entries = 1;//caloClusters_shape_xAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    errors.push_back(sum / total_entries);

    //statistical errors
    errors.push_back(0.0); 

    //effective sigma errors
    errors.push_back(0.0); 

    return errors;
}

vector<float> YErrors(int xBin, int yBin, int etaBin, int logetBin, float dPhiVal=0.0){
    vector<float> errors;
    errors.push_back(0);  //no errors

    //width errors
    float sum = 0;//caloClusters_shape_yAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    float total_entries = 1;//caloClusters_shape_yAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    errors.push_back(sum / total_entries);

    //statistical errors
    int nEntries = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->GetEntries();
    float stdDev = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->GetStdDev();

    errors.push_back(stdDev / (pow(nEntries, 0.25)));

    //effective sigma errors
    int dEtaBins_total = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin] -> GetNbinsX();
        //nEntries = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->GetEntries();

    double lower_percentile = 15.85, upper_percentile = 84.15;
    float lower_point = 0.0, upper_point = 0.0;
    int lower_point_idx = 0, upper_point_idx = 0;
    bool lower_point_found = false, upper_point_found = false;

    //caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->GetQuantiles(1,&lower_point, &lower_percentile);
    //caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->GetQuantiles(1,&upper_point, &upper_percentile);
    float sum_temp = 0;
    float dEta_Integral = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin]->Integral();
    float lower_percentage = (lower_percentile * dEta_Integral) / 100.0,
          upper_percentage = (upper_percentile * dEta_Integral) / 100.0;

    for(int dEtabin_idx=1; dEtabin_idx<dEtaBins_total+1; dEtabin_idx++){
                    
        float curr_bin_content = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin] -> GetBinContent(dEtabin_idx);
        sum_temp += curr_bin_content;

        if(sum_temp >= lower_percentage && !lower_point_found){
            lower_point_idx = dEtabin_idx;
            lower_point_found = true;
        }
        if(sum_temp >= upper_percentage && !upper_point_found){
            upper_point_idx = dEtabin_idx;
            upper_point_found = true;
        }
        if(lower_point_found && upper_point_found) break;
    } //y bin loop

    lower_point = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin] -> GetBinCenter(lower_point_idx);
    upper_point = caloClusters_eBins_dEtaDist[etaBin][logetBin][xBin] -> GetBinCenter(upper_point_idx);

    float mean = (upper_point + lower_point) / 2;
    float width = upper_point - lower_point;
    float sigma_eff = upper_point - lower_point;
    
    errors.push_back(fabs(2 * sigma_eff / (pow(nEntries, 0.25))));


    return errors;
}

vector<float> XErrors_DPhi(int xBin, int yBin, int dPhiBin){
    vector<float> errors;
    errors.push_back(0.0);  //no errors

    //width errors
    //float sum = caloClusters_shape_xAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    //float total_entries = caloClusters_shape_xAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    errors.push_back(0);
    //errors.push_back(sum / total_entries);

    //statistical errors
    errors.push_back(0.0); 

    //effective sigma errors
    errors.push_back(0.0); 

    return errors;
}

vector<float> YErrors_DPhi(int xBin, int yBin, int dPhiBin, float dPhiVal=0.0){
    vector<float> errors;
    errors.push_back(0);  //no errors

    //width errors
    //float sum = caloClusters_shape_yAxisWeight_Numerator[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);
    //float total_entries = caloClusters_shape_yAxisWeight_Denom[etaBin][logetBin] -> Integral(xBin, xBin, yBin, yBin);

    errors.push_back(0);
    //errors.push_back(sum / total_entries);

    //statistical errors
    int nEntries = caloClusters_dPhiDist[dPhiBin][xBin]->GetEntries();
    float stdDev = caloClusters_dPhiDist[dPhiBin][xBin]->GetStdDev();

    errors.push_back(stdDev / (pow(nEntries, 0.25)));

    //effective sigma errors
    int dEtaBins_total = caloClusters_dPhiDist[dPhiBin][xBin] -> GetNbinsX();

    double lower_percentile = 15.85, upper_percentile = 84.15;
    float lower_point = 0.0, upper_point = 0.0;
    int lower_point_idx = 0, upper_point_idx = 0;
    bool lower_point_found = false, upper_point_found = false;

    float sum_temp = 0;
    float dEta_Integral = caloClusters_dPhiDist[dPhiBin][xBin]->Integral();
    float lower_percentage = (lower_percentile * dEta_Integral) / 100.0,
          upper_percentage = (upper_percentile * dEta_Integral) / 100.0;

    for(int dEtabin_idx=1; dEtabin_idx<dEtaBins_total+1; dEtabin_idx++){
                    
        float curr_bin_content = caloClusters_dPhiDist[dPhiBin][xBin] -> GetBinContent(dEtabin_idx);
        sum_temp += curr_bin_content;

        if(sum_temp >= lower_percentage && !lower_point_found){
            lower_point_idx = dEtabin_idx;
            lower_point_found = true;
        }
        if(sum_temp >= upper_percentage && !upper_point_found){
            upper_point_idx = dEtabin_idx;
            upper_point_found = true;
        }
        if(lower_point_found && upper_point_found) break;
    } //y bin loop

    lower_point = caloClusters_dPhiDist[dPhiBin][xBin] -> GetBinCenter(lower_point_idx);
    upper_point = caloClusters_dPhiDist[dPhiBin][xBin] -> GetBinCenter(upper_point_idx);

    float mean = (upper_point + lower_point) / 2;
    float width = upper_point - lower_point;
    float sigma_eff = upper_point - lower_point;
    
    errors.push_back(fabs(2 * sigma_eff / (pow(nEntries, 0.25))));


    return errors;
}

TF1* fitMinimization(TF1*& function, TGraphErrors*& points, string fit_name, int seed_eta_idx, int loge_idx, string fit_label, string& fitStatus, bool& fitSuccess){

    fitSuccess = false;
    double NDF = 0, 
        chiSqperNDF = 999,
        fitIteration = 0,
        chiSq = 999;
    fitStatus = "";

    while((!fitSuccess) && fitIteration < 100){
        fitStatus.clear();

        points -> Fit((fit_name).c_str(),"RNQ");

        fitStatus = gMinuit->fCstatu;
        fitStatus = trimString(fitStatus);
        chiSq = function->GetChisquare();
        NDF = function->GetNDF();

        if(NDF != 0) chiSqperNDF = chiSq / NDF;

        if((fitStatus.compare("CONVERGED") == 0))
        { fitSuccess = true; }

        fitIteration++;
    }

        cout<<fit_label<<":\tEta: "<<seed_eta_idx<<"\tloget: "<<loge_idx<<"\titeration: "<<fitIteration<<endl;
        cout<<fit_label<<":\tIterations: "<<fitIteration<<" Status: "<<fitStatus<<" Chi^2: "<<chiSq << " NDF: "<<NDF<<endl<<endl;
    fit_outfile<<fit_label<<":\tIterations: "<<fitIteration<<" Status: "<<fitStatus<<" Chi^2: "<<chiSq << " NDF: "<<NDF<<endl;

    return function;

}

void separationOptimization(bool fit, int error=0, bool ringRejection=0, bool binRejection=0){
//optimization of parabola separation parameters

    TGraphErrors *envelope_plot_upper[15][15];
    TF1 *curve_fits_upper[15][15];

    TGraphErrors *envelope_plot_lower[15][15];
    TF1 *curve_fits_lower[15][15];

    TLegend* legends[15][15]; 

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    localParam_outfile << "SEED_ETA\tLOG(E)\tPARABOLA\tPARAMETER\tVALUE\tFROM_FIT"<<endl;

    //loop over seed eta bins - 15    
    for(int seed_eta_idx = 0; seed_eta_idx < 15; seed_eta_idx++){ 
        //loop over loget bins - 15
        for(int loge_idx = 0; loge_idx < 15; loge_idx++){  

            int envelope_percent = 98,
                x_bins_total = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetNbinsX(),
                y_bins_total = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetNbinsY();

            int point_counter = 0;

            vector<float> x_points_upper, x_points_lower, y_points_upper, y_points_lower, 
                          x_errors_upper, y_errors_upper, x_errors_lower, y_errors_lower;

            c1->Clear();

            int total_clus = caloClusters_shape_eBins[seed_eta_idx][loge_idx] -> Integral();
            float total_clus_et = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> Integral();
            float total_clus_e = caloClusters_shape_eBins_eWeight[seed_eta_idx][loge_idx] -> Integral();

            //loop over x bins
            for(int xbin_idx = 1; xbin_idx < x_bins_total+1; xbin_idx++){
                float xbin_integral = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);
                int n_entries = caloClusters_shape_eBins[seed_eta_idx][loge_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);

                if(binRejection){ if(n_entries <= 4) continue; }
                
                int ybin_idx_upper = -1, ybin_idx_lower = -1;
                float sum_temp_upper = 0, sum_temp_lower = 0;

                float xbin_percentage = ((float)envelope_percent * xbin_integral) / 100.0;

                if(xbin_integral <= 0.0) continue;

                float xBin_val = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> GetBinCenter(xbin_idx);

                if(ringRejection){
                    if(xBin_val > -0.05 && xBin_val < 0.05 && loge_idx <8) continue; //skip the ring
                    if(seed_eta_idx >= 6 && xBin_val > -0.065 && xBin_val < 0.065 && loge_idx <8) continue;
                    if(seed_eta_idx >= 8 && xBin_val > -0.1 && xBin_val < 0.1 && loge_idx <8) continue;
                    if(seed_eta_idx >= 11 && xBin_val > -0.15 && xBin_val < 0.15) continue;
                }
                //loop over y bins - upper parabola point
                for(int ybin_idx=1; ybin_idx<y_bins_total+1; ybin_idx++){
                    
                    float curr_bin_content = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_upper += curr_bin_content;

                    if(sum_temp_upper >= xbin_percentage){
                        ybin_idx_upper = ybin_idx;
                        break;
                    }
                } //y bin loop
                //loop over y bins - lower parabola point
                for(int ybin_idx = y_bins_total + 1; ybin_idx > 0; ybin_idx--){
                    
                    float curr_bin_content = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetBinContent(xbin_idx, ybin_idx);
                    sum_temp_lower += curr_bin_content;

                    if(sum_temp_lower >= xbin_percentage){
                        ybin_idx_lower = ybin_idx;
                        break;
                    }
                } //y bin loop

                float yBin_val_up = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_upper);
                float yBin_val_low = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_lower);

                if(ybin_idx_upper == ybin_idx_lower) continue;
                if(ybin_idx_upper - ybin_idx_lower == 1) continue;

                vector<float> x_up_err   = XErrors(xbin_idx, ybin_idx_upper, seed_eta_idx, loge_idx),
                              x_down_err = XErrors(xbin_idx, ybin_idx_lower, seed_eta_idx, loge_idx),
                              y_up_err   = YErrors(xbin_idx, ybin_idx_upper, seed_eta_idx, loge_idx, xBin_val),
                              y_down_err = YErrors(xbin_idx, ybin_idx_lower, seed_eta_idx, loge_idx, xBin_val);

                if (ybin_idx_upper != -1){
                    float envelope_x_upper = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                          envelope_y_upper = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_upper);

                    x_points_upper.push_back(envelope_x_upper);
                    x_errors_upper.push_back(x_up_err[error]);
                    
                    y_points_upper.push_back(envelope_y_upper);
                    y_errors_upper.push_back(y_up_err[error]);
                }
                if (ybin_idx_lower != -1){
                    float envelope_x_lower = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> GetBinCenter(xbin_idx),
                          envelope_y_lower = caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> GetBinCenter(ybin_idx_lower);

                    x_points_lower.push_back(envelope_x_lower);
                    x_errors_lower.push_back(x_down_err[error]);

                    y_points_lower.push_back(envelope_y_lower);
                    y_errors_lower.push_back(y_down_err[error]);
                }
            } //x bin loop

            //fill graphs for fitting
            envelope_plot_upper[seed_eta_idx][loge_idx] = new TGraphErrors(x_points_upper.size(), &x_points_upper[0], &y_points_upper[0], &x_errors_upper[0], &y_errors_upper[0]);
            envelope_plot_lower[seed_eta_idx][loge_idx] = new TGraphErrors(x_points_lower.size(), &x_points_lower[0], &y_points_lower[0], &x_errors_lower[0], &y_errors_lower[0]);//x_bins_total);

            if(envelope_plot_upper[seed_eta_idx][loge_idx]->GetN() <= 7 || envelope_plot_lower[seed_eta_idx][loge_idx]->GetN() <= 7){
                //empty bins retain phase 1 parameters (not considered for averaged parameters)
                localParam_outfile << mid_etas[seed_eta_idx] << "\t" 
                                    << mid_loges[loge_idx] << "\t"
                                    << "UP" << "\t"
                                    << "PHASE1" << "\t"
                                    << phase1_sepPars[0] << "\t"
                                    << phase1_sepPars[1] << "\t"
                                    << phase1_sepPars[2] << "\t"
                                    << phase1_sepPars[3] << endl
                                    << mid_etas[seed_eta_idx] << "\t" 
                                    << mid_loges[loge_idx] << "\t"
                                    << "LOW" << "\t"
                                    << "PHASE1" << "\t"
                                    << phase1_sepPars[0] << "\t"
                                    << phase1_sepPars[1] << "\t"
                                    << phase1_sepPars[2] << "\t"
                                    << phase1_sepPars[3] << endl;
                continue;
            }
            curr_logET = mid_loges.at(loge_idx);
            curr_seedEta = mid_etas.at(seed_eta_idx);

            //fitting
            if(fit){
                fitPlot[seed_eta_idx][loge_idx] = true;

                curve_fits_upper[seed_eta_idx][loge_idx] = new TF1(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),fitFunc_upper,-0.6,0.6,4);
                curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameters(original_w00,original_w01,original_w10,original_w11);

                curve_fits_upper[seed_eta_idx][loge_idx] -> SetParName(0, "W00");
                curve_fits_upper[seed_eta_idx][loge_idx] -> SetParName(1, "W01");
                curve_fits_upper[seed_eta_idx][loge_idx] -> SetParName(2, "W10");
                curve_fits_upper[seed_eta_idx][loge_idx] -> SetParName(3, "W11");

                fit_outfile<<"Eta: "<<seed_eta_idx<<" LogE: "<<loge_idx<<endl;

                string fitStatus = "";
                bool fitSuccess = false;
                curve_fits_upper[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_upper[seed_eta_idx][loge_idx], envelope_plot_upper[seed_eta_idx][loge_idx], ("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "UPPER_FIT", fitStatus, fitSuccess);

                if(!fitSuccess){
                    
                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(0, 0.0);
                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(1, 0.0);
                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(2, 0.0);
                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(3, 0.0);

                    curve_fits_upper[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_upper[seed_eta_idx][loge_idx], envelope_plot_upper[seed_eta_idx][loge_idx], ("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "UPPER_FIT", fitStatus, fitSuccess);

                    int counter = 0;
                    for(float try_w00=-0.01; try_w00<=0.01; try_w00=try_w00+0.001){
                        for(float try_w01=-0.01; try_w01<=0.01; try_w01=try_w01+0.001){
                            for(float try_w10=-0.01; try_w10<=0.01; try_w10=try_w10+0.001){
                                for(float try_w11=-0.01; try_w11<=0.01; try_w11=try_w11+0.001){
                                    cout<<++counter<<endl;
                                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(0, try_w00);
                                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(1, try_w01);
                                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(2, try_w10);
                                    curve_fits_upper[seed_eta_idx][loge_idx] -> SetParameter(3, try_w11);

                                    curve_fits_upper[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_upper[seed_eta_idx][loge_idx], envelope_plot_upper[seed_eta_idx][loge_idx], ("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "UPPER_FIT", fitStatus, fitSuccess);
                                    if(fitSuccess) break;
                                }
                                if(fitSuccess) break;
                            }
                            if(fitSuccess) break;
                        }
                        if(fitSuccess) break;
                    }
                }

                curve_fits_lower[seed_eta_idx][loge_idx] = new TF1(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),fitFunc_lower,-0.6,0.6,4);
                curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameters(original_w00,original_w01,original_w10,original_w11,0.0,0.0);

                curve_fits_lower[seed_eta_idx][loge_idx] -> SetParName(0, "w00");
                curve_fits_lower[seed_eta_idx][loge_idx] -> SetParName(1, "w01");
                curve_fits_lower[seed_eta_idx][loge_idx] -> SetParName(2, "w10");
                curve_fits_lower[seed_eta_idx][loge_idx] -> SetParName(3, "w11");

                curve_fits_lower[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_lower[seed_eta_idx][loge_idx], envelope_plot_lower[seed_eta_idx][loge_idx], ("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "LOWER_FIT", fitStatus, fitSuccess);

                if(!fitSuccess){
                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(0, 0.0);
                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(1, 0.0);
                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(2, 0.0);
                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(3, 0.0);

                    curve_fits_lower[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_lower[seed_eta_idx][loge_idx], envelope_plot_lower[seed_eta_idx][loge_idx], ("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "LOWER_FIT", fitStatus, fitSuccess);
                    int counter = 0;
                    for(float try_w00=-0.01; try_w00<=0.01; try_w00=try_w00+0.001){
                        for(float try_w01=-0.01; try_w01<=0.01; try_w01=try_w01+0.001){
                            for(float try_w10=-0.01; try_w10<=0.01; try_w10=try_w10+0.001){
                                for(float try_w11=-0.01; try_w11<=0.01; try_w11=try_w11+0.001){
                                    cout<<++counter<<endl;
                                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(0, try_w00);
                                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(1, try_w01);
                                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(2, try_w10);
                                    curve_fits_lower[seed_eta_idx][loge_idx] -> SetParameter(3, try_w11);

                                    curve_fits_lower[seed_eta_idx][loge_idx] = fitMinimization(curve_fits_lower[seed_eta_idx][loge_idx], envelope_plot_lower[seed_eta_idx][loge_idx], ("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),seed_eta_idx,loge_idx, "LOWER_FIT", fitStatus, fitSuccess);
                                    if(fitSuccess) break;
                                }
                                if(fitSuccess) break;
                            }
                            if(fitSuccess) break;
                        }
                        if(fitSuccess) break;
                    }
                } 
            }

            //plotting
            gStyle->SetOptStat(0);
            legends[seed_eta_idx][loge_idx] = new TLegend(-0.55,-0.14,-0.3,-0.08,"","");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> SetTitle((titles_etas[seed_eta_idx] + "     " + titles_loge[loge_idx]).c_str());
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> SetTitle("dPhi");
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> SetTitle("dEta");

            //fitting points plots
            envelope_plot_upper[seed_eta_idx][loge_idx] -> SetMarkerStyle(34);  
            envelope_plot_upper[seed_eta_idx][loge_idx] -> SetMarkerColor(2);  //red
            envelope_plot_upper[seed_eta_idx][loge_idx] -> SetMarkerSize(3);

            envelope_plot_lower[seed_eta_idx][loge_idx] -> SetMarkerStyle(47);
            envelope_plot_lower[seed_eta_idx][loge_idx] -> SetMarkerColor(6);  //pink
            envelope_plot_lower[seed_eta_idx][loge_idx] -> SetMarkerSize(3);

            curve_fits_upper[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            curve_fits_lower[seed_eta_idx][loge_idx] -> SetLineWidth(4);

            curve_fits_lower[seed_eta_idx][loge_idx] -> SetLineColor(6);  //pink


            legends[seed_eta_idx][loge_idx] -> AddEntry(curve_fits_upper[seed_eta_idx][loge_idx], "Upper - Fit","l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(curve_fits_lower[seed_eta_idx][loge_idx], "Lower - Fit","l");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> Draw("COLZ");
            envelope_plot_upper[seed_eta_idx][loge_idx] -> Draw("P, sames");
            envelope_plot_lower[seed_eta_idx][loge_idx] -> Draw("P, sames");

            if(fit){

                curve_fits_upper[seed_eta_idx][loge_idx] -> Draw("SAME");
                curve_fits_lower[seed_eta_idx][loge_idx] -> Draw("SAME");
            }

            c1->SaveAs(("envelope/envelope_fit_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".png").c_str());
            c1->SaveAs(("envelope/envelope_fit_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".pdf").c_str());

            //aggregate parameter values
            if(fit){

                localParam_outfile << mid_etas[seed_eta_idx] << "\t" 
                                    << mid_loges[loge_idx] << "\t"
                                    << "UP" << "\t"
                                    << "FIT" << "\t"
                                    << curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(0) << "\t"
                                    << curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(1) << "\t"
                                    << curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(2) << "\t"
                                    << curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(3) << endl
                                    << mid_etas[seed_eta_idx] << "\t" 
                                    << mid_loges[loge_idx] << "\t"
                                    << "LOW" << "\t"
                                    << "FIT" << "\t"
                                    << curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(0) << "\t"
                                    << curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(1) << "\t"
                                    << curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(2) << "\t"
                                    << curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(3) << endl;

                vector<float> weight;
                weight.push_back(1); //0 UNWEIGHTED
                weight.push_back(total_clus); //1
                weight.push_back(total_clus_et); //2
                weight.push_back(total_clus_e); //3
                weight.push_back(mid_etas[seed_eta_idx]); //4

                for(int aveIdx = 0; aveIdx < weightNum; aveIdx++){
                    fit_w00_upper_weighted.at(aveIdx).push_back((curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(0)) * weight.at(aveIdx));
                    fit_w01_upper_weighted.at(aveIdx).push_back((curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(1)) * weight.at(aveIdx));
                    fit_w10_upper_weighted.at(aveIdx).push_back((curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(2)) * weight.at(aveIdx));
                    fit_w11_upper_weighted.at(aveIdx).push_back((curve_fits_upper[seed_eta_idx][loge_idx] -> GetParameter(3)) * weight.at(aveIdx));

                    fit_w00_lower_weighted.at(aveIdx).push_back((curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(0)) * weight.at(aveIdx));
                    fit_w01_lower_weighted.at(aveIdx).push_back((curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(1)) * weight.at(aveIdx));
                    fit_w10_lower_weighted.at(aveIdx).push_back((curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(2)) * weight.at(aveIdx));
                    fit_w11_lower_weighted.at(aveIdx).push_back((curve_fits_lower[seed_eta_idx][loge_idx] -> GetParameter(3)) * weight.at(aveIdx));

                    fit_w00_upper_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w01_upper_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w10_upper_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w11_upper_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w00_lower_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w01_lower_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w10_lower_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                    fit_w11_lower_weighted_denom.at(aveIdx).push_back(weight.at(aveIdx));
                }
            }
        } //loget index
    } //seed eta index
}

void dPhiOptimization(){
//optimization of dynamic dPhi window parameters
    TGraphErrors *dphi_envelope_plot[4];
    TF1 *dphi_fits[4];
    TF1 *dphi_original[4];

    TLegend* dphi_legends[4]; 

    TCanvas *c2 = new TCanvas("c2","c2",3600,2400);

    float xrange_min = -0.4,
          xrange_max = 1.5;

    //loop over dphi function areas
    for(int dphi_idx = 0; dphi_idx<4; dphi_idx++){

        int envelope_percent = 98,
            x_bins_total     = cluster_dPhi_vs_loget[dphi_idx] -> GetNbinsX(),
            y_bins_total     = cluster_dPhi_vs_loget[dphi_idx] -> GetNbinsY();

        int point_counter = 0;

        c2->Clear();
        
        vector<float> x_points, y_points, 
                      x_errors, y_errors;

        //loop over x bins
        for(int xbin_idx=1; xbin_idx<x_bins_total+1; xbin_idx++){
            float xbin_integral = cluster_dPhi_vs_loget[dphi_idx] -> Integral(xbin_idx, xbin_idx, 1, y_bins_total+1);

            float envelope_biny = 0, 
                  sum_temp = 0;

            float xbin_percentage = ((float)envelope_percent * xbin_integral) / 100.0;

            float xBin_val = cluster_dPhi_vs_loget[dphi_idx] -> GetXaxis() -> GetBinCenter(xbin_idx);

            //skip if the bin has no contents
            if(xbin_integral <= 0.0) continue;


            int n_entries = caloClusters_dPhiDist[dphi_idx][xbin_idx]->GetEntries();
            if(n_entries <= 1) continue;

            //loop over y bins
            for(int ybin_idx=1; ybin_idx<y_bins_total+1; ybin_idx++){
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

            vector<float> x_err = XErrors_DPhi(xbin_idx, envelope_biny, dphi_idx),
                          y_err = YErrors_DPhi(xbin_idx, envelope_biny, dphi_idx);

            
            x_points.push_back(envelope_x);
            x_errors.push_back(x_err[2]);

            y_points.push_back(envelope_y);
            y_errors.push_back(y_err[2]);

        } //x bin loop
        dphi_envelope_plot[dphi_idx] = new TGraphErrors(x_points.size(), &x_points[0], &y_points[0], &x_errors[0], &y_errors[0]);
            

        //fitting
        dphi_fits[dphi_idx] = new TF1(("dphiFit"+to_string(dphi_idx)).c_str(),fitFunc_dPhi,xrange_min,xrange_max,4);
        dphi_fits[dphi_idx] -> SetParameters(yoffset_orig[dphi_idx], scale_orig[dphi_idx], xoffset_orig[dphi_idx], width_orig[dphi_idx]);

        dphi_fits[dphi_idx] -> SetParName(0, "yoffset");
        dphi_fits[dphi_idx] -> SetParName(1, "scale");
        dphi_fits[dphi_idx] -> SetParName(2, "xoffset");
        dphi_fits[dphi_idx] -> SetParName(3, "width");


        string fitStatus = "";
        bool fitSuccess = false;
        fitMinimization(dphi_fits[dphi_idx], dphi_envelope_plot[dphi_idx], ("dphiFit"+to_string(dphi_idx)).c_str(),0,0, ("DPhi_FIT_"+to_string(dphi_idx)).c_str(), fitStatus, fitSuccess);

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

        c2 -> SaveAs(("envelope/dphi_fit_"+to_string(dphi_idx)+".png").c_str());

        yoffset_opt.push_back(dphi_fits[dphi_idx]->GetParameter(0));
        scale_opt.push_back(dphi_fits[dphi_idx]->GetParameter(1));
        xoffset_opt.push_back(dphi_fits[dphi_idx]->GetParameter(2));
        width_opt.push_back(dphi_fits[dphi_idx]->GetParameter(3));

    } //dphi index
}

void averageSeparationParams(){
//aggregate the values of the parabola separation fits
    aveParam_outfile << "AVERAGE_WEIGHT\tW00\tW01\tW10\tW11"<<endl;

    for(int aveIdx = 0; aveIdx < weightNum; aveIdx++){
        float w00_upper_weighted = (accumulate(fit_w00_upper_weighted.at(aveIdx).begin(),fit_w00_upper_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w00_upper_weighted_denom.at(aveIdx).begin(),fit_w00_upper_weighted_denom.at(aveIdx).end(),0.0)),
              w01_upper_weighted = (accumulate(fit_w01_upper_weighted.at(aveIdx).begin(),fit_w01_upper_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w01_upper_weighted_denom.at(aveIdx).begin(),fit_w01_upper_weighted_denom.at(aveIdx).end(),0.0)),
              w10_upper_weighted = (accumulate(fit_w10_upper_weighted.at(aveIdx).begin(),fit_w10_upper_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w10_upper_weighted_denom.at(aveIdx).begin(),fit_w10_upper_weighted_denom.at(aveIdx).end(),0.0)),
              w11_upper_weighted = (accumulate(fit_w11_upper_weighted.at(aveIdx).begin(),fit_w11_upper_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w11_upper_weighted_denom.at(aveIdx).begin(),fit_w11_upper_weighted_denom.at(aveIdx).end(),0.0)),
              w00_lower_weighted = (accumulate(fit_w00_lower_weighted.at(aveIdx).begin(),fit_w00_lower_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w00_lower_weighted_denom.at(aveIdx).begin(),fit_w00_lower_weighted_denom.at(aveIdx).end(),0.0)),
              w01_lower_weighted = (accumulate(fit_w01_lower_weighted.at(aveIdx).begin(),fit_w01_lower_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w01_lower_weighted_denom.at(aveIdx).begin(),fit_w01_lower_weighted_denom.at(aveIdx).end(),0.0)),
              w10_lower_weighted = (accumulate(fit_w10_lower_weighted.at(aveIdx).begin(),fit_w10_lower_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w10_lower_weighted_denom.at(aveIdx).begin(),fit_w10_lower_weighted_denom.at(aveIdx).end(),0.0)),
              w11_lower_weighted = (accumulate(fit_w11_lower_weighted.at(aveIdx).begin(),fit_w11_lower_weighted.at(aveIdx).end(),0.0)) / (accumulate(fit_w11_lower_weighted_denom.at(aveIdx).begin(),fit_w11_lower_weighted_denom.at(aveIdx).end(),0.0));

        float final_w00_weighted =   ( (accumulate(fit_w00_upper_weighted.at(aveIdx).begin(),fit_w00_upper_weighted.at(aveIdx).end(),0.0)) + (accumulate(fit_w00_lower_weighted.at(aveIdx).begin(),fit_w00_lower_weighted.at(aveIdx).end(),0.0)) ) 
                                   / ( (accumulate(fit_w00_upper_weighted_denom.at(aveIdx).begin(),fit_w00_upper_weighted_denom.at(aveIdx).end(),0.0)) + (accumulate(fit_w00_lower_weighted_denom.at(aveIdx).begin(),fit_w00_lower_weighted_denom.at(aveIdx).end(),0.0)) ),
              final_w01_weighted =   ( (accumulate(fit_w01_upper_weighted.at(aveIdx).begin(),fit_w01_upper_weighted.at(aveIdx).end(),0.0)) + (accumulate(fit_w01_lower_weighted.at(aveIdx).begin(),fit_w01_lower_weighted.at(aveIdx).end(),0.0)) ) 
                                   / ( (accumulate(fit_w01_upper_weighted_denom.at(aveIdx).begin(),fit_w01_upper_weighted_denom.at(aveIdx).end(),0.0)) + (accumulate(fit_w01_lower_weighted_denom.at(aveIdx).begin(),fit_w01_lower_weighted_denom.at(aveIdx).end(),0.0)) ),
              final_w10_weighted =   ( (accumulate(fit_w10_upper_weighted.at(aveIdx).begin(),fit_w10_upper_weighted.at(aveIdx).end(),0.0)) + (accumulate(fit_w10_lower_weighted.at(aveIdx).begin(),fit_w10_lower_weighted.at(aveIdx).end(),0.0)) ) 
                                   / ( (accumulate(fit_w10_upper_weighted_denom.at(aveIdx).begin(),fit_w10_upper_weighted_denom.at(aveIdx).end(),0.0)) + (accumulate(fit_w10_lower_weighted_denom.at(aveIdx).begin(),fit_w10_lower_weighted_denom.at(aveIdx).end(),0.0)) ),
              final_w11_weighted =   ( (accumulate(fit_w11_upper_weighted.at(aveIdx).begin(),fit_w11_upper_weighted.at(aveIdx).end(),0.0)) + (accumulate(fit_w11_lower_weighted.at(aveIdx).begin(),fit_w11_lower_weighted.at(aveIdx).end(),0.0)) ) 
                                   / ( (accumulate(fit_w11_upper_weighted_denom.at(aveIdx).begin(),fit_w11_upper_weighted_denom.at(aveIdx).end(),0.0)) + (accumulate(fit_w11_lower_weighted_denom.at(aveIdx).begin(),fit_w11_lower_weighted_denom.at(aveIdx).end(),0.0)) );

        aveParam_outfile << weightTitles[aveIdx] << final_w00_weighted << final_w01_weighted << final_w10_weighted << final_w11_weighted << endl;

    }
}

void final_params_dPhi(){
//aggregate the values of the dphi fits

    dPhiParam_outfile << "ETA_BIN\tYOFFSET\tSCALE\tXOFFSET\tWIDTH"<<endl;

    string eta_bins[] = {"EB", "EE_0", "EE_1", "EE_2"};

    for(int i = 0; i < 4; i++){
        dPhiParam_outfile << eta_bins[i] << "\t"
                          << yoffset_opt[i] << "\t"
                          << scale_opt[i] << "\t"
                          << xoffset_opt[i] << "\t"
                          << width_opt[i] << endl;
    }

}

void ReadInfile(string inputFile){
    TFile *hist_infile = TFile::Open(inputFile.c_str());

    double seedEtaStep = (maxSeedEta - minSeedEta) / seedEtaBins;
    double logEStep = (maxLogE - minLogE) / logEBins;

    filename_etas.resize(seedEtaBins);
    filename_loges.resize(logEBins);

    double seedEtaVal = minSeedEta;
    for(int seedEtaIdx = 0; seedEtaIdx < seedEtaBins; seedEtaIdx++){
        titles_etas.push_back((to_string(seedEtaVal) + " < |#eta| #leq " + to_string(seedEtaVal + seedEtaStep)).c_str());
        mid_etas.push_back((2*seedEtaVal) / 2.0);
        seedEtaVal+=seedEtaStep;

        std::ostringstream ss;
        ss << std::setw(3) << std::setfill('0') << to_string(seedEtaIdx);
        filename_etas[seedEtaIdx] = ss.str();
    }

    double logEVal = minLogE;
    for(int logEIdx = 0; logEIdx < logEBins; logEIdx++){
        titles_loge.push_back((to_string(logEVal) + " #leq log_{10}(E) < " + to_string(logEVal + logEStep)).c_str());
        mid_loges.push_back((2*logEVal) / 2.0);
        logEVal+=logEStep;

        std::ostringstream ss;
        ss << std::setw(3) << std::setfill('0') << to_string(logEIdx);
        filename_loges[logEIdx] = ss.str();
    }

    fitPlot.resize(seedEtaBins, vector<bool>(logEBins, false));
    caloClusters_shape_eBins.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_shape_eBins_etWeight.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_shape_eBins_eWeight.resize(seedEtaBins, vector<TH2F*>(logEBins));
    caloClusters_eBins_dEtaDist.resize(seedEtaBins, vector<vector<TH1F*>>(logEBins, vector<TH1F*>(caloClusterShapeDEtaDistBins)));

    cluster_dPhi_vs_loget.resize(dPhiWindowEtaBins);
    caloClusters_dPhiDist.resize(dPhiWindowEtaBins, vector<TH1F*>(dPhiWindowETDistBins));

    fit_w00_up.resize(seedEtaBins, vector<float>(logEBins));
    fit_w01_up.resize(seedEtaBins, vector<float>(logEBins));
    fit_w10_up.resize(seedEtaBins, vector<float>(logEBins));
    fit_w11_up.resize(seedEtaBins, vector<float>(logEBins));
    fit_w00_low.resize(seedEtaBins, vector<float>(logEBins));
    fit_w01_low.resize(seedEtaBins, vector<float>(logEBins));
    fit_w10_low.resize(seedEtaBins, vector<float>(logEBins));
    fit_w11_low.resize(seedEtaBins, vector<float>(logEBins));

    for(int seedEtaIdx = 0; seedEtaIdx < seedEtaBins; seedEtaIdx++){
        for(int logEIdx = 0; logEIdx < logEBins; logEIdx++){
            hist_infile->GetObject(("caloClusters_shape_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),caloClusters_shape_eBins[seedEtaIdx][logEIdx]);
            hist_infile->GetObject(("caloClusters_shape_eBins_etWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),caloClusters_shape_eBins_etWeight[seedEtaIdx][logEIdx]);
            hist_infile->GetObject(("caloClusters_shape_eBins_eWeight_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)).c_str(),caloClusters_shape_eBins_eWeight[seedEtaIdx][logEIdx]);
            
            for(int dEtaDistBins = 0; dEtaDistBins < caloClusterShapeDEtaDistBins; dEtaDistBins++){
                hist_infile->GetObject(("ET_vs_dEta_eBins_"+to_string(seedEtaIdx)+"_"+to_string(logEIdx)+"_"+to_string(dEtaDistBins)).c_str(),caloClusters_eBins_dEtaDist[seedEtaIdx][logEIdx][dEtaDistBins]);
            }
        }
    }
    for(int ii = 0; ii < dPhiWindowEtaBins; ii++){
        hist_infile->GetObject(("caloClusters_dPhi_vs_logET_etaBin_"+to_string(ii)).c_str(),cluster_dPhi_vs_loget[ii]);

        for(int x = 0; x < dPhiWindowETDistBins; x++){
            hist_infile->GetObject(("ET_vs_dPhi_"+to_string(ii)+"_"+to_string(x)).c_str(),caloClusters_dPhiDist[ii][x]);
        }
    }

    curr_logET = 0.0;
    curr_seedEta = 0.0;
}

void EnvelopeOptimizer(string inputFile, string outputFile, string localParamOutFile, string aveParamOutFile, string dPhiParamOutFile, int errorType, bool ringRejection, bool binRejection){
//main program

    fit_outfile.open("Output/fit_status.txt", std::ofstream::out | std::ofstream::trunc);
    aveParam_outfile.open(aveParamOutFile, std::ofstream::out | std::ofstream::trunc);
    localParam_outfile.open(localParamOutFile, std::ofstream::out | std::ofstream::trunc);
    dPhiParam_outfile.open(dPhiParamOutFile, std::ofstream::out | std::ofstream::trunc);
    
    ReadInfile(inputFile);
    separationOptimization(1,errorType,ringRejection,binRejection);
                  //(bool fit, int error=0, bool ringRejection=0, bool binRejection=0)
    dPhiOptimization();
    averageSeparationParams();
    final_params_dPhi();
    //fullParabola_opt(1,2,1,1,1);
    //dPhi_opt();
    //final_params_dPhi();
    //final_params();
    //curveAverage();
    //plotRegionalCurves();
    //parabolaCompareHeatMap_Individual();
    //printParams();

    //savePlot();
 
    fit_outfile.close();
    aveParam_outfile.close();
    localParam_outfile.close();
    dPhiParam_outfile.close();
}