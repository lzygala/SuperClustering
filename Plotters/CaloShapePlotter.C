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

#include<sys/stat.h>
#include<errno.h>

    int seedEtaBins = 15,
        logEBins = 15,
        caloClusterShapeDEtaDistBins = 120,
        dPhiWindowEtaBins = 4,
        dPhiWindowETDistBins = 200;

    double minSeedEta = 0.0, 
           maxSeedEta = 3.0, 
           minLogE = -1.0, 
           maxLogE = 2.0;
    
    double maximum_val = 0.0174;
    double minimum_val = -0.0087;

    vector<double> midEtas,     //[seedEtaBins]
                   midLogEs;    //[logEBins]

    vector<TH2F*> cluster_dPhi_vs_loget;            //[dPhiWindowEtaBins]
    vector<vector<TH1F*>> caloClusters_dPhiDist;    //[dPhiWindowEtaBins][dPhiWindowETDistBins]

    vector<vector<TH2F*>> caloClusters_shape_eBins_etWeight;    //[seedEtaBins][logEBins]
    vector<vector<TH2F*>> caloClusters_shape_eBins_eWeight;     //[seedEtaBins][logEBins]
    vector<vector<TH2F*>> caloClusters_shape_eBins;             //[seedEtaBins][logEBins]
    vector<vector<vector<TH1F*>>> caloClusters_eBins_dEtaDist;  //[seedEtaBins][logEBins][caloClusterShapeDEtaDistBins]

    vector<vector<bool>> fitPlot;   //[seedEtaBins][logEBins]

    float curr_logE;
    float curr_seedEta;

    vector<string> averageTitles;

    vector<float> new_w00_averaged;
    vector<float> new_w01_averaged;
    vector<float> new_w10_averaged;
    vector<float> new_w11_averaged;

    vector<vector<vector<float>>> new_upperParams_local; //[seedEtaBins][logEBins][4]
    vector<vector<vector<float>>> new_lowerParams_local; //[seedEtaBins][logEBins][4]

    vector<string> dPhiSections;
    vector<float> new_yoffset,
                  new_scale,
                  new_xoffset,
                  new_width;

    //current values
    float phase1_w00 = -0.00571429;
    float phase1_w01 = -0.002;
    float phase1_w10 = 0.0135714;
    float phase1_w11 = 0.001;

    //original dPhi values
    vector<float> phase1_yoffset = {0.07151, 0.05058, -0.0913, -0.6246},
                  phase1_scale   = {0.5656,  0.7131,  44.04,   13.17},
                  phase1_xoffset = {0.2931,  0.01668, -5.326,  -7.037},
                  phase1_width   = {0.2976,  0.4114,  1.184,   2.836};

    vector<string> titles_etas, titles_loge, filename_etas, filename_loges;

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

Double_t upper_parabola(Double_t *x, Double_t *par){

    float eta_max = par[0];
    float logET = par[1];

    float w00 = par[2];
    float w01 = par[3];
    float w10 = par[4];
    float w11 = par[5];

    float p00 = -0.107537;
    float p01 = 0.590969;
    float p02 = -0.076494;
    float p10 = -0.0268843;
    float p11 = 0.147742;
    float p12 = -0.0191235;
        
    float c_upper = (p00 * pow(eta_max*sin(eta_max), 2)) + (p01 * eta_max*sin(eta_max)) + p02;
    float d_upper =( w10*eta_max*sin(eta_max) ) + (w11 / sqrt(1.1 + logET));
    float d_lower =( w00*eta_max*sin(eta_max) ) + (w01 / sqrt(1.1 + logET));
    float b_upper = d_upper - 0.5*(d_lower + d_upper);
    float a_upper = ((1 / (4 * c_upper))) - fabs(b_upper);
    Double_t upper_curve = (std::max((1 / (4 * a_upper)),0.0f)*pow(x[0],2)) + std::max(b_upper, 0.0087f) + 0.0087;

    return upper_curve;
}

Double_t lower_parabola(Double_t *x, Double_t *par){

    float eta_max = par[0];
    float logET = par[1];

    float w00 = par[2];
    float w01 = par[3];
    float w10 = par[4];
    float w11 = par[5];

    float p00 = -0.107537;
    float p01 = 0.590969;
    float p02 = -0.076494;
    float p10 = -0.0268843;
    float p11 = 0.147742;
    float p12 = -0.0191235;

    float c_lower = (p10 * pow(eta_max * sin(eta_max), 2)) + (p11 * eta_max * sin(eta_max)) + p12;
    float d_upper = (w10 * eta_max * sin(eta_max)) + (w11 / sqrt(1.1 + logET));
    float d_lower = (w00 * eta_max * sin(eta_max)) + (w01 / sqrt(1.1 + logET));
    float b_lower = d_lower - 0.5 * (d_lower + d_upper);
    float a_lower = ((1 / (4 * c_lower))) - fabs(b_lower);
    float lower_curve = (std::max((1 / (4 * a_lower)), 0.0f) * pow(x[0], 2)) + std::min(b_lower, -0.0087f);

    return lower_curve;
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

float dPhi_window(float yoffset, float scale, float xoffset, float width, float cutoff, float saturation, float logET){

    width = 1.0 / width;

    float maxdphi = yoffset + scale / (1 + std::exp((logET - xoffset) * width));
    maxdphi = std::min(maxdphi, cutoff);
    maxdphi = std::max(maxdphi, saturation);

    return maxdphi;
}

void plotRegionalCurves(int aveIdx = 0){
//plot Phase 1 curves, local curves + parameters, optimized curves + parameters
    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);
    TLegend* legends[seedEtaBins][logEBins]; 

    TF1 *upperParabolaAveraged[seedEtaBins][logEBins],      *lowerParabolaAveraged[seedEtaBins][logEBins],
        *upperParabolaLocal[seedEtaBins][logEBins],         *lowerParabolaLocal[seedEtaBins][logEBins], 
        *upperParabolaPhase1[seedEtaBins][logEBins],        *lowerParabolaPhase1[seedEtaBins][logEBins];

    TLine *minimum[seedEtaBins][logEBins], *maximum[seedEtaBins][logEBins];

    for(int seed_eta_idx = 0; seed_eta_idx < seedEtaBins; seed_eta_idx++){
        for(int loge_idx = 0; loge_idx < logEBins; loge_idx++){ 
            if(fitPlot[seed_eta_idx][loge_idx] == false) continue;

            minimum[seed_eta_idx][loge_idx] = new TLine(-0.1, minimum_val, 0.1, minimum_val);
            maximum[seed_eta_idx][loge_idx] = new TLine(-0.1, maximum_val, 0.1, maximum_val);

            upperParabolaLocal[seed_eta_idx][loge_idx] = new TF1(("upperParabolaLocal"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),upper_parabola,-0.6,0.6,6);
            upperParabolaLocal[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx], 
                                                                         new_upperParams_local[seed_eta_idx][loge_idx][0], //w00
                                                                         new_upperParams_local[seed_eta_idx][loge_idx][1], //w01
                                                                         new_upperParams_local[seed_eta_idx][loge_idx][2], //w10
                                                                         new_upperParams_local[seed_eta_idx][loge_idx][3]); //w11

            lowerParabolaLocal[seed_eta_idx][loge_idx] = new TF1(("lowerParabolaLocal"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),lower_parabola,-0.6,0.6,6);
            lowerParabolaLocal[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx],
                                                                         new_lowerParams_local[seed_eta_idx][loge_idx][0], //w00
                                                                         new_lowerParams_local[seed_eta_idx][loge_idx][1], //w01
                                                                         new_lowerParams_local[seed_eta_idx][loge_idx][2], //w10
                                                                         new_lowerParams_local[seed_eta_idx][loge_idx][3]); //w11

            upperParabolaAveraged[seed_eta_idx][loge_idx] = new TF1(("upperParabolaAveraged"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),upper_parabola,-0.6,0.6,6);
            upperParabolaAveraged[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx],
                                                                             new_w00_averaged.at(aveIdx), 
                                                                             new_w01_averaged.at(aveIdx), 
                                                                             new_w10_averaged.at(aveIdx), 
                                                                             new_w11_averaged.at(aveIdx));

            lowerParabolaAveraged[seed_eta_idx][loge_idx] = new TF1(("lowerParabolaAveraged"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),lower_parabola,-0.6,0.6,6);
            lowerParabolaAveraged[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx],
                                                                             new_w00_averaged.at(aveIdx), 
                                                                             new_w01_averaged.at(aveIdx), 
                                                                             new_w10_averaged.at(aveIdx), 
                                                                             new_w11_averaged.at(aveIdx));

            upperParabolaPhase1[seed_eta_idx][loge_idx] = new TF1(("upperParabolaPhase1"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),upper_parabola,-0.6,0.6,6);
            upperParabolaPhase1[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx],
                                                                          phase1_w00, 
                                                                          phase1_w01, 
                                                                          phase1_w10, 
                                                                          phase1_w11);

            lowerParabolaPhase1[seed_eta_idx][loge_idx] = new TF1(("lowerParabolaPhase1"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),lower_parabola,-0.6,0.6,6);
            lowerParabolaPhase1[seed_eta_idx][loge_idx] -> SetParameters(midEtas[seed_eta_idx], midLogEs[loge_idx], 
                                                                          phase1_w00, 
                                                                          phase1_w01, 
                                                                          phase1_w10, 
                                                                          phase1_w11);

            gStyle->SetOptStat(0);
            legends[seed_eta_idx][loge_idx] = new TLegend(-0.55,-0.14,-0.3,-0.05,"","");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> SetTitle((titles_etas[seed_eta_idx] + "     " + titles_loge[loge_idx]).c_str());
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> SetTitle("dPhi");
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> SetTitle("dEta");

            upperParabolaLocal[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            lowerParabolaLocal[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            upperParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            lowerParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            upperParabolaPhase1[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            lowerParabolaPhase1[seed_eta_idx][loge_idx] -> SetLineWidth(4);

            upperParabolaLocal[seed_eta_idx][loge_idx] -> SetLineColor(kRed);
            lowerParabolaLocal[seed_eta_idx][loge_idx] -> SetLineColor(kRed);
            upperParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineColor(kBlue);
            lowerParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineColor(kBlue);
            upperParabolaPhase1[seed_eta_idx][loge_idx] -> SetLineColor(kGreen);
            lowerParabolaPhase1[seed_eta_idx][loge_idx] -> SetLineColor(kGreen);

            upperParabolaLocal[seed_eta_idx][loge_idx] -> SetLineStyle(10);
            lowerParabolaLocal[seed_eta_idx][loge_idx] -> SetLineStyle(10);
            upperParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineStyle(9);
            lowerParabolaAveraged[seed_eta_idx][loge_idx] -> SetLineStyle(9);

            minimum[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            maximum[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            minimum[seed_eta_idx][loge_idx] -> SetLineColor(kBlack);
            maximum[seed_eta_idx][loge_idx] -> SetLineColor(kBlack);


            legends[seed_eta_idx][loge_idx] -> AddEntry(upperParabolaLocal[seed_eta_idx][loge_idx], "Local Region Parameters","l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(upperParabolaAveraged[seed_eta_idx][loge_idx], (averageTitles[aveIdx]+"Averaged Parameters").c_str(),"l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(upperParabolaPhase1[seed_eta_idx][loge_idx], "Phase 1 Parameters","l");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> Draw("COLZ");
            legends[seed_eta_idx][loge_idx] -> Draw("SAME");
            upperParabolaPhase1[seed_eta_idx][loge_idx] -> Draw("SAME");
            lowerParabolaPhase1[seed_eta_idx][loge_idx] -> Draw("SAME");
            upperParabolaAveraged[seed_eta_idx][loge_idx] -> Draw("SAME");
            lowerParabolaAveraged[seed_eta_idx][loge_idx] -> Draw("SAME");
            lowerParabolaLocal[seed_eta_idx][loge_idx] -> Draw("SAME");
            upperParabolaLocal[seed_eta_idx][loge_idx] -> Draw("SAME");
            minimum[seed_eta_idx][loge_idx] -> Draw("SAME");
            maximum[seed_eta_idx][loge_idx] -> Draw("SAME");

            string outdir = "Output/Plots/localRegionCurves/";
            setOutput(outdir);

            c1->SaveAs((outdir+"envFitCompare_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".png").c_str());
            c1->SaveAs((outdir+"envFitCompare_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".pdf").c_str());
        }
    }

}

void parabolaCompareHeatMap(){
//Compare local parabola values against various averaged parabola values
    int weightNum = new_w00_averaged.size();

    TH2F *parCom_LocalVFinal_Up[weightNum], *parCom_LocalVFinal_Low[weightNum];
    TGraphErrors *g_minPoints_Up[weightNum], *g_minPoints_Low[weightNum];

    TCanvas *c1 = new TCanvas("c1","c1",3600,3600);

    for(int aveIdx = 0; aveIdx < weightNum; aveIdx++){

        parCom_LocalVFinal_Up[aveIdx] = new TH2F(("parCom_LocalVFinal_Up_"+to_string(aveIdx)).c_str(), ("parCom_LocalVFinal_Up_"+to_string(aveIdx)).c_str(),15,-1.0,2.0,15,0.0,3.0);
        parCom_LocalVFinal_Low[aveIdx] = new TH2F(("parCom_LocalVFinal_Low_"+to_string(aveIdx)).c_str(), ("parCom_LocalVFinal_Low_"+to_string(aveIdx)).c_str(),15,-1.0,2.0,15,0.0,3.0);
        parCom_LocalVFinal_Width[aveIdx] = new TH2F(("parCom_LocalVFinal_Width_"+to_string(aveIdx)).c_str(), ("parCom_LocalVFinal_Width_"+to_string(aveIdx)).c_str(),15,-1.0,2.0,15,0.0,3.0);
        
        vector<float> minPoints_Up_x, minPoints_Up_y, minPoints_Low_x, minPoints_Low_y;

        for(int seed_eta_idx = 0; seed_eta_idx < seedEtaBins; seed_eta_idx++){ 
            
            for(int loge_idx = 0; loge_idx < logEBins; loge_idx++){ 
                if(fitPlot[seed_eta_idx][loge_idx] == false) continue;

                float local_up_dEta = upper_parabola_dEtaValue(new_upperParams_local[seed_eta_idx][loge_idx][0],
                                                               new_upperParams_local[seed_eta_idx][loge_idx][1],
                                                               new_upperParams_local[seed_eta_idx][loge_idx][2],
                                                               new_upperParams_local[seed_eta_idx][loge_idx][3],
                                                               midLogEs[loge_idx],
                                                               midEtas[seed_eta_idx],
                                                               0);//xVal) 
                float local_low_dEta = lower_parabola_dEtaValue(new_lowerParams_local[seed_eta_idx][loge_idx][0],
                                                                new_lowerParams_local[seed_eta_idx][loge_idx][1],
                                                                new_lowerParams_local[seed_eta_idx][loge_idx][2],
                                                                new_lowerParams_local[seed_eta_idx][loge_idx][3],
                                                                midLogEs[loge_idx],
                                                                midEtas[seed_eta_idx],
                                                                0);//xVal)   

                float averaged_up_dEta = upper_parabola_dEtaValue(new_w00_averaged.at(aveIdx),
                                                                  new_w01_averaged.at(aveIdx),
                                                                  new_w10_averaged.at(aveIdx),
                                                                  new_w11_averaged.at(aveIdx),
                                                                  midLogEs[loge_idx],
                                                                  midEtas[seed_eta_idx],
                                                                  0);//xVal) );
                float averaged_low_dEta = lower_parabola_dEtaValue(new_w00_averaged.at(aveIdx),
                                                                   new_w01_averaged.at(aveIdx),
                                                                   new_w10_averaged.at(aveIdx),
                                                                   new_w11_averaged.at(aveIdx),
                                                                   midLogEs[loge_idx],
                                                                   midEtas[seed_eta_idx],
                                                                   0);//xVal) ); 

                float plotVal_up = local_up_dEta / averaged_up_dEta;
                float plotVal_low = local_low_dEta / averaged_low_dEta;
                float plotVal_width = (local_up_dEta - local_low_dEta) / (averaged_up_dEta - averaged_low_dEta);

                parCom_LocalVFinal_Up[aveIdx] -> Fill(midLogEs[loge_idx], midEtas[seed_eta_idx], plotVal_up);
                parCom_LocalVFinal_Low[aveIdx] -> Fill(midLogEs[loge_idx], midEtas[seed_eta_idx], plotVal_low);
                parCom_LocalVFinal_Width[aveIdx] -> Fill(midLogEs[loge_idx], midEtas[seed_eta_idx], plotVal_width);

                float epsilon = 0.000000001;
                if(fabs(local_up_dEta - 0.0174) < epsilon && fabs(averaged_up_dEta - 0.0174) < epsilon && fabs(local_up_dEta - averaged_up_dEta) < epsilon){
                    
                    minPoints_Up_x.push_back(midLogEs[loge_idx]);
                    minPoints_Up_y.push_back(midEtas[seed_eta_idx]);
                }
                if(fabs(local_low_dEta - minimum_val) < epsilon && fabs(averaged_low_dEta - minimum_val) < epsilon && fabs(local_low_dEta - averaged_low_dEta) < epsilon){
                    minPoints_Low_x.push_back(midLogEs[loge_idx]);
                    minPoints_Low_y.push_back(midEtas[seed_eta_idx]);
                }
                
            }
        }
        
        g_minPoints_Up[aveIdx] = new TGraphErrors(minPoints_Up_x.size(), &minPoints_Up_x[0], &minPoints_Up_y[0]);
        g_minPoints_Low[aveIdx] = new TGraphErrors(minPoints_Low_x.size(), &minPoints_Low_x[0], &minPoints_Low_y[0]);

        g_minPoints_Up[aveIdx] -> SetMarkerStyle(29);
        g_minPoints_Low[aveIdx] -> SetMarkerStyle(29);

        g_minPoints_Up[aveIdx] -> SetMarkerSize(5);
        g_minPoints_Low[aveIdx] -> SetMarkerSize(5);

        parCom_LocalVFinal_Up[aveIdx] -> SetTitle("Upper Parabolas: Local Parameter Curve / Averaged Parameter Curve");
        parCom_LocalVFinal_Up[aveIdx] -> GetXaxis() -> SetTitle("log(E)");
        parCom_LocalVFinal_Up[aveIdx] -> GetYaxis() -> SetTitle("Seed Eta");
        parCom_LocalVFinal_Up[aveIdx] -> GetZaxis() -> SetRangeUser(0.0,5.0);
        
        parCom_LocalVFinal_Up[aveIdx] -> GetXaxis() -> SetNdivisions(15,0,0);
        parCom_LocalVFinal_Up[aveIdx] -> GetYaxis() -> SetNdivisions(15,0,0);

        parCom_LocalVFinal_Up[aveIdx] -> GetXaxis() -> SetLabelSize(0.02);
        parCom_LocalVFinal_Up[aveIdx] -> GetYaxis() -> SetLabelSize(0.02);

        parCom_LocalVFinal_Low[aveIdx] -> SetTitle("Lower Parabolas: Local Parameter Curve / Averaged Parameter Curve");
        parCom_LocalVFinal_Low[aveIdx] -> GetXaxis() -> SetTitle("log(E)");
        parCom_LocalVFinal_Low[aveIdx] -> GetYaxis() -> SetTitle("Seed Eta");
        parCom_LocalVFinal_Low[aveIdx] -> GetZaxis() -> SetRangeUser(0.0,5.0);
        
        parCom_LocalVFinal_Low[aveIdx] -> GetXaxis() -> SetNdivisions(15,0,0);
        parCom_LocalVFinal_Low[aveIdx] -> GetYaxis() -> SetNdivisions(15,0,0);
        parCom_LocalVFinal_Low[aveIdx] -> GetXaxis() -> SetLabelSize(0.02);
        parCom_LocalVFinal_Low[aveIdx] -> GetYaxis() -> SetLabelSize(0.02);


        parCom_LocalVFinal_Width[aveIdx] -> SetTitle("Parabola Widths: Local Parameter Curve / Averaged Parameter Curve");
        parCom_LocalVFinal_Width[aveIdx] -> GetXaxis() -> SetTitle("log(E)");
        parCom_LocalVFinal_Width[aveIdx] -> GetYaxis() -> SetTitle("Seed Eta");
        parCom_LocalVFinal_Width[aveIdx] -> GetZaxis() -> SetTitle("Local Parameters Width / Averaged Parameters Width");
        parCom_LocalVFinal_Width[aveIdx] -> GetZaxis() -> SetRangeUser(0.0,5.0);
        
        parCom_LocalVFinal_Width[aveIdx] -> GetXaxis() -> SetNdivisions(15,0,0);
        parCom_LocalVFinal_Width[aveIdx] -> GetYaxis() -> SetNdivisions(15,0,0);
        parCom_LocalVFinal_Width[aveIdx] -> GetXaxis() -> SetLabelSize(0.02);
        parCom_LocalVFinal_Width[aveIdx] -> GetYaxis() -> SetLabelSize(0.02);

        string outdir = "Output/Plots/";
        setOutput(outdir);

        c1->SetGrid();
        parCom_LocalVFinal_Up[aveIdx] -> Draw("COLZ");
        g_minPoints_Up[aveIdx] -> Draw("P, sames");
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Up_"+averageTitles[aveIdx]+".png").c_str());
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Up_"+averageTitles[aveIdx]+".pdf").c_str());


        c1->SetGrid();
        parCom_LocalVFinal_Low[aveIdx] -> Draw("COLZ");
        g_minPoints_Low[aveIdx] -> Draw("P, sames");
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Low_"+averageTitles[aveIdx]+".png").c_str());
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Low_"+averageTitles[aveIdx]+".pdf").c_str());

        c1->SetGrid();
        parCom_LocalVFinal_Width[aveIdx] -> Draw("COLZ");
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Width_"+averageTitles[aveIdx]+".png").c_str());
        c1 -> SaveAs((outdir+"parCom_LocalVFinal_Width_"+averageTitles[aveIdx]+".pdf").c_str());
    }

}

void Plot(int aveIdx = 0){

    TF1 *fit_curves_upper[seedEtaBins][logEBins],
        *fit_curves_lower[seedEtaBins][logEBins],
        *phase1_curves_upper[seedEtaBins][logEBins],
        *phase1_curves_lower[seedEtaBins][logEBins];

    TLine *fit_dPhi_left[seedEtaBins][logEBins],
          *fit_dPhi_right[seedEtaBins][logEBins],
          *phase1_dPhi_left[seedEtaBins][logEBins],
          *phase1_dPhi_right[seedEtaBins][logEBins];

    TLegend *legends[seedEtaBins][logEBins], *legends_dPhi[seedEtaBins][logEBins]; 

    TCanvas *c1 = new TCanvas("c1","c1",3600,2400);

    for(int seed_eta_idx = 0; seed_eta_idx < seedEtaBins; seed_eta_idx++){ 
        for(int loge_idx = 0; loge_idx < logEBins; loge_idx++){ 

            float absSeedEta = midEtas.at(seed_eta_idx);
            int etaBin = ((int)(absSeedEta >= 1.479) + (int)(absSeedEta >= 1.75) + (int)(absSeedEta >= 2.0));
            float saturation, cutoff;
            switch (etaBin) {
                case 0:  // EB
                    saturation = 0.14;
                    cutoff = 0.60;
                    break;
                case 1:  // 1.479 -> 1.75
                    saturation = 0.14;
                    cutoff = 0.55;
                    break;
                case 2:  // 1.75 -> 2.0
                    saturation = 0.12;
                    cutoff = 0.45;
                    break;
                case 3:  // 2.0 and up
                    saturation = 0.12;
                    cutoff = 0.30;
                    break;
            }
            float curr_logE = midLogEs.at(loge_idx);
            float curr_seedEta = midEtas.at(seed_eta_idx);

            float fit_maxdPhi = dPhi_window(new_yoffset[etaBin], new_scale[etaBin], new_xoffset[etaBin], new_width[etaBin], cutoff, saturation, curr_logE);
            float phase1_maxdPhi = dPhi_window(phase1_yoffset[etaBin], phase1_scale[etaBin], phase1_xoffset[etaBin], phase1_width[etaBin], cutoff, saturation, curr_logE);


            fit_dPhi_left[seed_eta_idx][loge_idx] = new TLine(-(fit_maxdPhi), -0.15, -(fit_maxdPhi), 0.2);
            fit_dPhi_right[seed_eta_idx][loge_idx] = new TLine((fit_maxdPhi), -0.15, (fit_maxdPhi), 0.2);
            phase1_dPhi_left[seed_eta_idx][loge_idx] = new TLine(-(phase1_maxdPhi), -0.15, -(phase1_maxdPhi), 0.2);
            phase1_dPhi_right[seed_eta_idx][loge_idx] = new TLine((phase1_maxdPhi), -0.15, (phase1_maxdPhi), 0.2);

            fit_curves_upper[seed_eta_idx][loge_idx] = new TF1(("parabolaFit_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),upper_parabola,-0.6,0.6,6);
            fit_curves_upper[seed_eta_idx][loge_idx] -> SetParameters(curr_seedEta, curr_logE, new_w00_averaged[aveIdx], new_w01_averaged[aveIdx], new_w10_averaged[aveIdx], new_w11_averaged[aveIdx]);
            fit_curves_lower[seed_eta_idx][loge_idx] = new TF1(("parabolaFit_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),lower_parabola,-0.6,0.6,6);
            fit_curves_lower[seed_eta_idx][loge_idx] -> SetParameters(curr_seedEta, curr_logE, new_w00_averaged[aveIdx], new_w01_averaged[aveIdx], new_w10_averaged[aveIdx], new_w11_averaged[aveIdx]);

            phase1_curves_upper[seed_eta_idx][loge_idx] = new TF1(("parabolaphase1_upper_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),upper_parabola,-0.6,0.6,6);
            phase1_curves_upper[seed_eta_idx][loge_idx] -> SetParameters(curr_seedEta, curr_logE, phase1_w00, phase1_w01, phase1_w10, phase1_w11);
            phase1_curves_lower[seed_eta_idx][loge_idx] = new TF1(("parabolaphase1_lower_"+to_string(seed_eta_idx)+"_"+to_string(loge_idx)).c_str(),lower_parabola,-0.6,0.6,6);
            phase1_curves_lower[seed_eta_idx][loge_idx] -> SetParameters(curr_seedEta, curr_logE, phase1_w00, phase1_w01, phase1_w10, phase1_w11);


            //plotting
            gStyle->SetOptStat(0);
            legends[seed_eta_idx][loge_idx] = new TLegend(-0.55,-0.14,-0.3,-0.08,"","");
            legends_dPhi[seed_eta_idx][loge_idx] = new TLegend(0.3,-0.14,0.55,-0.08,"","");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> SetTitle((titles_etas[seed_eta_idx] + "     " + titles_loge[loge_idx]).c_str());
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetXaxis() -> SetTitle("dPhi");
            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> GetYaxis() -> SetTitle("dEta");

            fit_curves_upper[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            fit_curves_lower[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            phase1_curves_upper[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            phase1_curves_lower[seed_eta_idx][loge_idx] -> SetLineWidth(4);

            fit_dPhi_left[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            fit_dPhi_right[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            phase1_dPhi_left[seed_eta_idx][loge_idx] -> SetLineWidth(4);
            phase1_dPhi_right[seed_eta_idx][loge_idx] -> SetLineWidth(4);

            fit_curves_lower[seed_eta_idx][loge_idx] -> SetLineColor(6);  //pink
            phase1_curves_upper[seed_eta_idx][loge_idx] -> SetLineColor(kGreen);  //green
            phase1_curves_lower[seed_eta_idx][loge_idx] -> SetLineColor(kCyan);  //cyan

            fit_dPhi_left[seed_eta_idx][loge_idx] -> SetLineColor(kBlack);
            fit_dPhi_right[seed_eta_idx][loge_idx] -> SetLineColor(kBlack);
            phase1_dPhi_left[seed_eta_idx][loge_idx] -> SetLineColor(kBlue);
            phase1_dPhi_right[seed_eta_idx][loge_idx] -> SetLineColor(kBlue);

            phase1_curves_upper[seed_eta_idx][loge_idx] -> SetLineStyle(9);  //dashed
            phase1_curves_lower[seed_eta_idx][loge_idx] -> SetLineStyle(9);  //dashed
            phase1_dPhi_left[seed_eta_idx][loge_idx] -> SetLineStyle(9);
            phase1_dPhi_right[seed_eta_idx][loge_idx] -> SetLineStyle(9);

            legends[seed_eta_idx][loge_idx] -> AddEntry(fit_curves_upper[seed_eta_idx][loge_idx], "Upper - Optimized Averaged","l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(fit_curves_lower[seed_eta_idx][loge_idx], "Lower - Optimized Averaged","l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(phase1_curves_upper[seed_eta_idx][loge_idx], "Upper - Phase 1","l");
            legends[seed_eta_idx][loge_idx] -> AddEntry(phase1_curves_lower[seed_eta_idx][loge_idx], "Lower - Phase 1","l");
            legends_dPhi[seed_eta_idx][loge_idx] -> AddEntry(fit_dPhi_left[seed_eta_idx][loge_idx], "DPhi Window - Optimized","l");
            legends_dPhi[seed_eta_idx][loge_idx] -> AddEntry(phase1_dPhi_left[seed_eta_idx][loge_idx], "DPhi Window - Phase 1","l");

            caloClusters_shape_eBins_etWeight[seed_eta_idx][loge_idx] -> Draw("COLZ");
            legends[seed_eta_idx][loge_idx] -> Draw("SAME");
            legends_dPhi[seed_eta_idx][loge_idx] -> Draw("SAME");

            fit_curves_upper[seed_eta_idx][loge_idx] -> Draw("SAME");
            fit_curves_lower[seed_eta_idx][loge_idx] -> Draw("SAME");
            phase1_curves_upper[seed_eta_idx][loge_idx] -> Draw("SAME");
            phase1_curves_lower[seed_eta_idx][loge_idx] -> Draw("SAME");

            fit_dPhi_left[seed_eta_idx][loge_idx] -> Draw("SAME");
            fit_dPhi_right[seed_eta_idx][loge_idx] -> Draw("SAME");
            phase1_dPhi_left[seed_eta_idx][loge_idx] -> Draw("SAME");
            phase1_dPhi_right[seed_eta_idx][loge_idx] -> Draw("SAME");
            
            string outdir = "Output/Plots/";
            setOutput(outdir);

            c1->SaveAs((outdir+"caloShape_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".png").c_str());
            c1->SaveAs((outdir+"caloShape_Eta_"+filename_etas[seed_eta_idx]+"_"+filename_loges[loge_idx]+".pdf").c_str());
        }
    }
}

void ReadInfile(string inputFile){
    cout<<"Reading input file"<<endl;
    TFile *hist_infile = TFile::Open(inputFile.c_str());

    double seedEtaStep = (maxSeedEta - minSeedEta) / seedEtaBins;
    double logEStep = (maxLogE - minLogE) / logEBins;

    filename_etas.resize(seedEtaBins);
    filename_loges.resize(logEBins);

    double seedEtaVal = minSeedEta;
    for(int seedEtaIdx = 0; seedEtaIdx < seedEtaBins; seedEtaIdx++){
        titles_etas.push_back((to_string(seedEtaVal) + " < |#eta| #leq " + to_string(seedEtaVal + seedEtaStep)).c_str());
        midEtas.push_back((2*seedEtaVal) / 2.0);
        seedEtaVal+=seedEtaStep;

        std::ostringstream ss;
        ss << std::setw(3) << std::setfill('0') << to_string(seedEtaIdx);
        filename_etas[seedEtaIdx] = ss.str();
    }

    double logEVal = minLogE;
    for(int logEIdx = 0; logEIdx < logEBins; logEIdx++){
        titles_loge.push_back((to_string(logEVal) + " #leq log_{10}(E) < " + to_string(logEVal + logEStep)).c_str());
        midLogEs.push_back((2*logEVal) / 2.0);
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

    new_upperParams_local.resize(seedEtaBins, vector<vector<float>>(logEBins, vector<float>(4)));
    new_lowerParams_local.resize(seedEtaBins, vector<vector<float>>(logEBins, vector<float>(4)));

    cluster_dPhi_vs_loget.resize(dPhiWindowEtaBins);
    caloClusters_dPhiDist.resize(dPhiWindowEtaBins, vector<TH1F*>(dPhiWindowETDistBins));

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

    //hist_infile->Close();
    cout<<"Closed input file"<<endl;
}

void ReadInParameters(std::string localParamFileName, std::string aveParamFileName, std::string dPhiParamFileName){

    ifstream aveParam_infile;
    ifstream localParam_infile;
    ifstream dPhiParam_infile;

    aveParam_infile.open(aveParamFileName, std::ofstream::in | std::ofstream::trunc);
    localParam_infile.open(localParamFileName, std::ofstream::in | std::ofstream::trunc);
    dPhiParam_infile.open(dPhiParamFileName, std::ofstream::in | std::ofstream::trunc);
    

    std::string curLine;

    //Averaged Parameter File Format:
    //AVE_WEIGHT\tW00\tW01\tW10\tW11

    std::getline(aveParam_infile, curLine); //top line
    while(!aveParam_infile.eof()){
        getline(aveParam_infile, curLine);
        istringstream iss(curLine);
        vector<string> tokens;
        string token;
        while(std::getline(iss, token, '\t'))
            tokens.push_back(token);
        averageTitles.push_back(tokens.at(0));
        new_w00_averaged.push_back(std::stof(tokens.at(1)));
        new_w01_averaged.push_back(std::stof(tokens.at(2)));
        new_w10_averaged.push_back(std::stof(tokens.at(3)));
        new_w11_averaged.push_back(std::stof(tokens.at(4)));
    }

    
    //Local Region Parameter File Format:
    //SEED_ETA\tLOG(E)\tPARABOLA\tFROM_FIT\tW00\tW01\tW10\tW11
    
    std::getline(localParam_infile, curLine); //top line
    while(!localParam_infile.eof()){
        getline(localParam_infile, curLine);
        istringstream iss(curLine);
        vector<string> tokens;
        string token;
        while(std::getline(iss, token, '\t'))
            tokens.push_back(token);

        std::vector<double>::iterator seedEta = std::find(midEtas.begin(), midEtas.end(), std::stod(tokens.at(0)));
        std::vector<double>::iterator logE = std::find(midLogEs.begin(), midLogEs.end(), std::stod(tokens.at(1)));
        int seedEtaIdx = std::distance(midEtas.begin(), seedEta);
        int logEIdx = std::distance(midLogEs.begin(), logE);
        
        if(tokens.at(3) == "FIT") fitPlot[seedEtaIdx][logEIdx] = true;
        else fitPlot[seedEtaIdx][logEIdx] = false;

        if(tokens.at(2) == "UP"){
            new_upperParams_local[seedEtaIdx][logEIdx][0] = std::stof(tokens.at(4));
            new_upperParams_local[seedEtaIdx][logEIdx][1] = std::stof(tokens.at(5));
            new_upperParams_local[seedEtaIdx][logEIdx][2] = std::stof(tokens.at(6));
            new_upperParams_local[seedEtaIdx][logEIdx][3] = std::stof(tokens.at(7));
        }
        if(tokens.at(2) == "LOW"){
            new_lowerParams_local[seedEtaIdx][logEIdx][0] = std::stof(tokens.at(4));
            new_lowerParams_local[seedEtaIdx][logEIdx][1] = std::stof(tokens.at(5));
            new_lowerParams_local[seedEtaIdx][logEIdx][2] = std::stof(tokens.at(6));
            new_lowerParams_local[seedEtaIdx][logEIdx][3] = std::stof(tokens.at(7));
        }
    }


    //DPhi Parameter File Format:
    //ETA_BIN\tYOFFSET\tSCALE\tXOFFSET\tWIDTH

    std::getline(dPhiParam_infile, curLine);
    while(!dPhiParam_infile.eof()){
        getline(dPhiParam_infile, curLine);
        istringstream iss(curLine);
        vector<string> tokens;
        string token;
        while(std::getline(iss, token, '\t'))
            tokens.push_back(token);
        dPhiSections.push_back(tokens.at(0));
        new_yoffset.push_back(std::stof(tokens.at(1)));
        new_scale.push_back(std::stof(tokens.at(2)));
        new_xoffset.push_back(std::stof(tokens.at(3)));
        new_width.push_back(std::stof(tokens.at(4)));
    }

    aveParam_infile.close();
    localParam_infile.close();
    dPhiParam_infile.close();
}

void CaloShapePlotter(std::string inputFile, std::string localParamFileName, std::string aveParamFileName, std::string dPhiParamFileName){

    ReadInfile(inputFile);
    ReadInParameters(localParamFileName, aveParamFileName, dPhiParamFileName);

    parabolaCompareHeatMap();
    Plot();
    plotRegionalCurves();

}