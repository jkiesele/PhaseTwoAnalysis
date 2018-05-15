#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"
#include "THStack.h"


void getRegion(int I)
{
    
    
  TString histName[] = {"detajj" ,"mht30","Mee","Mmumu","ele_mt","mu_mt","Jet1_pt","Jet2_pt","dphijj","metnolep"};
  //TString histName[] = {"detajj" ,"mht30","Jet1_pt","Jet2_pt","Jet1_eta","Jet2_eta","dphijj","metnolep"};
  
    /*  TString binRange[] = {"(100,0,400)","(100,-5,5)","(100,0,400)","(10,-5,5)","(100,0,5000)",
     "(100,0,10)","(100,0,3.1416)","(50,0,50)","(10,0,10)","(10,0,10)",
     "(10,0,10)","(10,0,10)","(10,0,10)","(100,0,400)","(100,-5,5)",
     "(4,0,4)","(21,-1,20)","(100,0,5)","(28,-6,22)","(10,0,5)",
     "(10,0,10)","(100,0,400)","(100,-5,5)","(4,0,4)","(21,-1,20)",
     "(100,0,5)","(10,0,5)","(28,-6,22)","(10,0,10)","(100,0,5000)",
     "(100,0,10)","(100,0,3.1416)","(100,0,500)","(100,0,3.1416)"};*/
    //SR BIN RANGE I = 1 TString binRange[] = {"(32,1.8,8.2)","(37,230,970)","(14,55,125)","(14,55,125)","(18,0,180)","(18,0,180)"};
   //  Z CR BIN RANGE I = 3,5 TString binRange[] = {"(8,4,8)","(20,200,1000)","(14,55,125)","(14,55,125)","(18,0,180)","(18,0,180)"};
    // W CR BIN RANGE I = 2,4
    TString binRange[] = {"(20,4.1,8.1)","(50,0,1000)","(14,55,125)","(14,55,125)","(18,0,180)","(18,0,180)","(100,0,400)","(100,0,400)","(20,0,3.1416)","(50,0,1000)"};
    //TString histLabelX[] = {"GenJet1_pt", "GenJet1_eta","GenJet2_pt","GenJet2_eta" ,"Gen_Mjj" , "Gen_detajj","Gen_dphijj" , "njets30","nlooseEle","ntightEle","nlooseMu","ntightMu","nlooseGamma","Jet1_pt" , "Jet1_eta" ,"Jet1_ID" ,"Jet1_genidx","Jet1_genjetDR","Jet1_parton","Jet1_JES" ,"hJet1_deepcsv","Jet2_pt","Jet2_eta" ,"Jet2_ID","Jet2_genidx" ,"Jet2_genjetDR","Jet2_JES", "Jet2_parton" , "Jet2_deepcsv" ,"Mjj" ,"detajj" ,"dphijj" , "metnolep","jetmetnolepmindphi" };
   // TString histLabelX[] = {"#Delta#eta_{jj}" ,"E_{T,miss (no lep)}","M_{ee}","M_{#mu#mu}","m_{T_{ele}}","m_{T_{#mu}}"};
    TString histLabelX[] = {"#Delta#eta_{jj}" ,"MHT","M_{ee}","M_{#mu#mu}","m_{T_{ele}}","m_{T_{#mu}}","P_{T_{j1}}","P_{T_{j2}}","#Delta#phi_{jj}"};
    TString Fillcolors[] = {"#cc0099","#99ccff"};
    TString Linecolors[] = {"#990066","#000099"};
    TString StackTitle = "Proba";
    // TString ProcessNames[] = {"VBFH","DYToLL-M-50_0J","DYToLL-M-50_1J","DYToLL-M-50_2J","DYToLL-M-50_3J","EWKWMinus2Jets_WToLNu_M-50","EWKWPlus2Jets_WToLNu_M-50","EWKZ2Jets_ZToLL_M-50","EWKZ2Jets_ZToNuNu","QCD_Mdijet-1000toInf","ST_tW_DR_14TeV_top","ST_tW_DR_14TeV_antitop","ST_tch_14TeV_antitop","ST_tch_14TeV_top","TT_TuneCUETP8M2T4","WToLNu_0J","WToLNu_1J","WToLNu_2J","WToLNu_3J","ZJetsToNuNu_HT-100To200","ZJetsToNuNu_HT-1200To2500","ZJetsToNuNu_HT-200To400","ZJetsToNuNu_HT-400To600","ZJetsToNuNu_HT-600To800"};
    
    
    
    
    TString ProcessNames[] = {"VBFH",
        "QCD_Mdijet-1000toInf",
        "EWKZ2Jets_ZToLL_M-50",
        "DYToLL-M-50_0J","DYToLL-M-50_1J","DYToLL-M-50_2J","DYToLL-M-50_3J",
        "ST_tW_DR_14TeV_top","ST_tW_DR_14TeV_antitop","ST_tch_14TeV_antitop","ST_tch_14TeV_top","TT_TuneCUETP8M2T4",
        "EWKWMinus2Jets_WToLNu_M-50","EWKWPlus2Jets_WToLNu_M-50",
        "EWKZ2Jets_ZToNuNu",
        "WToLNu_0J","WToLNu_1J","WToLNu_2J","WToLNu_3J",
        "ZJetsToNuNu_HT-100To200","ZJetsToNuNu_HT-1200To2500","ZJetsToNuNu_HT-200To400","ZJetsToNuNu_HT-400To600","ZJetsToNuNu_HT-600To800"};
    
    double N_gen_noPU[] = {480596,
        3958281,
        300000,
        2692863,2279465,2434045,726997,
        466371, 462190, 493128, 539267, 4743000,
        190404, 210952,
        244200,
        117171, 221561, 363169, 584737,
        138083, 97546, 97884, 128084, 104821
    };
    
    double N_gen_PU200[] = { 115660,
        3654010,
        251150,
        0, 2142302, 136936,
        453087, 500910, 123318, 416119, 2776312,
        151816, 144716,
        150600,
        223354, 69321, 236271, 277278,
        108564, 78895, 71915, 73684, 4166};
        
        
        
        
        
        
        
        
    
    
    
    
    
    TString LabelNames[] = {"VBFH","#DY#"};
    TString PU[] = {"noPU","PU200"};
    TString selection[] = {"Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5 ",
        "Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5   && nlooseEle==0 && nlooseMu==0",
        "Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5 && ntightEle==1 && nlooseMu==0",
        "Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5 && ntightEle==2 && nlooseMu==0",
        "Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5 && nlooseEle==0 && ntightMu==1",
        "Jet1_pt>80 && Jet2_pt>40 && Jet1_eta*Jet2_eta<0 && Mjj>1300 && detajj>4 && dphijj<1.5 && nlooseEle==0 && ntightMu==2 ",
        "Mjj>1300",
        "Mjj>1300&&jetmetnolepmindphi>0.5"};
    TString Regions[] = {"","0e0mu","1e0mu","2e0mu","0e1mu","0e2mu","mjj","mjj_jetmetnolepmindphi"};
    enum  {VBFH, QCD, EWK_Zll, QCD_Zll, Top, EWK_W_lnu, EWK_Z_nunu, QCD_W_lnu, QCD_Z_nunu};
    TString Yields_names [] = {"VBFH", "QCD", "Zll (EWK)", "Zll (QCD)", "Top", "W->lnu (EWK)", "Z->nunu (EWK)","W->lnu (QCD)", "Z->nunu (QCD)"};
    
    TH1F *histPlot,*dataPlot;
    
    int fileSize = 2;
    int plotSize = 2;
    
    TTree* LightTree = NULL;
    
    for (int regions = I; regions<I+1; regions++)
    {
        ofstream Yield("Yield_R"+Regions[regions]+".txt");
        Yield<<"Printing event Yields"<<endl;
        Yield<<"Format: N_sel, N_gen, w, w*N_sel, w*N_sel/N_gen"<<endl;
        TFile *outfile = TFile::Open("plots"+Regions[regions]+".root","RECREATE"); 
        for (int pu = 0; pu<2; pu++)
            //  for (int pu = 0; pu<1; pu++)
        {
            double Process_Yields[9];
            double Process_genYields[9];
            Yield<<PU[pu]<<endl;
            
            
            //for (int variable = 0; variable <34; variable++)
            for (int variable = 0; variable <10; variable++)
            {
                for (int yields = 0; yields<9;yields++)
                {
                    Process_Yields[yields] = 0;
                    Process_genYields[yields] = 0;
                }
                
                THStack *A = new THStack("test","test");
                TH1F *background;
                TLegend* leg = 0;
                double leg_xl = 0.65, leg_xr = 0.85, leg_yb = 0.55, leg_yt = 0.85 ;
                leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt);
                
                for(int i = 0;i < 24; i++)
                {
                    
                    //setLegendProperties(leg,"Control Region");
                    
                    
                    //TFile *_file = new TFile("HistosFile_VBFHtest_noPU.root");
                    TFile *_file = new TFile("/afs/cern.ch/work/a/amagnan/UPSGAna/180226/LooseVBFsel/"+PU[pu]+"/"+ProcessNames[i]+"/HistosFile_"+ProcessNames[i]+"_"+PU[pu]+".root");
                    _file->GetObject("LightTree",LightTree);
                    if (LightTree!=NULL){
                        
                        
                        if (i==0)
                        {
                            
                            
                            
                            double weight = 0;
                            // int a0  = LightTree->Draw("GenJet1_pt>>dataPlot(100,0,1000)","metnolep>0");
                            if((ProcessNames[i]+"_"+PU[pu]) == "VBFH_noPU") weight = 4.278*35900.0/480596;
                            if((ProcessNames[i]+"_"+PU[pu]) == "VBFH_PU200") weight = 4.278*35900.0/115660;
                            // cout<<"vbf weight ="<<weight<<endl;
                            int counters1  = LightTree->Draw(histName[variable]+">>SignalPlot"+PU[pu]+binRange[variable],selection[regions]);
                            dataPlot= (TH1F*)gDirectory->Get("SignalPlot"+PU[pu]);
                            
                            for (int j = 0; j<dataPlot->GetNbinsX();j++)
                            {
                                dataPlot->SetBinContent(j,dataPlot->GetBinContent(j)*weight);
                                
                            }
                            dataPlot -> SetLineColor(kBlack);
                            leg->AddEntry(dataPlot,"Signal","LP");
                            if (variable==0) Process_Yields[VBFH]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            
                            
                        }
                        else
                        {
                            
                            int counters2  = LightTree->Draw(histName[variable]+">>"+ProcessNames[i]+"_"+PU[pu]+binRange[variable],selection[regions]);
                            histPlot= (TH1F*)gDirectory->Get(ProcessNames[i]+"_"+PU[pu]);
                            /*   if ((i>=1)&&(i<=4)){
                             histPlot -> SetFillColor(kGreen);
                             histPlot -> SetLineColor(kGreen-2);
                             if (i==1) leg->AddEntry(histPlot,"DY","F");
                             }
                             if ((i>=5)&&(i<=8)){
                             histPlot -> SetFillColor(kBlue);
                             histPlot -> SetLineColor(kBlue-2);
                             leg->AddEntry(histPlot,ProcessNames[i],"F");
                             }
                             
                             if (i==9){
                             histPlot -> SetFillColor(kBlack+30);
                             histPlot -> SetLineColor(kBlack+28);
                             leg->AddEntry(histPlot,ProcessNames[i],"F");
                             }
                             
                             if ((i>=10)&&(i<=13)){
                             histPlot -> SetFillColor(kOrange);
                             histPlot -> SetLineColor(kOrange-2);
                             if (i==10) leg->AddEntry(histPlot,"top/antitop","F");
                             }
                             if (i==14){
                             histPlot -> SetFillColor(kBlue+10);
                             histPlot -> SetLineColor(kBlue+12);
                             leg->AddEntry(histPlot,ProcessNames[i],"F");
                             }
                             if ((i>=15)&&(i<=18)){
                             histPlot -> SetFillColor(kMagenta);
                             histPlot -> SetLineColor(kMagenta-2);
                             if (i==15) leg->AddEntry(histPlot,"WToLNu","F");
                             }
                             
                             if ((i>=19)&&(i<=23)){
                             histPlot -> SetFillColor(kGreen+15);
                             histPlot -> SetLineColor(kGreen+12);
                             if (i==19) leg->AddEntry(histPlot,"ZJetsToNuNu","F");
                             }*/
                            double weight =0 ;
                            
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_0J_noPU") weight = 3728*35900.0/2692863;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_1J_noPU") weight = 1096*35900.0/2279465;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_2J_noPU") weight = 384*35900.0/2434045;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_3J_noPU") weight = 165.6*35900.0/726997;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKWMinus2Jets_WToLNu_M-50_noPU") weight = 23.6*35900.0/190404;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKWPlus2Jets_WToLNu_M-50_noPU") weight = 30.37*35900.0/210952;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKZ2Jets_ZToLL_M-50_noPU") weight = 4.38*35900.0/300000;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKZ2Jets_ZToNuNu_noPU") weight = 11.05*35900.0/244200;
                            if((ProcessNames[i]+"_"+PU[pu]) == "QCD_Mdijet-1000toInf_noPU") weight = 1135*35900.0/3958281;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tW_DR_14TeV_antitop_noPU") weight = 45.06*35900.0/462190;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tW_DR_14TeV_top_noPU") weight = 45.06*35900.0/466371;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tch_14TeV_antitop_noPU") weight = 29.2*35900.0/493128;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tch_14TeV_top_noPU") weight = 48.03*35900.0/539267;
                            if((ProcessNames[i]+"_"+PU[pu]) == "TT_TuneCUETP8M2T4_noPU") weight = 984.5*35900.0/4933256;
                            if((ProcessNames[i]+"_"+PU[pu]) == "VBFH_noPU") weight = 4.278*35900.0/480596;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_0J_noPU") weight = 38870*35900.0/117171;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_1J_noPU") weight = 10330*35900.0/221561;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_2J_noPU") weight = 3314*35900.0/363169;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_3J_noPU") weight = 1891*35900.0/584737;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-100To200_noPU") weight = 304.5*35900.0/138083;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-1200To2500_noPU") weight = 0.345*35900.0/97546;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-200To400_noPU") weight = 85.88*35900.0/97884;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-2500ToInf_noPU") weight = 0.0069*35900.0/0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-2500ToInf_noPU") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-400To600_noPU") weight = 12.35*35900.0/128084;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-600To800_noPU") weight = 3.029*35900.0/104821;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-800To1200_noPU") weight = 1.400*35900.0/68672;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-800To1200_noPU") weight = 0;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_0J_PU200") weight = 3728*35900.0/0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_0J_PU200") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_1J_PU200") weight = 1096*35900.0/2142302;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_2J_PU200") weight = 384*35900.0/2659054;
                            if((ProcessNames[i]+"_"+PU[pu]) == "DYToLL-M-50_3J_PU200") weight = 165.6*35900.0/136936;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKWMinus2Jets_WToLNu_M-50_PU200") weight = 23.6*35900.0/151816;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKWPlus2Jets_WToLNu_M-50_PU200") weight = 30.37*35900.0/144716;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKZ2Jets_ZToLL_M-50_PU200") weight = 4.38*35900.0/251150;
                            if((ProcessNames[i]+"_"+PU[pu]) == "EWKZ2Jets_ZToNuNu_PU200") weight = 11.05*35900.0/150600;
                            if((ProcessNames[i]+"_"+PU[pu]) == "QCD_Mdijet-1000toInf_PU200") weight = 1135*35900.0/3654010;
                            // if((ProcessNames[i]+"_"+PU[pu]) == "QCD_Mdijet-1000toInf_PU200") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tW_DR_14TeV_antitop_PU200") weight = 45.06*35900.0/500910;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ST_tW_DR_14TeV_antitop_PU200") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tW_DR_14TeV_top_PU200") weight = 45.06*35900.0/453087;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tch_14TeV_antitop_PU200") weight = 29.2*35900.0/123318;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ST_tch_14TeV_top_PU200") weight = 48.03*35900.0/416119;
                            if((ProcessNames[i]+"_"+PU[pu]) == "TT_TuneCUETP8M2T4_PU200") weight = 984.5*35900.0/2776312;
                            if((ProcessNames[i]+"_"+PU[pu]) == "VBFH_PU200") weight = 4.278*35900.0/115660;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_0J_PU200") weight = 38870*35900.0/223354;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_1J_PU200") weight = 10330*35900.0/69321;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_2J_PU200") weight = 3314*35900.0/236271;
                            if((ProcessNames[i]+"_"+PU[pu]) == "WToLNu_3J_PU200") weight = 1891*35900.0/277278;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-100To200_PU200") weight = 304.5*35900.0/108564;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-1200To2500_PU200") weight = 0.345*35900.0/78895;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-200To400_PU200") weight = 85.88*35900.0/71915;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-2500ToInf_PU200") weight = 0.0069*35900.0/0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-2500ToInf_PU200") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-400To600_PU200") weight = 12.35*35900.0/73684;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-600To800_PU200") weight = 3.029*35900.0/4166;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-600To800_PU200") weight = 0;
                            if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-800To1200_PU200") weight = 1.400*35900.0/49011;
                            //if((ProcessNames[i]+"_"+PU[pu]) == "ZJetsToNuNu_HT-800To1200_PU200") weight = 0;
                            
                            if (variable==0) {
			      cout<<ProcessNames[i]+"_"+PU[pu]<<" = "<< weight<<endl;
			    }
			    if (weight > 1000) {
			      weight = 0;
			      if (variable==0)std::cout << " -- Weight too large, ignoring this sample... setting weight to 0" << std::endl;
			    }
                            
                            
                            
                            if ((i>=1)&&(i<2)){
                                histPlot -> SetFillColor(kWhite);
                                histPlot -> SetLineColor(kWhite-1);
                                if (i==1) leg->AddEntry(histPlot,"QCD","F");
                                if (variable==0) Process_Yields[QCD]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=2)&&(i<3)){
                                histPlot -> SetFillColor(TColor::GetColor("#91ABC4"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==2) leg->AddEntry(histPlot,"EWK Zll","F");
                                if (variable==0) Process_Yields[EWK_Zll]+= 1.0*weight* LightTree->GetEntries(selection[regions]);
                            }
                            
                            if ((i>=3)&&(i<7)){
                                histPlot -> SetFillColor(TColor::GetColor("#9A9EAB"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==3) leg->AddEntry(histPlot,"QCD Zll","F");
                                if (variable==0) Process_Yields[QCD_Zll]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=7)&&(i<12)){
                                histPlot -> SetFillColor(TColor::GetColor("#CF3721"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==7) leg->AddEntry(histPlot,"Top","F");
                                if (variable==0) Process_Yields[Top]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=12)&&(i<14)){
                                histPlot -> SetFillColor(TColor::GetColor("#0066CC"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==12) leg->AddEntry(histPlot,"EWK W+jets->lnu","F");
                                if (variable==0) Process_Yields[EWK_W_lnu]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=14)&&(i<15)){
                                histPlot -> SetFillColor(TColor::GetColor("#00CCCC"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==14) leg->AddEntry(histPlot,"EWK Z->nunu","F");
                                if (variable==0) Process_Yields[EWK_Z_nunu]+= 1.0*weight* LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=15)&&(i<19)){
                                histPlot -> SetFillColor(TColor::GetColor("#E19D07"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==15) leg->AddEntry(histPlot,"QCD W+jets->lnu","F");
                                if (variable==0) Process_Yields[QCD_W_lnu]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            if ((i>=19)&&(i<24)){
                                histPlot -> SetFillColor(TColor::GetColor("#4D975D"));
                                histPlot -> SetLineColor(kBlack);
                                if (i==19) leg->AddEntry(histPlot,"QCD Z->nunu","F");
                                if (variable==0) Process_Yields[QCD_Z_nunu]+= 1.0*weight*LightTree->GetEntries(selection[regions]);
                            }
                            
                            
                            
                            
                            
                            
                            
                            
                            for (int j = 0; j<histPlot->GetNbinsX();j++)
                            {
                                histPlot->SetBinContent(j,histPlot->GetBinContent(j)*weight);
                                
                            }
                            histPlot -> SetFillStyle(1001);
                            histPlot->SetLineStyle(0);
                            histPlot->SetMarkerStyle(20);
                            histPlot->GetXaxis()->SetLabelFont(42);
                            histPlot->GetXaxis()->SetLabelOffset(0.007);
                            histPlot->GetXaxis()->SetLabelSize(0.05);
                            histPlot->GetXaxis()->SetTitleSize(0.06);
                            histPlot->GetXaxis()->SetTitleOffset(0.9);
                            histPlot->GetXaxis()->SetTitleFont(42);
                            histPlot->GetYaxis()->SetLabelFont(42);
                            histPlot->GetYaxis()->SetLabelOffset(0.007);
                            histPlot->GetYaxis()->SetLabelSize(0.05);
                            histPlot->GetYaxis()->SetTitleSize(0.06);
                            histPlot->GetYaxis()->SetTitleOffset(1.25);
                            histPlot->GetYaxis()->SetTitleFont(42);
                            histPlot->GetZaxis()->SetLabelFont(42);
                            histPlot->GetZaxis()->SetLabelOffset(0.007);
                            histPlot->GetZaxis()->SetLabelSize(0.05);
                            histPlot->GetZaxis()->SetTitleSize(0.06);
                            histPlot->GetZaxis()->SetTitleFont(42);
                            histPlot->GetXaxis()->SetTitle(histLabelX[variable]+"(GeV)");
                            
                            A->Add(histPlot);
                            if (i==1) background = (TH1F*) histPlot->Clone("totalBkg");
                            else background->Add(histPlot);
                            
                            
                        }
                    }
                    
                  //  _file->Close();
                }
                
                
                TCanvas *c1 = new TCanvas("c1", "myPlots",0,67,600,600);
                gStyle->SetOptFit(1);
                gStyle->SetOptStat(0);
                gStyle->SetOptTitle(0);
                c1->Range(-102.5,-10.38415,847.5,69.4939);
                c1->SetFillColor(0);
                c1->SetBorderMode(0);
                c1->SetBorderSize(2);
                c1->SetTickx(1);
                c1->SetTicky(1);
                c1->SetLeftMargin(0.15);
                c1->SetRightMargin(0.05);
                c1->SetTopMargin(0.05);
                c1->SetBottomMargin(0.13);
                c1->SetFrameFillStyle(0);
                c1->SetFrameBorderMode(0);
                c1->SetFrameFillStyle(0);
                c1->SetFrameBorderMode(0);
                
                
                dataPlot->SetFillColor(kBlack);
                dataPlot->SetLineColor(kBlack);
                dataPlot->SetFillStyle(0);
                dataPlot->SetMarkerStyle(20);
                if (histName[variable]!="detajj")
                    dataPlot->GetXaxis()->SetTitle(histLabelX[variable]+"(GeV)");
                else
                    dataPlot->GetXaxis()->SetTitle(histLabelX[variable]);
                    
                dataPlot->GetXaxis()->SetLabelFont(42);
                dataPlot->GetXaxis()->SetLabelSize(0.05);
                dataPlot->GetXaxis()->SetTitleSize(0.06);
                dataPlot->GetXaxis()->SetTitleOffset(0.9);
                dataPlot->GetXaxis()->SetTitleFont(42);
                dataPlot->GetYaxis()->SetTitle("Events");
                dataPlot->GetYaxis()->SetLabelFont(42);
                dataPlot->GetYaxis()->SetLabelSize(0.05);
                dataPlot->GetYaxis()->SetTitleSize(0.06);
                dataPlot->GetYaxis()->SetTitleOffset(1.25);
                dataPlot->GetYaxis()->SetTitleFont(42);
                dataPlot->GetZaxis()->SetLabelFont(42);
                dataPlot->GetZaxis()->SetLabelSize(0.035);
                dataPlot->GetZaxis()->SetTitleSize(0.035);
                dataPlot->GetZaxis()->SetTitleFont(42);
                // dataPlot->Draw("e1 goff");
                
                dataPlot->Draw("e1 goff");
                dataPlot->GetYaxis()->SetRangeUser(0,1.1*A->GetMaximum());
                A->Draw("same hist goff");
                dataPlot->Draw("same e1 goff");
                gPad->RedrawAxis();
                leg->Draw("goff");
                
                TLatex *   tex = new TLatex(0.95,0.96,"35.9 fb^{-1} (14 TeV)");
                tex->SetNDC();
                tex->SetTextAlign(31);
                tex->SetTextFont(42);
                tex->SetTextSize(0.03);
                tex->SetLineWidth(2);
                tex->Draw("goff");
                tex = new TLatex(0.15,0.96,"CMS");
                tex->SetNDC();
                tex->SetTextFont(61);
                tex->SetTextSize(0.0375);
                tex->SetLineWidth(2);
                tex->Draw("goff");
                tex = new TLatex(0.23,0.96,"Preliminary "+PU[pu]);
                tex->SetNDC();
                tex->SetTextFont(52);
                tex->SetTextSize(0.0285);
                tex->SetLineWidth(2);
                tex->Draw("goff");
                
                outfile->cd();
                dataPlot->Write("Signal_"+histName[variable]+"_"+PU[pu]);
                background->Write("Background_"+histName[variable]+"_"+PU[pu]);
                
                c1->SaveAs("plots"+Regions[regions]+"/"+histName[variable]+"_"+PU[pu]+"_auto.pdf");
                c1->SaveAs("plots"+Regions[regions]+"/"+histName[variable]+"_"+PU[pu]+"_auto.C");
                c1->SaveAs("plots"+Regions[regions]+"/"+histName[variable]+"_"+PU[pu]+"_auto.root");
                
                if (variable==0){
                    for (int yields = VBFH;yields<=QCD_Z_nunu; yields++)
                        Yield<<Yields_names[yields]<<" :="<<Process_Yields[yields]<<endl;
                }
                
                //cout<<"save"<<endl;
            }//loop on variables
            Yield<<"  "<<endl;
        }//loop on PU
        Yield.close();
	outfile->Close();
    }//loop on regions

}//getRegion method

int plotter(){

  getRegion(1);
  return 0;

}
