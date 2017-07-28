/*
 * TDRStyle.h
 *
 *  Created on: 28 Jul 2017
 *      Author: jkiesele
 */

#ifndef TDRSTYLE_H_
#define TDRSTYLE_H_

#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

static inline TCanvas * setCanvas(bool split=false){

    int yrange=600;
    if(split)
        yrange=800;
    TCanvas * can = new TCanvas("can", "can", 800,yrange);
    if (split){
        can->Divide(1, 2)                                                ;
        can->GetPad(1)->SetPad("Top", "", 0., 0.25, 1.0, 1.0, 0, -1, 0)  ;
        can->GetPad(1)->SetTopMargin(0.069)                              ;
        can->GetPad(1)->SetBottomMargin(0.0184)                          ;
        can->GetPad(1)->SetRightMargin(0.046)                            ;
        can->GetPad(1)->SetLeftMargin(0.138)                             ;
        can->GetPad(1)->SetTicks(1, 1)                                   ;
        can->GetPad(2)->SetPad("Bottom", "", 0., 0., 1.0, 0.25, 0, -1, 0);
        can->GetPad(2)->SetTopMargin(0.0092)                             ;
        can->GetPad(2)->SetBottomMargin(0.368)                           ;
        can->GetPad(2)->SetRightMargin(0.046)                            ;
        can->GetPad(2)->SetLeftMargin(0.138)                             ;
        can->GetPad(2)->SetTicks(1, 1)                                   ;
    }
    else{
        can->GetPad(0)->SetTopMargin(0.069)    ;
        can->GetPad(0)->SetRightMargin(0.046)  ;
        can->GetPad(0)->SetLeftMargin(0.138)   ;
        can->GetPad(0)->SetBottomMargin(0.15)  ;
        can->GetPad(0)->SetTicks(1, 1)         ;
    }
    can->cd(1);

    return can;

}


static inline void formatHisto(TH1 * hist){
    hist->GetXaxis()->SetTitleSize(0.0462874993682)   ;
    hist->GetXaxis()->SetTitleOffset(1.0)             ;
    hist->GetXaxis()->SetLabelSize(0.0462874993682)   ;
    hist->GetXaxis()->SetLabelOffset(0.0100567853078) ;

    hist->GetYaxis()->SetTitleSize(0.0462874993682)   ;
    hist->GetYaxis()->SetTitleOffset(1.32249999046)   ;
    hist->GetYaxis()->SetLabelSize(0.0462874993682)   ;
    hist->GetYaxis()->SetLabelOffset(0.005)           ;

    hist->GetZaxis()->SetTitleSize(0.0462874993682)   ;
    hist->GetZaxis()->SetTitleOffset(1.14999997616)   ;
    hist->GetZaxis()->SetLabelSize(0.0462874993682)   ;
}
static inline void formatRatio(TH1* h){
    h->GetXaxis()->SetTitleSize(0.138862490654)      ;
    h->GetXaxis()->SetTitleOffset(1.0)               ;
    h->GetXaxis()->SetLabelSize(0.138862490654)      ;
    h->GetXaxis()->SetLabelOffset(0.0150851774961)   ;

    h->GetYaxis()->SetTitleSize(0.138862490654)      ;
    h->GetYaxis()->SetLabelSize(0.138862490654)      ;
    h->GetYaxis()->SetTitleOffset(0.440833330154)    ;
    h->GetYaxis()->SetTitleOffset(0.440833330154)    ;

    h->GetYaxis()->SetNdivisions(505)                ;
    h->GetYaxis()->SetRangeUser(0.4, 1.6)            ;
}

static inline void drawCMS(bool onTop=false){
    TString text="Phase-2 Simulation";
    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(62);
    latex->SetTextSize(0.05175);
    if(onTop)
        latex->DrawLatex(0.18, 0.94, "CMS");
    else
        latex->DrawLatex(0.18, 0.85, "CMS");
    latex->SetTextSize(0.0414);
    latex->SetTextFont(52);
    if(onTop)
        latex->DrawLatex(0.26625,  0.94, text);
    else
        latex->DrawLatex(0.26625, 0.85, text);
}
static inline void drawEnPu(TString pileup="", TString lumi=""){
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.046);
    latex->SetTextColor(1)  ;
    latex->SetTextFont(42)  ;
    latex->SetTextAlign(31) ;
    TString tex = "14 TeV";
    if (pileup.Length()) tex += ", "+pileup+" PU";
    if (lumi.Length()) tex = lumi +", "+tex;
    latex->DrawLatex(0.95, 0.94, tex);

}

#endif /* TDRSTYLE_H_ */
