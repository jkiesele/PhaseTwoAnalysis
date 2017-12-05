#ifndef _minievent_h_
#define _minievent_h_
// -*- C++ -*-
//
// Package:     PhaseTwoAnalysis/NTupler
// Class:       MiniEvent
// Description: Define the structure of ntuples

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    ngl=0; ngj=0; ngp=0;
    g_nw=0;ng=0;
    nle=0; nme=0; nte = 0; nlm=0; ntm=0; nj=0; nmet=0; nlp=0; ntp=0;
    lumi=0;
    nvtx=0;
    run=0;
    event=0;
  }
  static constexpr int maxpart=50;
  static constexpr int maxjets=200;
  static constexpr int maxweights=500;

  Int_t run,event,lumi;

  //gen level event
  Int_t ng,ngj,ngl,ngp,g_nw;
  Float_t g_w[maxweights];
  Float_t gl_p[maxpart], gl_px[maxpart], gl_py[maxpart], gl_pz[maxpart], gl_nrj[maxpart], gl_pt[maxpart], gl_eta[maxpart], gl_phi[maxpart], gl_mass[maxpart], gl_relIso[maxpart];
  Int_t gl_pid[maxpart], gl_ch[maxpart], gl_st[maxpart];
  Float_t gj_pt[maxjets], gj_eta[maxjets], gj_phi[maxjets], gj_mass[maxjets];
  Float_t gp_p[maxpart], gp_px[maxpart], gp_py[maxpart], gp_pz[maxpart], gp_nrj[maxpart], gp_pt[maxpart], gp_eta[maxpart], gp_phi[maxpart];
  Int_t gp_st[maxpart];

  //reco level event
  Int_t nvtx;
  Float_t v_pt2[maxjets];
  Int_t nle, nme, nte, nlm, ntm, nj, nmet, nlp, ntp;
  Int_t le_ch[maxpart], le_g[maxpart];
  Float_t le_pt[maxpart], le_eta[maxpart], le_phi[maxpart], le_mass[maxpart], le_relIso[maxpart], le_bdt[maxpart];
  Int_t me_ch[maxpart], me_g[maxpart];
  Float_t me_pt[maxpart], me_eta[maxpart], me_phi[maxpart], me_mass[maxpart], me_relIso[maxpart], me_bdt[maxpart];
  Int_t te_ch[maxpart], te_g[maxpart];
  Float_t te_pt[maxpart], te_eta[maxpart], te_phi[maxpart], te_mass[maxpart], te_relIso[maxpart], te_bdt[maxpart];
  Int_t lm_ch[maxpart], lm_g[maxpart];
  Float_t lm_pt[maxpart], lm_eta[maxpart], lm_phi[maxpart], lm_mass[maxpart], lm_relIso[maxpart];
  Int_t tm_ch[maxpart], tm_g[maxpart];
  Float_t tm_pt[maxpart], tm_eta[maxpart], tm_phi[maxpart], tm_mass[maxpart], tm_relIso[maxpart];
  Int_t j_id[maxjets], j_g[maxjets], j_mvav2[maxjets], j_deepcsv[maxjets], j_flav[maxjets], j_hadflav[maxjets], j_pid[maxjets];
  Float_t j_pt[maxjets], j_eta[maxjets], j_phi[maxjets], j_mass[maxjets];
  Float_t met_pt[maxpart], met_eta[maxpart], met_phi[maxpart];
  Int_t lp_g[maxpart], tp_g[maxpart], lp_isEB[maxpart], tp_isEB[maxpart];
  Float_t lp_pt[maxpart], lp_eta[maxpart], lp_phi[maxpart], lp_nrj[maxpart], lp_pt_multi[maxpart], lp_eta_multi[maxpart], lp_phi_multi[maxpart], lp_nrj_multi[maxpart], lp_bdt[maxpart];
  Float_t tp_pt[maxpart], tp_eta[maxpart], tp_phi[maxpart], tp_nrj[maxpart], tp_pt_multi[maxpart], tp_eta_multi[maxpart], tp_phi_multi[maxpart], tp_nrj_multi[maxpart], tp_bdt[maxpart];

};

void createMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_genPhotons_, TTree *t_looseElecs_, TTree *t_mediumElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_, TTree *t_loosePhotons_, TTree *t_tightPhotons_, MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_genPhotons_, TTree *t_looseElecs_, TTree *t_mediumElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_, TTree *t_loosePhotons_, TTree *t_tightPhotons_, MiniEvent_t &ev);

#endif
