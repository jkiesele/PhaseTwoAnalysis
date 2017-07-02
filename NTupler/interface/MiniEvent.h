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
    ngl=0; ngj=0;  
    nle=0; nte = 0; nlm=0; ntm=0; nj=0; nmet=0;
  }

  Int_t run,event,lumi;


  //gen level event
  Int_t ng,ngj,ngl;
  Float_t gl_p[50], gl_px[50], gl_py[50], gl_pz[50], gl_nrj[50], gl_pt[50], gl_eta[50], gl_phi[50], gl_mass[50], gl_relIso[50];
  Int_t gl_pid[50], gl_ch[50], gl_st[50];
  Float_t gj_pt[200], gj_eta[200], gj_phi[200], gj_mass[200];

  //reco level event
  Int_t nvtx;
  Float_t v_pt2[200];
  Int_t nle, nte, nlm, ntm, nj, nmet;
  Int_t le_ch[50], le_g[50];
  Float_t le_pt[50], le_eta[50], le_phi[50], le_mass[50], le_relIso[50];
  Int_t te_ch[50], te_g[50];
  Float_t te_pt[50], te_eta[50], te_phi[50], te_mass[50], te_relIso[50];
  Int_t lm_ch[50], lm_g[50];
  Float_t lm_pt[50], lm_eta[50], lm_phi[50], lm_mass[50], lm_relIso[50];
  Int_t tm_ch[50], tm_g[50];
  Float_t tm_pt[50], tm_eta[50], tm_phi[50], tm_mass[50], tm_relIso[50];
  Int_t j_id[200], j_g[200], j_csvv2[200], j_deepcsv[200], j_flav[200], j_hadflav[200], j_pid[200];
  Float_t j_pt[200], j_eta[200], j_phi[200], j_mass[200];
  Float_t met_pt[10], met_eta[10], met_phi[10];

};

void createMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_looseElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_,MiniEvent_t &ev);

#endif
