#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    ngl=0; ngj=0; ngpf=0; 
    nl=0; nj=0; nmet=0; npf=0;
  }

  Int_t run,event,lumi;


  //gen level event
  Int_t ng,ngj,ngl,ngpf;
  Float_t gl_pt[50], gl_eta[50], gl_phi[50], gl_mass[50], gl_relIso[50];
  Float_t gj_pt[200], gj_eta[200], gj_phi[200], gj_mass[200];
  Int_t gj_isb[200];
  Int_t gpf_id[5000];
  Float_t gpf_pt[5000], gpf_eta[5000], gpf_phi[5000], gpf_mass[5000];

  //reco level event
  Int_t nvtx;
  Int_t nl, nj, npf, nmet;
  Int_t l_id[50], l_pid[50], l_g[50];
  Float_t l_pt[50], l_eta[50], l_phi[50], l_mass[50], l_relIso[50];
  Int_t j_g[200], j_flav[200], j_hadflav[200], j_pid[200];
  Float_t j_pt[200], j_eta[200], j_phi[200], j_mass[200], j_csvv2[200], j_deepcsv[200];
  Int_t pf_j[5000], pf_jnpf[5000], pf_jnhppf[5000], pf_id[5000], pf_muon[5000];
  Float_t pf_pt[5000], pf_eta[5000], pf_phi[5000], pf_mass[5000];
  Bool_t pf_hp[5000];
  Float_t met_pt[10],met_phi[10];

};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev,bool full=false);

#endif
