#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->Branch("run",        &ev.run,       "run/I");
  t->Branch("event",      &ev.event,     "event/I");
  t->Branch("lumi",       &ev.lumi,      "lumi/I");

  //gen level event
  t->Branch("ngl",        &ev.ngl,       "ngl/I");
  t->Branch("gl_pid",     ev.gl_pid,     "gl_pid[ngl]/I");
  t->Branch("gl_pt",      ev.gl_pt,      "gl_pt[ngl]/F");
  t->Branch("gl_eta",     ev.gl_eta,     "gl_eta[ngl]/F");
  t->Branch("gl_phi",     ev.gl_phi,     "gl_phi[ngl]/F");
  t->Branch("gl_mass",    ev.gl_mass,    "gl_mass[ngl]/F");
  t->Branch("gl_relIso",  ev.gl_relIso,  "gl_relIso[ngl]/F");

  t->Branch("ngj",        &ev.ngj,       "ngj/I");
  t->Branch("gj_pt",      ev.gj_pt,      "gj_pt[ngj]/F");
  t->Branch("gj_eta",     ev.gj_eta,     "gj_eta[ngj]/F");
  t->Branch("gj_phi",     ev.gj_phi,     "gj_phi[ngj]/F");
  t->Branch("gj_mass",    ev.gj_mass,    "gj_mass[ngj]/F");

  //reco level event
  t->Branch("nvtx",       &ev.nvtx,      "nvtx/I");

  t->Branch("nl",         &ev.nl,        "nl/I");
  t->Branch("l_id",       ev.l_id,       "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,      "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,        "l_g[nl]/I");
  t->Branch("l_pt",       ev.l_pt,       "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,      "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,      "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,     "l_mass[nl]/F");
  t->Branch("l_relIso",   ev.l_relIso,   "l_relIso[nl]/F");

  t->Branch("nj",         &ev.nj,        "nj/I");
  t->Branch("j_id",       ev.j_id,       "j_id[nj]/I");
  t->Branch("j_g",        ev.j_g,        "j_g[nj]/I");
  t->Branch("j_pt",       ev.j_pt,       "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,      "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,      "j_phi[nj]/F");
  t->Branch("j_mass",     ev.j_mass,     "j_mass[nj]/F");
  t->Branch("j_csvv2",    ev.j_csvv2,    "j_csvv2[nj]/F");
  t->Branch("j_deepcsv",  ev.j_deepcsv,  "j_deepcsvl[nj]/F");
  t->Branch("j_flav",     ev.j_flav,     "j_flav[nj]/I");
  t->Branch("j_hadflav",  ev.j_hadflav,  "j_hadflav[nj]/I");
  t->Branch("j_pid",      ev.j_pid,      "j_pid[nj]/I");

  t->Branch("npf",        &ev.npf,       "npf/I");
  t->Branch("pf_pid",     ev.pf_pid,     "pf_pid[npf]/I");
  t->Branch("pf_pt",      ev.pf_pt,      "pf_pt[npf]/F");
  t->Branch("pf_eta",     ev.pf_eta,     "pf_eta[npf]/F");
  t->Branch("pf_phi",     ev.pf_phi,     "pf_phi[npf]/F");
  t->Branch("pf_mass",    ev.pf_mass,    "pf_mass[npf]/F");
  t->Branch("pf_relIso",  ev.pf_relIso,  "pf_relIso[npf]/F");
  t->Branch("pf_hp",      ev.pf_hp,      "pf_hp[npf]/O");

  t->Branch("nmet",       &ev.nmet,      "nmet/I");
  t->Branch("met_pt",     ev.met_pt,     "met_pt[nmet]/F");
  t->Branch("met_phi",    ev.met_phi,    "met_phi[nmet]/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->SetBranchAddress("run",         &ev.run);
  t->SetBranchAddress("event",       &ev.event);
  t->SetBranchAddress("lumi",        &ev.lumi);

  //gen level event 
  t->SetBranchAddress("ngl",         &ev.ngl);
  t->SetBranchAddress("gl_pid",      ev.gl_pid);
  t->SetBranchAddress("gl_pt",       ev.gl_pt);
  t->SetBranchAddress("gl_eta",      ev.gl_eta);
  t->SetBranchAddress("gl_phi",      ev.gl_phi);
  t->SetBranchAddress("gl_mass",     ev.gl_mass);
  t->SetBranchAddress("gl_relIso",   ev.gl_relIso);

  t->SetBranchAddress("ngj",         &ev.ngj);
  t->SetBranchAddress("gj_pt",       ev.gj_pt);
  t->SetBranchAddress("gj_eta",      ev.gj_eta);
  t->SetBranchAddress("gj_phi",      ev.gj_phi);
  t->SetBranchAddress("gj_mass",     ev.gj_mass);

  //reco level event
  t->SetBranchAddress("nvtx",        &ev.nvtx);

  t->SetBranchAddress("nl",          &ev.nl);
  t->SetBranchAddress("l_id",        ev.l_id);
  t->SetBranchAddress("l_pid",       ev.l_pid);
  t->SetBranchAddress("l_g",         ev.l_g);
  t->SetBranchAddress("l_pt",        ev.l_pt);
  t->SetBranchAddress("l_eta",       ev.l_eta);
  t->SetBranchAddress("l_phi",       ev.l_phi);
  t->SetBranchAddress("l_mass",      ev.l_mass);
  t->SetBranchAddress("l_relIso",           ev.l_relIso);

  t->SetBranchAddress("nj",          &ev.nj);
  t->SetBranchAddress("j_id",        ev.j_id);
  t->SetBranchAddress("j_g",         ev.j_g);
  t->SetBranchAddress("j_pt",        ev.j_pt);
  t->SetBranchAddress("j_eta",       ev.j_eta);
  t->SetBranchAddress("j_phi",       ev.j_phi);
  t->SetBranchAddress("j_mass",      ev.j_mass);
  t->SetBranchAddress("j_csvv2",     ev.j_csvv2);
  t->SetBranchAddress("j_deepcsv",   ev.j_deepcsv);
  t->SetBranchAddress("j_flav",      ev.j_flav);
  t->SetBranchAddress("j_hadflav",   ev.j_hadflav);
  t->SetBranchAddress("j_pid",       ev.j_pid);

  t->SetBranchAddress("npf",       &ev.npf);
  t->SetBranchAddress("pf_pid",    ev.pf_pid);
  t->SetBranchAddress("pf_pt",     ev.pf_pt);
  t->SetBranchAddress("pf_eta",    ev.pf_eta);
  t->SetBranchAddress("pf_phi",    ev.pf_phi);
  t->SetBranchAddress("pf_mass",   ev.pf_mass);
  t->SetBranchAddress("pf_relIso", ev.pf_relIso);
  t->SetBranchAddress("pf_hp",     ev.pf_hp);

  t->SetBranchAddress("nmet",        &ev.nmet);
  t->SetBranchAddress("met_pt",      ev.met_pt);
  t->SetBranchAddress("met_phi",     ev.met_phi);
}
