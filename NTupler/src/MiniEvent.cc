#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

void createMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_genPhotons_,
		TTree *t_looseElecs_, TTree *t_mediumElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_, TTree *t_loosePhotons_, TTree *t_tightPhotons_, MiniEvent_t &ev)
{
  //event header
  t_event_->Branch("Run",               &ev.run,        "Run/I");
  t_event_->Branch("Event",             &ev.event,      "Event/I");
  t_event_->Branch("Lumi",              &ev.lumi,       "Lumi/I");

  //gen level event
  t_event_->Branch("Weight_size",       &ev.g_nw,       "Weight_size/I");
  t_event_->Branch("Weight",            ev.g_w,         "Weight[Weight_size]/F");

  t_genParts_->Branch("Particle_size",  &ev.ngl,        "Particle_size/I");
  t_genParts_->Branch("PID",            ev.gl_pid,      "PID[Particle_size]/I");
  t_genParts_->Branch("Charge",         ev.gl_ch,       "Charge[Particle_size]/I");
  t_genParts_->Branch("Status",         ev.gl_st,       "Status[Particle_size]/I");
  t_genParts_->Branch("P",              ev.gl_p,        "P[Particle_size]/F");
  t_genParts_->Branch("Px",             ev.gl_px,       "Px[Particle_size]/F");
  t_genParts_->Branch("Py",             ev.gl_py,       "Py[Particle_size]/F");
  t_genParts_->Branch("Pz",             ev.gl_pz,       "Pz[Particle_size]/F");
  t_genParts_->Branch("E",              ev.gl_nrj,      "E[Particle_size]/F");
  t_genParts_->Branch("PT",             ev.gl_pt,       "PT[Particle_size]/F");
  t_genParts_->Branch("Eta",            ev.gl_eta,      "Eta[Particle_size]/F");
  t_genParts_->Branch("Phi",            ev.gl_phi,      "Phi[Particle_size]/F");
  t_genParts_->Branch("Mass",           ev.gl_mass,     "Mass[Particle_size]/F");
  t_genParts_->Branch("IsolationVar",   ev.gl_relIso,   "IsolationVar/F");

  t_genJets_->Branch("GenJet_size",     &ev.ngj,        "GenJet_size/I");
  t_genJets_->Branch("PT",              ev.gj_pt,       "PT[GenJet_size]/F");
  t_genJets_->Branch("Eta",             ev.gj_eta,      "Eta[GenJet_size]/F");
  t_genJets_->Branch("Phi",             ev.gj_phi,      "Phi[GenJet_size]/F");
  t_genJets_->Branch("Mass",            ev.gj_mass,     "Mass[GenJet_size]/F");

  t_genPhotons_->Branch("GenPhoton_size",  &ev.ngp,        "GenPhoton_size/I");
  t_genPhotons_->Branch("Status",          ev.gp_st,       "Status[GenPhoton_size]/I");
  t_genPhotons_->Branch("P",               ev.gp_p,        "P[GenPhoton_size]/F");
  t_genPhotons_->Branch("Px",              ev.gp_px,       "Px[GenPhoton_size]/F");
  t_genPhotons_->Branch("Py",              ev.gp_py,       "Py[GenPhoton_size]/F");
  t_genPhotons_->Branch("Pz",              ev.gp_pz,       "Pz[GenPhoton_size]/F");
  t_genPhotons_->Branch("E",               ev.gp_nrj,      "E[GenPhoton_size]/F");
  t_genPhotons_->Branch("PT",              ev.gp_pt,       "PT[GenPhoton_size]/F");
  t_genPhotons_->Branch("Eta",             ev.gp_eta,      "Eta[GenPhoton_size]/F");
  t_genPhotons_->Branch("Phi",             ev.gp_phi,      "Phi[GenPhoton_size]/F");
   
  //reco level event
  t_vertices_->Branch("Vertex_size",    &ev.nvtx,       "Vertex_size/I");
  t_vertices_->Branch("SumPT2",         &ev.v_pt2,      "SumPT2[Vertex_size]/F");

  t_looseElecs_->Branch("ElectronLoose_size", &ev.nle,  "ElectronLoose_size/I");
  t_looseElecs_->Branch("Charge",       ev.le_ch,       "Charge[ElectronLoose_size]/I");
  t_looseElecs_->Branch("Particle",     ev.le_g,        "Particle[ElectronLoose_size]/I");
  t_looseElecs_->Branch("PT",           ev.le_pt,       "PT[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Eta",          ev.le_eta,      "Eta[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Phi",          ev.le_phi,      "Phi[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Mass",         ev.le_mass,     "Mass[ElectronLoose_size]/F");
  t_looseElecs_->Branch("IsolationVar", ev.le_relIso,   "IsolationVar[ElectronLoose_size]/F");

  t_mediumElecs_->Branch("ElectronMedium_size", &ev.nme,  "ElectronMedium_size/I");
  t_mediumElecs_->Branch("Charge",               ev.me_ch,       "Charge[ElectronMedium_size]/I");
  t_mediumElecs_->Branch("Particle",             ev.me_g,        "Particle[ElectronMedium_size]/I");
  t_mediumElecs_->Branch("PT",                   ev.me_pt,       "PT[ElectronMedium_size]/F");
  t_mediumElecs_->Branch("Eta",                  ev.me_eta,      "Eta[ElectronMedium_size]/F");
  t_mediumElecs_->Branch("Phi",                  ev.me_phi,      "Phi[ElectronMedium_size]/F");
  t_mediumElecs_->Branch("Mass",                 ev.me_mass,     "Mass[ElectronMedium_size]/F");
  t_mediumElecs_->Branch("IsolationVar",         ev.me_relIso,   "IsolationVar[ElectronMedium_size]/F");

  t_tightElecs_->Branch("ElectronTight_size", &ev.nte,  "ElectronTight_size/I");
  t_tightElecs_->Branch("Charge",       ev.te_ch,       "Charge[ElectronTight_size]/I");
  t_tightElecs_->Branch("Particle",     ev.te_g,        "Particle[ElectronTight_size]/I");
  t_tightElecs_->Branch("PT",           ev.te_pt,       "PT[ElectronTight_size]/F");
  t_tightElecs_->Branch("Eta",          ev.te_eta,      "Eta[ElectronTight_size]/F");
  t_tightElecs_->Branch("Phi",          ev.te_phi,      "Phi[ElectronTight_size]/F");
  t_tightElecs_->Branch("Mass",         ev.te_mass,     "Mass[ElectronTight_size]/F");
  t_tightElecs_->Branch("IsolationVar", ev.te_relIso,   "IsolationVar[ElectronTight_size]/F");

  t_looseMuons_->Branch("MuonLoose_size", &ev.nlm,      "MuonLoose_size/I");
  t_looseMuons_->Branch("Charge",       ev.lm_ch,       "Charge[MuonLoose_size]/I");
  t_looseMuons_->Branch("Particle",     ev.lm_g,        "Particle[MuonLoose_size]/I");
  t_looseMuons_->Branch("PT",           ev.lm_pt,       "PT[MuonLoose_size]/F");
  t_looseMuons_->Branch("Eta",          ev.lm_eta,      "Eta[MuonLoose_size]/F");
  t_looseMuons_->Branch("Phi",          ev.lm_phi,      "Phi[MuonLoose_size]/F");
  t_looseMuons_->Branch("Mass",         ev.lm_mass,     "Mass[MuonLoose_size]/F");
  t_looseMuons_->Branch("IsolationVar", ev.lm_relIso,   "IsolationVar[MuonLoose_size]/F");

  t_tightMuons_->Branch("MuonTight_size", &ev.ntm,      "MuonTight_size/I");
  t_tightMuons_->Branch("Charge",       ev.tm_ch,       "Charge[MuonTight_size]/I");
  t_tightMuons_->Branch("Particle",     ev.tm_g,        "Particle[MuonTight_size]/I");
  t_tightMuons_->Branch("PT",           ev.tm_pt,       "PT[MuonTight_size]/F");
  t_tightMuons_->Branch("Eta",          ev.tm_eta,      "Eta[MuonTight_size]/F");
  t_tightMuons_->Branch("Phi",          ev.tm_phi,      "Phi[MuonTight_size]/F");
  t_tightMuons_->Branch("Mass",         ev.tm_mass,     "Mass[MuonTight_size]/F");
  t_tightMuons_->Branch("IsolationVar", ev.tm_relIso,   "IsolationVar[MuonTight_size]/F");

  t_puppiJets_->Branch("JetPUPPI_size", &ev.nj,         "JetPUPPI_size/I");
  t_puppiJets_->Branch("ID",            ev.j_id,        "ID[JetPUPPI_size]/I");
  t_puppiJets_->Branch("GenJet",        ev.j_g,         "GenJet[JetPUPPI_size]/I");
  t_puppiJets_->Branch("PT",            ev.j_pt,        "PT[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Eta",           ev.j_eta,       "Eta[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Phi",           ev.j_phi,       "Phi[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Mass",          ev.j_mass,      "Mass[JetPUPPI_size]/F");
  t_puppiJets_->Branch("MVAv2",         ev.j_mvav2,     "MVAv2[JetPUPPI_size]/I");
  t_puppiJets_->Branch("DeepCSV",       ev.j_deepcsv,   "DeepCSV[JetPUPPI_size]/I");
  t_puppiJets_->Branch("PartonFlavor",  ev.j_flav,      "PartonFlavor[JetPUPPI_size]/I");
  t_puppiJets_->Branch("HadronFlavor",  ev.j_hadflav,   "HadronFlavor[JetPUPPI_size]/I");
  t_puppiJets_->Branch("GenPartonPID",  ev.j_pid,       "GenPartonPID[JetPUPPI_size]/I");

  t_puppiMET_->Branch("PuppiMissingET_size", &ev.nmet,  "PuppiMissingET_size/I");
  t_puppiMET_->Branch("MET",            ev.met_pt,      "MET[PuppiMissingET_size]/F");
  t_puppiMET_->Branch("Phi",            ev.met_phi,     "Phi[PuppiMissingET_size]/F");
  t_puppiMET_->Branch("Eta",            ev.met_eta,     "Eta[PuppiMissingET_size]/F");

  t_loosePhotons_->Branch("PhotonLoose_size", &ev.nlp,     "PhotonLoose_size/I");
  t_loosePhotons_->Branch("Particle",     ev.lp_g,         "Particle[PhotonLoose_size]/I");
  t_loosePhotons_->Branch("IsEB",         ev.lp_isEB,      "IsEB[PhotonLoose_size]/I");
  t_loosePhotons_->Branch("PT",           ev.lp_pt,        "PT[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Eta",          ev.lp_eta,       "Eta[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Phi",          ev.lp_phi,       "Phi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("E",            ev.lp_nrj,       "E[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("PT_multi",     ev.lp_pt_multi,  "PT_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Eta_multi",    ev.lp_eta_multi, "Eta_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Phi_multi",    ev.lp_phi_multi, "Phi_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("E_multi",      ev.lp_nrj_multi, "E_multi[PhotonLoose_size]/F");

  t_tightPhotons_->Branch("PhotonTight_size", &ev.ntp,     "PhotonTight_size/I");
  t_tightPhotons_->Branch("Particle",     ev.tp_g,         "Particle[PhotonTight_size]/I");
  t_tightPhotons_->Branch("IsEB",         ev.tp_isEB,      "IsEB[PhotonTight_size]/I");
  t_tightPhotons_->Branch("PT",           ev.tp_pt,        "PT[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Eta",          ev.tp_eta,       "Eta[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Phi",          ev.tp_phi,       "Phi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("E",            ev.tp_nrj,       "E[PhotonTight_size]/F");
  t_tightPhotons_->Branch("PT_multi",     ev.tp_pt_multi,  "PT_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Eta_multi",    ev.tp_eta_multi, "Eta_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Phi_multi",    ev.tp_phi_multi, "Phi_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("E_multi",      ev.tp_nrj_multi, "E_multi[PhotonTight_size]/F");
}

void attachToMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_genPhotons_, TTree *t_looseElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_, TTree *t_loosePhotons_, TTree *t_tightPhotons_, MiniEvent_t &ev)
{
  //event header
  t_event_->SetBranchAddress("Run",               &ev.run);
  t_event_->SetBranchAddress("Event",             &ev.event);
  t_event_->SetBranchAddress("Lumi",              &ev.lumi);

  //gen level event
  t_event_->SetBranchAddress("Weight_size",       &ev.g_nw);
  t_event_->SetBranchAddress("Weight",            ev.g_w);

  t_genParts_->SetBranchAddress("Particle_size",  &ev.ngl);
  t_genParts_->SetBranchAddress("PID",            ev.gl_pid);
  t_genParts_->SetBranchAddress("Charge",         ev.gl_ch);
  t_genParts_->SetBranchAddress("Status",         ev.gl_st);
  t_genParts_->SetBranchAddress("P",              ev.gl_p);
  t_genParts_->SetBranchAddress("Px",             ev.gl_px);
  t_genParts_->SetBranchAddress("Py",             ev.gl_py);
  t_genParts_->SetBranchAddress("Pz",             ev.gl_pz);
  t_genParts_->SetBranchAddress("E",              ev.gl_nrj);
  t_genParts_->SetBranchAddress("PT",             ev.gl_pt);
  t_genParts_->SetBranchAddress("Eta",            ev.gl_eta);
  t_genParts_->SetBranchAddress("Phi",            ev.gl_phi);
  t_genParts_->SetBranchAddress("Mass",           ev.gl_mass);
  t_genParts_->SetBranchAddress("IsolationVar",   ev.gl_relIso);

  t_genJets_->SetBranchAddress("GenJet_size",     &ev.ngj);
  t_genJets_->SetBranchAddress("PT",              ev.gj_pt);
  t_genJets_->SetBranchAddress("Eta",             ev.gj_eta);
  t_genJets_->SetBranchAddress("Phi",             ev.gj_phi);
  t_genJets_->SetBranchAddress("Mass",            ev.gj_mass);

  t_genPhotons_->SetBranchAddress("GenPhoton_size",  &ev.ngp);
  t_genPhotons_->SetBranchAddress("Status",          ev.gp_st);
  t_genPhotons_->SetBranchAddress("P",               ev.gp_p);
  t_genPhotons_->SetBranchAddress("Px",              ev.gp_px);
  t_genPhotons_->SetBranchAddress("Py",              ev.gp_py);
  t_genPhotons_->SetBranchAddress("Pz",              ev.gp_pz);
  t_genPhotons_->SetBranchAddress("E",               ev.gp_nrj);
  t_genPhotons_->SetBranchAddress("PT",              ev.gp_pt);
  t_genPhotons_->SetBranchAddress("Eta",             ev.gp_eta);
  t_genPhotons_->SetBranchAddress("Phi",             ev.gp_phi);
   
  //reco level event
  t_vertices_->SetBranchAddress("Vertex_size",    &ev.nvtx);
  t_vertices_->SetBranchAddress("SumPT2",         &ev.v_pt2);

  t_looseElecs_->SetBranchAddress("ElectronLoose_size", &ev.nle);
  t_looseElecs_->SetBranchAddress("Charge",       ev.le_ch);
  t_looseElecs_->SetBranchAddress("Particle",     ev.le_g);
  t_looseElecs_->SetBranchAddress("PT",           ev.le_pt);
  t_looseElecs_->SetBranchAddress("Eta",          ev.le_eta);
  t_looseElecs_->SetBranchAddress("Phi",          ev.le_phi);
  t_looseElecs_->SetBranchAddress("Mass",         ev.le_mass);
  t_looseElecs_->SetBranchAddress("IsolationVar", ev.le_relIso);

  t_tightElecs_->SetBranchAddress("ElectronTight_size", &ev.nte);
  t_tightElecs_->SetBranchAddress("Charge",       ev.te_ch);
  t_tightElecs_->SetBranchAddress("Particle",     ev.te_g);
  t_tightElecs_->SetBranchAddress("PT",           ev.te_pt);
  t_tightElecs_->SetBranchAddress("Eta",          ev.te_eta);
  t_tightElecs_->SetBranchAddress("Phi",          ev.te_phi);
  t_tightElecs_->SetBranchAddress("Mass",         ev.te_mass);
  t_tightElecs_->SetBranchAddress("IsolationVar", ev.te_relIso);
  
  t_mediumElecs_->SetBranchAddress("ElectronMedium_size", &ev.nme);
  t_mediumElecs_->SetBranchAddress("Charge",               ev.me_ch);
  t_mediumElecs_->SetBranchAddress("Particle",             ev.me_g);
  t_mediumElecs_->SetBranchAddress("PT",                   ev.me_pt);
  t_mediumElecs_->SetBranchAddress("Eta",                  ev.me_eta);
  t_mediumElecs_->SetBranchAddress("Phi",                  ev.me_phi);
  t_mediumElecs_->SetBranchAddress("Mass",                 ev.me_mass);
  t_mediumElecs_->SetBranchAddress("IsolationVar",         ev.me_relIso);

  t_looseMuons_->SetBranchAddress("MuonLoose_size", &ev.nlm);
  t_looseMuons_->SetBranchAddress("Charge",       ev.lm_ch);
  t_looseMuons_->SetBranchAddress("Particle",     ev.lm_g);
  t_looseMuons_->SetBranchAddress("PT",           ev.lm_pt);
  t_looseMuons_->SetBranchAddress("Eta",          ev.lm_eta);
  t_looseMuons_->SetBranchAddress("Phi",          ev.lm_phi);
  t_looseMuons_->SetBranchAddress("Mass",         ev.lm_mass);
  t_looseMuons_->SetBranchAddress("IsolationVar", ev.lm_relIso);

  t_tightMuons_->SetBranchAddress("MuonTight_size", &ev.ntm);
  t_tightMuons_->SetBranchAddress("Charge",       ev.tm_ch);
  t_tightMuons_->SetBranchAddress("Particle",     ev.tm_g);
  t_tightMuons_->SetBranchAddress("PT",           ev.tm_pt);
  t_tightMuons_->SetBranchAddress("Eta",          ev.tm_eta);
  t_tightMuons_->SetBranchAddress("Phi",          ev.tm_phi);
  t_tightMuons_->SetBranchAddress("Mass",         ev.tm_mass);
  t_tightMuons_->SetBranchAddress("IsolationVar", ev.tm_relIso);

  t_puppiJets_->SetBranchAddress("JetPUPPI_size", &ev.nj);
  t_puppiJets_->SetBranchAddress("ID",            ev.j_id);
  t_puppiJets_->SetBranchAddress("GenJet",        ev.j_g);
  t_puppiJets_->SetBranchAddress("PT",            ev.j_pt);
  t_puppiJets_->SetBranchAddress("Eta",           ev.j_eta);
  t_puppiJets_->SetBranchAddress("Phi",           ev.j_phi);
  t_puppiJets_->SetBranchAddress("Mass",          ev.j_mass);
  t_puppiJets_->SetBranchAddress("MVAv2",         ev.j_mvav2);
  t_puppiJets_->SetBranchAddress("DeepCSV",       ev.j_deepcsv);
  t_puppiJets_->SetBranchAddress("PartonFlavor",  ev.j_flav);
  t_puppiJets_->SetBranchAddress("HadronFlavor",  ev.j_hadflav);
  t_puppiJets_->SetBranchAddress("GenPartonPID",  ev.j_pid);

  t_puppiMET_->SetBranchAddress("PuppiMissingET_size", &ev.nmet);
  t_puppiMET_->SetBranchAddress("MET",            ev.met_pt);
  t_puppiMET_->SetBranchAddress("Phi",            ev.met_phi);
  t_puppiMET_->SetBranchAddress("Eta",            ev.met_eta);

  t_loosePhotons_->SetBranchAddress("PhotonLoose_size", &ev.nlp);
  t_loosePhotons_->SetBranchAddress("Particle",     ev.lp_g);
  t_loosePhotons_->SetBranchAddress("IsEB",         ev.lp_isEB);
  t_loosePhotons_->SetBranchAddress("PT",           ev.lp_pt);
  t_loosePhotons_->SetBranchAddress("Eta",          ev.lp_eta);
  t_loosePhotons_->SetBranchAddress("Phi",          ev.lp_phi);
  t_loosePhotons_->SetBranchAddress("E",            ev.lp_nrj);
  t_loosePhotons_->SetBranchAddress("PT_multi",     ev.lp_pt_multi);
  t_loosePhotons_->SetBranchAddress("Eta_multi",    ev.lp_eta_multi);
  t_loosePhotons_->SetBranchAddress("Phi_multi",    ev.lp_phi_multi);
  t_loosePhotons_->SetBranchAddress("E_multi",      ev.lp_nrj_multi);

  t_tightPhotons_->SetBranchAddress("PhotonTight_size", &ev.ntp);
  t_tightPhotons_->SetBranchAddress("Particle",     ev.tp_g);
  t_tightPhotons_->SetBranchAddress("IsEB",         ev.tp_isEB);
  t_tightPhotons_->SetBranchAddress("PT",           ev.tp_pt);
  t_tightPhotons_->SetBranchAddress("Eta",          ev.tp_eta);
  t_tightPhotons_->SetBranchAddress("Phi",          ev.tp_phi);
  t_tightPhotons_->SetBranchAddress("E",            ev.tp_nrj);
  t_tightPhotons_->SetBranchAddress("PT_multi",     ev.tp_pt_multi);
  t_tightPhotons_->SetBranchAddress("Eta_multi",    ev.tp_eta_multi);
  t_tightPhotons_->SetBranchAddress("Phi_multi",    ev.tp_phi_multi);
  t_tightPhotons_->SetBranchAddress("E_multi",      ev.tp_nrj_multi);
}
