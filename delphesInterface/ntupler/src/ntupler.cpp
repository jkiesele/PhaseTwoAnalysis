/*
 * ntupler.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/ntupler.h"
#include "interface/scaleFactors.h"
//dirty hack
#include "../../../NTupler/src/MiniEvent.cc"
#include "TDirectory.h"
#include "TH1F.h"

void ntupler::analyze(size_t childid /* this info can be used for printouts */){

	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"JetPUPPI");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
	d_ana::dBranchHandler<Photon>      photon(tree(),"Photon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");
	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;


	//create output
	TString chilidstr="";
	chilidstr+=childid;
	TFile * outfile= new TFile(getOutDir()+"/p2ntuple_"+(TString)getLegendName()+"_"+chilidstr+".root","RECREATE");
	TDirectory *counterdir = outfile->mkdir("weightCounter");
	counterdir->cd();
	TH1F * h_event_weight = new TH1F("Event_weight","Event_weight",1,0,1);

	outfile->cd();
	TDirectory *ntupledir = outfile->mkdir("ntuple");
	ntupledir->cd();

	MiniEvent_t ev_;
	TTree * t_event_        = new TTree("Event","Event");
	TTree * t_genParts_     = new TTree("Particle","Particle");
	TTree * t_genPhotons_   = new TTree("GenPhoton","GenPhoton");
	TTree * t_vertices_     = new TTree("Vertex","Vertex");
	TTree * t_genJets_      = new TTree("GenJet","GenJet");
	TTree * t_looseElecs_   = new TTree("ElectronLoose","ElectronLoose");
	TTree * t_mediumElecs_  = new TTree("ElectronMedium","ElectronMedium");
	TTree * t_tightElecs_   = new TTree("ElectronTight","ElectronTight");
	TTree * t_looseMuons_   = new TTree("MuonLoose","MuonLoose");
	TTree * t_tightMuons_   = new TTree("MuonTight","MuonTight");
	TTree * t_puppiJets_    = new TTree("JetPUPPI","JetPUPPI");
	TTree * t_puppiMET_     = new TTree("PuppiMissingET","PuppiMissingET");
	TTree * t_loosePhotons_ = new TTree("PhotonLoose","PhotonLoose");
	TTree * t_tightPhotons_ = new TTree("PhotonTight","PhotonTight");
	createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_genPhotons_, t_looseElecs_,
			t_mediumElecs_,t_tightElecs_, t_looseMuons_, t_tightMuons_, t_puppiJets_, t_puppiMET_, t_loosePhotons_,
			t_tightPhotons_, ev_);


	//load effective corrections for delphes samples vs fullSim
	scaleFactors
	tightelecsf,medelecsf,looseelecsf,
	tightmuonsf,loosemuonsf,
	jetsf,
	tightphotonsf,loosephotonsf,
	metsf;

	TString basepath=getenv("CMSSW_BASE");
	basepath+="/src/PhaseTwoAnalysis/delphesInterface/ntupler/data/";

	tightelecsf.loadTH2D  (basepath+"ElectronTight_PTEta.root","FullSimOverDelphes");
	medelecsf.loadTH2D    (basepath+"ElectronMedium_PTEta.root","FullSimOverDelphes");
	//looseelecsf.loadTH2D  (cmsswbase+"bla.root","histo");
    //
	tightmuonsf.loadTH2D  (basepath+"MuonTight_PTEta.root","FullSimOverDelphes");
	//loosemuonsf.loadTH2D  (cmsswbase+"bla.root","histo");
    //
	//jetsf.loadTH2D        (cmsswbase+"bla.root","histo");
    //
	tightphotonsf.loadTH2D(basepath+"PhotonTight_PTEta.root","FullSimOverDelphes");
	//loosephotonsf.loadTH2D(cmsswbase+"bla.root","histo");
    //
	//metsf.loadTH2D        (cmsswbase+"bla.root","histo");

	for(size_t eventno=0;eventno<nevents;eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);

		if(event.size()){
			h_event_weight->Fill(0.,(double)event.at(0)->Weight);
		}
		else{
			h_event_weight->Fill(0.,1);
		}


		//skimming

		std::vector<Photon*>selectedphotons;
		for(size_t i=0;i<photon.size();i++){
			if(photon.at(i)->PT<20)continue;
			if(photon.at(i)->IsolationVarRhoCorr / photon.at(i)->E > 0.25)
				continue;
			selectedphotons.push_back(photon.at(i));
		}
		if(selectedphotons.size()<1)continue;
		std::vector<Electron*>selectedelectrons;
		for(size_t i=0;i<elecs.size();i++){
			if(elecs.at(i)->PT<15)continue;
			selectedelectrons.push_back(elecs.at(i));
		}
		if(muontight.size()+selectedelectrons.size()<1)continue;
		std::vector<Jet*>selectedjets;
		for(size_t i=0;i<jet.size();i++){
			if(jet.at(i)->PT<20)continue;
			selectedjets.push_back(jet.at(i));
		}
		if(selectedjets.size()<1)continue;


		//filling
		if(event.size()){
			ev_.event = event.at(0)->Number;
			ev_.g_nw = 1;
			ev_.g_w[0] = event.at(0)->Weight;
		}
		else{
			ev_.event = 0;
			ev_.g_nw = 1;
			ev_.g_w[0] = 1;
		}

		ev_.ntp=0;
		for(size_t i=0;i<selectedphotons.size();i++){
			if(ev_.ntp>=MiniEvent_t::maxpart)break;
			ev_.tp_eta[ev_.ntp]=selectedphotons.at(i)->Eta;
			ev_.tp_pt [ev_.ntp]=selectedphotons.at(i)->PT;
			ev_.tp_phi[ev_.ntp]=selectedphotons.at(i)->Phi;
			ev_.tp_nrj[ev_.ntp]=selectedphotons.at(i)->E;
			ev_.tp_sf[ev_.ntp]=tightphotonsf.getSF(selectedphotons.at(i)->Eta,selectedphotons.at(i)->PT);
			ev_.ntp++;
		}

		ev_.ntm=0;
		for(size_t i=0;i<muontight.size();i++){
			if(ev_.ntm>=MiniEvent_t::maxpart)break;
			ev_.tm_pt    [ev_.ntm] =muontight.at(i)->PT;
			ev_.tm_eta   [ev_.ntm]=muontight.at(i)->Eta;
			ev_.tm_phi   [ev_.ntm]=muontight.at(i)->Phi;
			ev_.tm_mass  [ev_.ntm]=0.105;
			ev_.tm_relIso[ev_.ntm]=muontight.at(i)->IsolationVarRhoCorr/muontight.at(i)->PT;
			ev_.tm_sf[ev_.ntm]=tightmuonsf.getSF(muontight.at(i)->Eta,muontight.at(i)->PT);
			ev_.ntm++;
		}
		ev_.nte=0;
		ev_.nme=0;
		for(size_t i=0;i<selectedelectrons.size();i++){
			if(ev_.nme>=MiniEvent_t::maxpart)break;

			ev_.me_pt    [ev_.nme] =selectedelectrons.at(i)->PT;
			ev_.me_eta   [ev_.nme]=selectedelectrons.at(i)->Eta;
			ev_.me_phi   [ev_.nme]=selectedelectrons.at(i)->Phi;
			ev_.me_mass  [ev_.nme]=0.00051;
			ev_.me_relIso[ev_.nme]=selectedelectrons.at(i)->IsolationVarRhoCorr /selectedelectrons.at(i)->PT ;
			ev_.me_sf[ev_.nme]=medelecsf.getSF(selectedelectrons.at(i)->Eta,selectedelectrons.at(i)->PT);
			ev_.nme++;

			ev_.te_pt    [ev_.nte] =selectedelectrons.at(i)->PT;
			ev_.te_eta   [ev_.nte]=selectedelectrons.at(i)->Eta;
			ev_.te_phi   [ev_.nte]=selectedelectrons.at(i)->Phi;
			ev_.te_mass  [ev_.nte]=0.00051;
			ev_.te_relIso[ev_.nte]=selectedelectrons.at(i)->IsolationVarRhoCorr /selectedelectrons.at(i)->PT ;
			ev_.te_sf[ev_.nte]=tightelecsf.getSF(selectedelectrons.at(i)->Eta,selectedelectrons.at(i)->PT);
			ev_.nte++;


		}
		ev_.nj=0;
		for(size_t i=0;i<selectedjets.size();i++){
			if(ev_.nj>=MiniEvent_t::maxjets)break;
			ev_.j_pt  [ev_.nj] =selectedjets.at(i)->PT;
			ev_.j_eta [ev_.nj]=selectedjets.at(i)->Eta;
			ev_.j_phi [ev_.nj]=selectedjets.at(i)->Phi;
			ev_.j_mass[ev_.nj]=selectedjets.at(i)->Mass;

			ev_.j_hadflav[ev_.nj]=selectedjets.at(i)->Flavor;

			ev_.j_deepcsv[ev_.nj]=0;
			ev_.j_mvav2[ev_.nj]=0;
			if(selectedjets.at(i)->BTag){
				ev_.j_deepcsv[ev_.nj]=0b00000111;
				ev_.j_mvav2[ev_.nj]=0b00000111;
			}
			ev_.j_sf[ev_.nj]=jetsf.getSF(selectedjets.at(i)->Eta,selectedjets.at(i)->PT);
			ev_.nj++;
		}

		ev_.nmet=0;
		for(size_t i=0;i<met.size();i++){
			if(ev_.nmet>=MiniEvent_t::maxpart) break;
			ev_.met_eta[ev_.nmet]=met.at(i)->Eta ;
			ev_.met_pt [ev_.nmet]=met.at(i)->MET ;
			ev_.met_phi[ev_.nmet]=met.at(i)->Phi ;
			ev_.met_sf[ev_.nmet]=metsf.getSF(0,selectedjets.at(i)->PT);
			ev_.nmet++;
		}


		t_event_->Fill();
		t_genParts_->Fill();
		t_genPhotons_->Fill();
		t_vertices_->Fill();
		t_genJets_->Fill();
		t_looseElecs_->Fill();
		t_tightElecs_->Fill();
		t_looseMuons_->Fill();
		t_tightMuons_->Fill();
		t_puppiJets_->Fill();
		t_puppiMET_->Fill();
		t_loosePhotons_->Fill();
		t_tightPhotons_->Fill();

	}

	counterdir->cd();
	h_event_weight->Write();

	ntupledir->cd();
	t_event_        ->Write();
	t_genParts_     ->Write();
	t_genPhotons_   ->Write();
	t_vertices_     ->Write();
	t_genJets_      ->Write();
	t_looseElecs_   ->Write();
	t_mediumElecs_   ->Write();
	t_tightElecs_   ->Write();
	t_looseMuons_   ->Write();
	t_tightMuons_   ->Write();
	t_puppiJets_    ->Write();
	t_puppiMET_     ->Write();
	t_loosePhotons_ ->Write();
	t_tightPhotons_ ->Write();

	outfile->Close();
	/*
	 * Must be called in the end, takes care of thread-safe writeout and
	 * call-back to the parent process
	 */
	processEndFunction();
}



void ntupler::postProcess(){

	/* empty */

}



