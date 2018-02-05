/*
 * ntupler.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/minivalidator.h"
#include "interface/scaleFactors.h"
//dirty hack
#include "TDirectory.h"
#include "TH1F.h"
#include "TLorentzVector.h"

void minivalidator::analyze(size_t childid /* this info can be used for printouts */){

        d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
        d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
        d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
        d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
        d_ana::dBranchHandler<Jet>         jet(tree(),"JetPUPPI");
        d_ana::dBranchHandler<Jet>         taujet(tree(),"Jet");
        d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
        d_ana::dBranchHandler<Photon>      photon(tree(),"Photon");
        d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");
        d_ana::dBranchHandler<MissingET>   puppimet(tree(),"PuppiMissingET");
        d_ana::dBranchHandler<MissingET>   genpumet(tree(),"GenPileUpMissingET");
        d_ana::dBranchHandler<MissingET>   genmet(tree(),"GenMissingET");


        size_t nevents=tree()->entries();
        if(isTestMode())
                nevents/=100;


        //create output
        TString chilidstr="";
        chilidstr+=childid;
        TFile * outfile= new TFile(getOutDir()+"/p2val_"+(TString)getLegendName()+"_"+chilidstr+".root","RECREATE");
        TDirectory *counterdir = outfile->mkdir("weightCounter");
        counterdir->cd();
        TH1F * h_event_weight = new TH1F("Event_weight","Event_weight",1,0,1);

        outfile->cd();
        TDirectory *muondir = outfile->mkdir("muonTight");
        muondir->cd();
        TH1F * h_muontight_all_pt = new TH1F("h_muontight_all_pt","Muon PT",200,0,200);
        TH1F * h_muontight_all_eta = new TH1F("h_muontight_all_eta","Muon eta",200,-5,5);
        TH1F * h_muontight_all_phi = new TH1F("h_muontight_all_phi","Muon phi",200,-3.5,3.5);
        TH1F * h_muontight_all_tof = new TH1F("h_muontight_all_tof","Muon tof",200,0,10);
        TH1F * h_muontight_all_charge = new TH1F("h_muontight_all_charge","Muon charge",2,-1,1);
        TH1F * h_muontight_all_IsolationVar = new TH1F("h_muontight_all_IsolationVar","Muon IsolationVar",100,0,20);
        TH1F * h_muontight_all_IsolationVarRhoCorr = new TH1F("h_muontight_all_IsolationVarRhoCorr","Muon IsolationVarRhoCorr",100,0,20);
        TH1F * h_muontight_all_SumPtCharged = new TH1F("h_muontight_all_SumPtCharged","Muon SumPtCharged",100,0,20);
        TH1F * h_muontight_all_SumPtNeutral = new TH1F("h_muontight_all_SumPtNeutral","Muon SumPtNeutral",100,0,20);
        TH1F * h_muontight_all_SumPtChargedPU = new TH1F("h_muontight_all_SumPtChargedPU","Muon SumPtChargedPU",100,0,20);
        TH1F * h_muontight_all_SumPt = new TH1F("h_muontight_all_SumPt","Muon SumPt",100,0,20);

        TH1F * h_muontight_iso_pt = new TH1F("h_muontight_iso_pt","Muon PT",200,0,200);
        TH1F * h_muontight_iso_eta = new TH1F("h_muontight_iso_eta","Muon eta",200,-5,5);
        TH1F * h_muontight_iso_phi = new TH1F("h_muontight_iso_phi","Muon phi",200,-3.5,3.5);

        TH2F * h_muontight_all_IsolationVarRhoCorr_pt = new TH2F ("h_muontight_all_IsolationVarRhoCorr_pt","",200,0,20,200,0,20);

        outfile->cd();
        TDirectory *zmmdir = outfile->mkdir("zmm");
        zmmdir->cd();
        TH1F * h_zmm_m1_pt = new TH1F("h_zmm_m1_pt","Muon PT 1",200,0,200);
        TH1F * h_zmm_m2_pt = new TH1F("h_zmm_m2_pt","Muon PT 2",200,0,200);
        TH1F * h_zmm_mass  = new TH1F("h_zmm_mass","ZMM mass",200,60,120);
        TH1F * h_zmm_met_pt = new TH1F("h_zmm_met_pt","MET ZMM sel",200,0,200);
        TH1F * h_zmm_puppimet_pt = new TH1F("h_zmm_puppimet_pt","MET no sel",200,0,200);
        TH1F * h_zmm_genpumet_pt = new TH1F("h_zmm_genpumet_pt","MET no sel",200,0,200);
        TH1F * h_zmm_genmet_pt = new TH1F("h_zmm_genmet_pt","MET no sel",200,0,200);

        TDirectory *metdir = outfile->mkdir("met");
        metdir->cd();
        TH1F * h_met_pt = new TH1F("h_met_pt","MET no sel",200,0,200);
        TH1F * h_puppimet_pt = new TH1F("h_puppimet_pt","MET no sel",200,0,200);
        TH1F * h_genpumet_pt = new TH1F("h_genpumet_pt","MET no sel",200,0,200);
        TH1F * h_genmet_pt = new TH1F("h_genmet_pt","MET no sel",200,0,200);

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

                if(event.size()<1)continue;

                h_event_weight->Fill(0.,(double)event.at(0)->Weight);

                // No selection
                h_met_pt->Fill(met.at(0)->MET);
                h_puppimet_pt->Fill(puppimet.at(0)->MET);
                h_genpumet_pt->Fill(genpumet.at(0)->MET);
                h_genmet_pt->Fill(genmet.at(0)->MET);

                // Playing at genlevel



                // Lets look for muons  
                Double_t maxpt1=0, maxpt2=0;
                size_t index1=0, index2=0;
                size_t ngoodmuons=0;

                for(size_t i=0;i<muontight.size();i++){
                        h_muontight_all_pt->Fill(muontight.at(i)->PT);
                        h_muontight_all_eta->Fill(muontight.at(i)->Eta);
                        h_muontight_all_phi->Fill(muontight.at(i)->Phi);
                        h_muontight_all_charge->Fill(muontight.at(i)->Charge);
                        h_muontight_all_tof->Fill(muontight.at(i)->T);
                        h_muontight_all_IsolationVar->Fill(muontight.at(i)->IsolationVar);
                        h_muontight_all_IsolationVarRhoCorr->Fill(muontight.at(i)->IsolationVarRhoCorr);
                        h_muontight_all_SumPtCharged->Fill(muontight.at(i)->SumPtCharged);
                        h_muontight_all_SumPtNeutral->Fill(muontight.at(i)->SumPtNeutral);
                        h_muontight_all_SumPtChargedPU->Fill(muontight.at(i)->SumPtChargedPU);
                        h_muontight_all_SumPt->Fill(muontight.at(i)->SumPt);
                  
                        h_muontight_all_IsolationVarRhoCorr_pt ->Fill(muontight.at(i)->IsolationVar,muontight.at(i)->PT);

                        if(muontight.at(i)->IsolationVarRhoCorr<0.1) { // this is just a guess  
                                h_muontight_iso_pt->Fill(muontight.at(i)->PT);
                                h_muontight_iso_eta->Fill(muontight.at(i)->Eta);
                                h_muontight_iso_phi->Fill(muontight.at(i)->Phi);                

                                if(muontight.at(i)->PT>maxpt1) { maxpt2=maxpt1; index2=index1; maxpt1=muontight.at(i)->PT; index1=i;}
                                else if (muontight.at(i)->PT>maxpt2) {maxpt2=muontight.at(i)->PT; index2=i;} 
                                ngoodmuons++;
                        }
                }
                if(ngoodmuons>=2) {
                        if( muontight.at(index1)->PT>20 && muontight.at(index2)->PT > 20 && 
                                        (muontight.at(index1)->Charge!=muontight.at(index2)->Charge) ) { 
                        h_zmm_m1_pt->Fill(muontight.at(index1)->PT);
                        h_zmm_m2_pt->Fill(muontight.at(index2)->PT);

                        TLorentzVector m1, m2;
                        m1.SetPtEtaPhiM(muontight.at(index1)->PT,muontight.at(index1)->Eta,muontight.at(index1)->Phi,0.105); 
                        m2.SetPtEtaPhiM(muontight.at(index2)->PT,muontight.at(index2)->Eta,muontight.at(index2)->Phi,0.105);

                        h_zmm_mass->Fill( (m1+m2).M() );

                        h_zmm_met_pt->Fill(met.at(0)->MET);
                        h_zmm_puppimet_pt->Fill(puppimet.at(0)->MET);
                        h_zmm_genpumet_pt->Fill(genpumet.at(0)->MET);
                        h_zmm_genmet_pt->Fill(genmet.at(0)->MET);


//                        Double_t mass2=2*muontight.at(index1)->PT*muontight.at(index2)->PT*( cosh(muontight.at(index1)->Eta-muontight.at(index2)->Eta) - cos(muontight.at(index1)->Phi-muontight.at(index2)->Phi) );
//                        std::cout<<muontight.at(index1)->PT<<"   "<<m1.Pt()<<std::endl;
//                        std::cout<<sqrt(mass2)<<"   "<<(m1+m2).M()<<std::endl;

                  }
                }
        }   



counterdir->cd();
h_event_weight->Write();

muondir->cd();
h_muontight_all_pt->Write();
h_muontight_all_eta    ->Write();
h_muontight_all_phi    ->Write();
h_muontight_all_charge ->Write();
h_muontight_all_tof    ->Write();
h_muontight_all_IsolationVar        ->Write();
h_muontight_all_IsolationVarRhoCorr ->Write();
h_muontight_all_SumPtCharged        ->Write();
h_muontight_all_SumPtNeutral        ->Write();
h_muontight_all_SumPtChargedPU      ->Write();
h_muontight_all_SumPt               ->Write();

h_muontight_all_IsolationVarRhoCorr_pt->Write();

h_muontight_iso_pt->Write();
h_muontight_iso_eta    ->Write();
h_muontight_iso_phi    ->Write();

zmmdir->cd();
h_zmm_m1_pt->Write();
h_zmm_m2_pt->Write();
h_zmm_mass->Write();
h_zmm_met_pt->Write();
h_zmm_puppimet_pt->Write();
h_zmm_genpumet_pt->Write();
h_zmm_genmet_pt->Write();

metdir->cd();
h_met_pt->Write();
h_puppimet_pt->Write();
h_genpumet_pt->Write();
h_genmet_pt->Write();

outfile->Close();
/*
 * Must be called in the end, takes care of thread-safe writeout and
 * call-back to the parent process
 */
processEndFunction();
}



void minivalidator::postProcess(){

        /* empty */

}



