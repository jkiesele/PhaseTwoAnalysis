/*
 * LinkDef.h
 *
 *  Created on: 26 Aug 2016
 *      Author: jkiesele
 */

#ifndef CLASSDEF_LINKDEF_H_
#define CLASSDEF_LINKDEF_H_



#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector<Event>+;
#pragma link C++ class std::vector<LHCOEvent>+;
#pragma link C++ class std::vector<LHEFEvent>+;
#pragma link C++ class std::vector<LHEFWeight>+;
#pragma link C++ class std::vector<HepMCEvent>+;
#pragma link C++ class std::vector<GenParticle>+;
#pragma link C++ class std::vector<Vertex>+;
#pragma link C++ class std::vector<MissingET>+;
#pragma link C++ class std::vector<ScalarHT>+;
#pragma link C++ class std::vector<Rho>+;
#pragma link C++ class std::vector<Weight>+;
#pragma link C++ class std::vector<Photon>+;
#pragma link C++ class std::vector<Electron>+;
#pragma link C++ class std::vector<Muon>+;
#pragma link C++ class std::vector<Jet>+;
#pragma link C++ class std::vector<Track>+;
#pragma link C++ class std::vector<Tower>+;
#pragma link C++ class std::vector<HectorHit>+;
#pragma link C++ class std::vector<Candidate>+;

#endif


#endif /* CLASSDEF_LINKDEF_H_ */
