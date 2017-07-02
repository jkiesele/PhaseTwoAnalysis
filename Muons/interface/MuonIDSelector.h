#ifndef MuonIDSelector_h
#define MuonIDSelector_h
/* \class MuonIDSelector
 *
 */
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

struct MuonIDSelector {
  MuonIDSelector( bool tightSelect ) : 
    tightSelect_( tightSelect ) { }
  bool operator()(reco::Muon muon, reco::Vertex vtx) const { 
    bool isLoose  = (fabs(muon.eta()) < 2.4 && muon::isLooseMuon(muon)) || (fabs(muon.eta()) > 2.4 && isME0MuonSel(muon, 3, 4, 3, 4, 0.5));
    // bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.3));
    bool isTight  = (fabs(muon.eta()) < 2.4 && muon::isTightMuon(muon,vtx) || (fabs(muon.eta()) > 2.4 && isME0MuonSel(muon, 3, 4, 3, 4, 0.1));

    return ( tightSelct_ ? isTight : isLoose); 

  }

private:
  bool tightSelect_;
  
  bool isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = fabs(chamber->x - segment->x);
          deltaY   = fabs(chamber->y - segment->y);
          pullX    = fabs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = fabs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = fabs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

};

#endif
