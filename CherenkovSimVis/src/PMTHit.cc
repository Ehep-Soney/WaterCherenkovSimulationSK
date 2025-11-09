// PMTHit.cc
#include "PMTHit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

PMTHit::PMTHit()
  : G4VHit(),
    fPMTID(-1),
    fTotalPE(0),
    fPosition(G4ThreeVector())
{}

PMTHit::PMTHit(int pmtID)
  : G4VHit(),
    fPMTID(pmtID),
    fTotalPE(0),
    fPosition(G4ThreeVector())
{}

void PMTHit::Print() const {
    G4cout << "PMTHit: ID=" << fPMTID
           << " totalPE=" << fTotalPE
           << " pos=" << fPosition / m << " m";
    if (!fTimes.empty()) {
        G4cout << " times(ns):";
        for (auto t : fTimes) G4cout << " " << (t / ns);
    }
    G4cout << G4endl;
}