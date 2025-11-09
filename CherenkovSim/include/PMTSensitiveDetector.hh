// PMTSensitiveDetector.hh
// No changes needed to header - all fixes are in the .cc implementation
#pragma once

#include "G4VSensitiveDetector.hh"
#include "PMTHit.hh"
#include "SimplePMTResponse.hh"
#include "G4THitsCollection.hh"
#include "G4TouchableHistory.hh"
#include <map>

class PMTSensitiveDetector : public G4VSensitiveDetector {
public:
    PMTSensitiveDetector(const G4String& name = "PMTSD",
                         const G4String& collectionName = "PMTHitsCollection");
    ~PMTSensitiveDetector() override = default;
    
    // Geant4 interface
    void Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* hist) override;
    void EndOfEvent(G4HCofThisEvent* hce) override;
    
    // Configure response
    void SetPMTResponse(const SimplePMTResponse& resp) { fResponse = resp; }

private:
    PMTHitsCollection* fHitsCollection;
    std::map<int, int> fHitMap; // pmtID -> index in hitsCollection
    int fHCID;
    SimplePMTResponse fResponse;
};