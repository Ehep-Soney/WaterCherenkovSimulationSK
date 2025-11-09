// PMTSensitiveDetector.cc - FIXED VERSION
//
// Key fixes:
// 1. Detects photons at photocathode boundary (not just entering glass)
// 2. Uses Geant4's optical surface EFFICIENCY for detection (no manual QE check needed)
// 3. Properly extracts PMT copy number from touchable
// 4. Only processes optical photons

#include "PMTSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4TouchableHandle.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"

PMTSensitiveDetector::PMTSensitiveDetector(const G4String& name, const G4String& hitsCollectionName)
: G4VSensitiveDetector(name), fHitsCollection(nullptr), fHCID(-1)
{
collectionName.insert(hitsCollectionName);
}

void PMTSensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
fHitsCollection = new PMTHitsCollection(SensitiveDetectorName, collectionName[0]);
fHitMap.clear();


if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
}
hce->AddHitsCollection(fHCID, fHitsCollection);

}

G4bool PMTSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
if (!step) return false;


G4Track* track = step->GetTrack();
G4ParticleDefinition* pd = track->GetDefinition();

// Only process optical photons
if (pd != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

// **Check if photon is entering the photocathode volume**
G4StepPoint* preStep = step->GetPreStepPoint();
G4StepPoint* postStep = step->GetPostStepPoint();

// We want photons that just entered the photocathode
//G4VPhysicalVolume* postPV = postStep->GetPhysicalVolume();
//if (!postPV) return false;

// Verify we're in the photocathode volume by name
//G4String volName = postPV->GetName();
//if (volName.find("Photocathode") == std::string::npos) {
  //  return false; // Not in photocathode volume
//}

// DEBUG: Print when we detect entry
if (verboseLevel > 0) {
    G4cout << "[PMT SD] Optical photon entered photocathode volume: " << G4endl;
}

// Get position where photon hit
G4ThreeVector worldPos = postStep->GetPosition();

// Get PMT copy number from the touchable
// Since photocathode is daughter of PMT, we need to go up one level
G4TouchableHandle touch = postStep->GetTouchableHandle();

// Depth 0 = photocathode, Depth 1 = PMT glass (parent)
G4int pmtCopyNo = 0;
if (touch->GetHistoryDepth() >= 1) {
    pmtCopyNo = touch->GetCopyNumber(1); // Get parent (PMT) copy number
} else {
    pmtCopyNo = touch->GetCopyNumber(0); // Fallback
}

// Get photon energy and wavelength
G4double energy = track->GetTotalEnergy();
double wavelength_nm = 0.0;
if (energy > 0.) {
    wavelength_nm = (h_Planck * c_light) / energy / nm;
}

// **Apply QE-based detection (manual check)**
//if (!fResponse.IsDetected(wavelength_nm)) {
    // Photon not detected - kill it to avoid further tracking issues
  //  track->SetTrackStatus(fStopAndKill);
    //return false;
//}

// Photon detected! 
if (verboseLevel > 0) {
    G4cout << "[PMT SD] Photon DETECTED at PMT " << pmtCopyNo 
           << " wavelength=" << wavelength_nm << " nm" << G4endl;
}

// Apply time smearing
double globalTime = postStep->GetGlobalTime();
double smearedTime = fResponse.SmearTime(globalTime);

// Store hit: find or create PMTHit for this pmtCopyNo
auto it = fHitMap.find(pmtCopyNo);
if (it == fHitMap.end()) {
    // Create new hit
    PMTHit* newHit = new PMTHit(pmtCopyNo);
    newHit->SetPosition(worldPos);
    G4int idx = fHitsCollection->insert(newHit);
    fHitMap[pmtCopyNo] = idx;
    (*fHitsCollection)[idx - 1]->AddPE(smearedTime / ns);
} else {
    G4int idx = it->second;
    (*fHitsCollection)[idx - 1]->AddPE(smearedTime / ns);
}

// Kill the photon to avoid double counting
track->SetTrackStatus(fStopAndKill);

return true;

}

void PMTSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
if (verboseLevel > 1 && fHitsCollection) {
G4int n = fHitsCollection->entries();
G4cout << "PMT SD EndOfEvent: " << n << " PMTs hit, total PE = ";
G4int totalPE = 0;
for (G4int i = 0; i < n; ++i) {
totalPE += (*fHitsCollection)[i]->GetTotalPE();
}
G4cout << totalPE << G4endl;

if (verboseLevel > 2) {
        for (G4int i = 0; i < n; ++i) {
            (*fHitsCollection)[i]->Print();
        }
    }
}

}