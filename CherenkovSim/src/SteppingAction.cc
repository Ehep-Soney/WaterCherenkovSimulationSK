// src/SteppingAction.cc
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{
Reset();
}

void SteppingAction::Reset() {
fFoundTrack = false;
fEntryPoint.set(0, 0, 0);
fExitPoint.set(0, 0, 0);
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// We only care about the primary muon (Track ID 1)
if (step->GetTrack()->GetTrackID() != 1) {
return;
}


// --- SAFETY CHECK ---
// Get the physical volume of the current step's endpoint
G4VPhysicalVolume* postStepVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
// If the pointer is null, it means the particle is leaving the world.
// We must stop here to prevent a crash.
if (!postStepVolume) {
return;
    }

// --- END SAFETY CHECK ---//

// Get the logical volume that the particle is *currently in* (at the post-step point)
G4LogicalVolume* currentVolume = postStepVolume->GetLogicalVolume();

// We only care about steps that are inside our water volume
if (currentVolume->GetName() != "InnerDetectorLogical") {
    return;
}

// If we are here, the muon is inside the water.
if (!fFoundTrack) {
    // This is the first step inside the water.
    // The "pre-step point" of this step is the exact entry point.
    fEntryPoint = step->GetPreStepPoint()->GetPosition();
    fFoundTrack = true;
}

// We are inside the water, so we continuously update the exit point.
// The last "post-step point" we record will be the final exit point.
fExitPoint = step->GetPostStepPoint()->GetPosition();

}