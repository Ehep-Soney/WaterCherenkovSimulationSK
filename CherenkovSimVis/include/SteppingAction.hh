// include/SteppingAction.hh
#pragma once

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class EventAction; // Forward declaration

class SteppingAction : public G4UserSteppingAction {
public:
// Constructor needs to know about the EventAction to pass it data
SteppingAction(EventAction* eventAction);
~SteppingAction() override = default;

void UserSteppingAction(const G4Step* step) override;

// Called by EventAction at the beginning of an event to reset
void Reset();

// Getters for EventAction to retrieve the results
const G4ThreeVector& GetEntryPoint() const { return fEntryPoint; }
const G4ThreeVector& GetExitPoint() const { return fExitPoint; }
bool FoundTrackInWater() const { return fFoundTrack; }

private:
EventAction* fEventAction;
G4ThreeVector fEntryPoint;
G4ThreeVector fExitPoint;
bool fFoundTrack;
};