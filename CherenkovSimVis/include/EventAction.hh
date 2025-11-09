// EventAction.hh
#pragma once

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class EventAction : public G4UserEventAction {
public:
EventAction();
~EventAction() override = default;

void BeginOfEventAction(const G4Event* event) override;
void EndOfEventAction(const G4Event* event) override;

// simple accessor for last event's total PE
int GetLastEventTotalPE() const { return fLastEventPE; }

private:
int fLastEventPE;
};