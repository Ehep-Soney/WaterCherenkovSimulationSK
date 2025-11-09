// EventAction.hh - NEW ROOT-ENABLED VERSION
#pragma once

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include <vector> // We need this to use std::vector

// Forward declaration of our data manager. This is a C++ trick that
// tells the compiler "this class exists" without needing its full header.
class DataManager;

class EventAction : public G4UserEventAction {
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // This accessor is still useful, so we keep it.
    int GetLastEventTotalPE() const { return fLastEventPE; }

private:
    int fLastEventPE;

    // --- NEW MEMBERS for ROOT data handling ---

    // A pointer to hold the single instance of our DataManager.
    DataManager* fDataManager;

    // A set of vectors to temporarily store all the hit data for ONE event.
    // At the end of the event, we will pass these vectors to the DataManager.
    std::vector<int>    fPmtID;
    std::vector<double> fPmtX;
    std::vector<double> fPmtY;
    std::vector<double> fPmtZ;
    std::vector<int>    fPmtNumPEs;
};