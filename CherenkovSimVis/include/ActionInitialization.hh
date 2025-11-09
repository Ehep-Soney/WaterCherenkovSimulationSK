// include/ActionInitialization.hh

#pragma once

#include "G4VUserActionInitialization.hh"

class SteppingAction; // Forward declaration

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization() = default;
    ~ActionInitialization() override = default;

    // We only need the Build() method for the simple, single-threaded design
    void Build() const override;

    // This getter method is how the EventAction will find the SteppingAction
    SteppingAction* GetSteppingAction() const;

private:
    // This variable will hold the single instance of our SteppingAction
    mutable SteppingAction* fSteppingAction = nullptr;
};