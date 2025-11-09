// src/ActionInitialization.cc

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const
{
    // In the simple single-threaded model, we create one of each action.
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());

    // We create an EventAction
    auto eventAction = new EventAction();
    SetUserAction(eventAction);

    // We create a SteppingAction and pass it the EventAction pointer.
    // We also store the pointer in our member variable.
    fSteppingAction = new SteppingAction(eventAction);
    SetUserAction(fSteppingAction);
}

SteppingAction* ActionInitialization::GetSteppingAction() const
{
    // This getter method simply returns the stored pointer.
    return fSteppingAction;
}