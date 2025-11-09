// include/RunAction.hh - MODIFIED VERSION
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream> // Required for std::ofstream

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    // Getter method for the EventAction to access the output file stream
    std::ofstream* GetOutputFile() const { return fOutputFile; }

private:
    // A pointer to our output file stream object
    std::ofstream* fOutputFile = nullptr;
};