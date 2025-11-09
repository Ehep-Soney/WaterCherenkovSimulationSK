// src/RunAction.cc - NEW VERSION
#include "RunAction.hh"
#include "G4Run.hh"
#include "G4ios.hh"
#include <fstream>

RunAction::RunAction() : G4UserRunAction()
{
    // This message confirms that our RunAction was created successfully.
    G4cout << "RunAction constructor is called." << G4endl;
}

RunAction::~RunAction()
{
    // This is important for preventing memory leaks.
    // If we created a file object, we must delete it.
    delete fOutputFile;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    // Create a new file for this run.
    // This will overwrite any existing file with the same name.
    fOutputFile = new std::ofstream("muon_data.csv");

    // Check if the file was opened successfully
    if (fOutputFile && fOutputFile->is_open()) {
        // Write the header row for our CSV file.
        // This makes it easy to read with Python/pandas later.
        (*fOutputFile) << "path_length_cm,total_pe,cos_theta\n";
    } else {
        G4cerr << "Error: Could not open output file muon_data.csv" << G4endl;
    }
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
    // Check if the file object exists and is open before trying to close it
    if (fOutputFile && fOutputFile->is_open()) {
        fOutputFile->close();
        G4cout << "Output data saved to muon_data.csv" << G4endl;
    }

    G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;
}