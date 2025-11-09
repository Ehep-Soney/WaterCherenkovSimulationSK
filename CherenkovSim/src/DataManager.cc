// src/DataManager.cc
#include "DataManager.hh"

// We must include the full ROOT headers here in the .cc file
#include "TFile.h"
#include "TTree.h"
#include "G4ios.hh" // Geant4's cout wrapper

// This is the implementation of the singleton pattern.
// We create a single static instance of the class.
DataManager* DataManager::GetInstance() {
    static DataManager instance;
    return &instance;
}

// The constructor is where we set up the TTree.
DataManager::DataManager() {
    // Create a new TTree. The first argument is the TTree's name in the file,
    // and the second is its title.
    fTree = new TTree("hits", "PMT Hit Data Tree");

    // Now, we connect our member variables to the "branches" of the TTree.
    // This tells ROOT: "When I call Fill(), take the data from these C++
    // variables and write it to the file."
    
    // The syntax is: Branch("BranchName", &variableAddress)
    fTree->Branch("eventID", &fEventID);
    fTree->Branch("pmtID", &fPmtID);
    fTree->Branch("pmtX", &fPmtX);
    fTree->Branch("pmtY", &fPmtY);
    fTree->Branch("pmtZ", &fPmtZ);
    fTree->Branch("numPEs", &fPmtNumPEs);
}

// The destructor is where we clean up memory.
DataManager::~DataManager() {
    // TFile owns the TTree, so deleting the file will also delete the tree.
    delete fRootFile;
}

// Opens the ROOT file for writing.
void DataManager::OpenFile(const std::string& filename) {
    // Check if a file is already open
    if (fRootFile) {
        G4cout << "DataManager::OpenFile - WARNING: A file is already open. Closing it." << G4endl;
        CloseFile();
    }
    
    // Create a new TFile. "RECREATE" means if the file already exists, overwrite it.
    fRootFile = new TFile(filename.c_str(), "RECREATE");
    if (!fRootFile || fRootFile->IsZombie()) {
        G4cerr << "DataManager::OpenFile - ERROR: Could not open file: " << filename << G4endl;
        fRootFile = nullptr;
    } else {
        G4cout << "DataManager: Output file opened: " << filename << G4endl;
        // Associate the TTree with the currently open file
        fTree->SetDirectory(fRootFile);
    }
}

// Writes the TTree to the file and closes it.
void DataManager::CloseFile() {
    if (fRootFile) {
        // Change to the file's context before writing
        fRootFile->cd();
        // The "Write" command saves the TTree header to the file.
        fTree->Write();
        // Close the file
        fRootFile->Close();
        G4cout << "DataManager: Output file closed." << G4endl;
        fRootFile = nullptr;
    }
}

// This function takes the hit data from EventAction and copies it into
// our member variables, then calls TTree::Fill().
void DataManager::FillEventData(int eventID, const std::vector<int>& pmtIDs, 
                               const std::vector<double>& pmtXs, const std::vector<double>& pmtYs, 
                               const std::vector<double>& pmtZs, const std::vector<int>& pmtPEs) {
    if (!fRootFile) {
        G4cerr << "DataManager::FillEventData - ERROR: No file is open." << G4endl;
        return;
    }

    // Copy the data from the input vectors to our member variables
    fEventID = eventID;
    fPmtID = pmtIDs;
    fPmtX = pmtXs;
    fPmtY = pmtYs;
    fPmtZ = pmtZs;
    fPmtNumPEs = pmtPEs;

    // This is the command that actually writes the data for one event to the TTree.
    fTree->Fill();
}

// This function is called at the start of each event to clear the std::vectors.
// This is crucial to prevent data from one event leaking into the next.
void DataManager::ClearEventData() {
    fEventID = 0;
    fPmtID.clear();
    fPmtX.clear();
    fPmtY.clear();
    fPmtZ.clear();
    fPmtNumPEs.clear();
}