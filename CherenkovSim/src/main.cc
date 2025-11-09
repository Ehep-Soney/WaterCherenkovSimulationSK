// main.cc - FINAL ROOT-ENABLED VERSION
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "DataManager.hh" // <-- Include our new DataManager

// You may need to include TrackingAction.hh if you are using it
// #include "TrackingAction.hh"

int main(int argc, char** argv) {
    auto* runManager = new G4RunManager();
    runManager->SetUserInitialization(new DetectorConstruction());
    
    auto* physicsList = new FTFP_BERT();
    physicsList->RegisterPhysics(new G4OpticalPhysics());
    runManager->SetUserInitialization(physicsList);
    
    runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserAction(new EventAction());
    // runManager->SetUserAction(new TrackingAction()); // Uncomment if you use it
    runManager->Initialize();

    // Get the DataManager instance
    DataManager* dataManager = DataManager::GetInstance();
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (argc > 1) {
        // --- BATCH MODE ---
        // Open the ROOT file, run the macro, and close the file.
        
        dataManager->OpenFile("pmt_hits.root"); // <-- Open the ROOT file

        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
        
        dataManager->CloseFile(); // <-- Close the ROOT file

        G4cout << "----> Batch run finished. Data saved to pmt_hits.root <----" << G4endl;

    } else {
        // --- INTERACTIVE MODE ---
        // (No data will be saved in this mode)
        auto* visManager = new G4VisExecutive();
        visManager->Initialize();
        auto* ui = new G4UIExecutive(argc, argv);
        
        // Use a compatible visualization macro
        UImanager->ApplyCommand("/control/execute vis_compatible.mac");
        
        ui->SessionStart();
        delete ui;
        delete visManager;
    }

    delete runManager;
    return 0;
}