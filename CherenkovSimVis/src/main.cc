#include "G4RunManager.hh"
#include "EventAction.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"

int main(int argc, char** argv) {
// Create the run manager
auto* runManager = new G4RunManager();

// Register detector construction
runManager->SetUserInitialization(new DetectorConstruction());

// Register physics list with optical physics
auto* physicsList = new FTFP_BERT();
auto* opticalPhysics = new G4OpticalPhysics();

// Set verbosity for debugging
opticalPhysics->SetVerboseLevel(1);

physicsList->RegisterPhysics(opticalPhysics);
runManager->SetUserInitialization(physicsList);

// Register primary generator
runManager->SetUserAction(new PrimaryGeneratorAction());

// Register event action
runManager->SetUserAction(new EventAction());

// Initialize run manager BEFORE setting optical parameters
runManager->Initialize();

// Get UI manager for commands
G4UImanager* UImanager = G4UImanager::GetUIpointer();

// **Configure optical processes via UI commands (works for all Geant4 versions)**
UImanager->ApplyCommand("/process/optical/verbose 1");

// Cerenkov settings
UImanager->ApplyCommand("/process/optical/cerenkov/setMaxBetaChange 10.0");
UImanager->ApplyCommand("/process/optical/cerenkov/setTrackSecondariesFirst true");

// Boundary process settings (critical for PMT detection)
UImanager->ApplyCommand("/process/optical/boundary/setInvokeSD true");

// Scintillation (if needed)
UImanager->ApplyCommand("/process/optical/scintillation/setTrackSecondariesFirst false");

// Visualization
auto* visManager = new G4VisExecutive();
visManager->Initialize();

// Check if macro file provided
if (argc > 1) {
    // Batch mode with macro
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
} else {
    // Interactive mode
    auto* ui = new G4UIExecutive(argc, argv);
    
    // Optional: execute a visualization macro
    // UImanager->ApplyCommand("/control/execute vis.mac");
    
    ui->SessionStart();
    delete ui;
}

// Clean up
delete visManager;
delete runManager;

return 0;

}