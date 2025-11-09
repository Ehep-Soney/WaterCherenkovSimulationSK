// src/PrimaryGeneratorAction.cc - "LIGHT BULB" POINT SOURCE Version
#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4OpticalPhoton.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(new G4ParticleGun(1)) {
    fParticleGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fParticleGun; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // === EDIT THIS LINE BETWEEN RUNS ===
    const G4double z_source = 16.0 * m; // Try 0m, 8m, 12m, 16m etc.
    // ===================================
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, z_source));

    const G4int num_photons_per_event = 50000;
    for (G4int i = 0; i < num_photons_per_event; ++i) {
        fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
        G4double angle = G4UniformRand() * 360.0 * deg;
        fParticleGun->SetParticlePolarization(G4ThreeVector(std::cos(angle), std::sin(angle), 0));
        fParticleGun->SetParticleEnergy(3.0 * eV); // Blue light
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}