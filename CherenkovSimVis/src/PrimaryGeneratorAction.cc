// src/PrimaryGeneratorAction.cc - FIXED AND RANDOMIZED
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomTools.hh"
#include "G4PhysicalConstants.hh"
#include <cmath> // For std::sqrt, std::cos, std::sin

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr)
{
    fParticleGun = new G4ParticleGun(1);

    // Set particle type to mu- (this is done only once)
    auto particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
    fParticleGun->SetParticleDefinition(particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // --- Define the Generation Surface ---
    // A flat circular disk positioned just above the detector tank.
    // Tank radius is 16.9m, half-height is 18.1m.
    const G4double generation_radius = 16.9 * m;
    const G4double z_position = 20.0 * m; // Start 2m above the tank top

    // --- 1. Randomize Starting Position (uniformly on the disk) ---
    // We use sqrt(rand) to ensure uniform distribution over the area of the disk.
    G4double r = generation_radius * std::sqrt(G4UniformRand());
    G4double phi_pos = G4UniformRand() * 2. * CLHEP::pi;
    G4double x = r * std::cos(phi_pos);
    G4double y = r * std::sin(phi_pos);

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z_position));

    // --- 2. Randomize Direction (isotropic downward-going hemisphere) ---
    // This creates an isotropic flux pointing towards the detector.
    G4double cos_theta_dir = G4UniformRand(); // cos(theta) from 0 to 1
    G4double sin_theta_dir = std::sqrt(1. - cos_theta_dir * cos_theta_dir);
    G4double phi_dir = G4UniformRand() * 2. * CLHEP::pi;

    // Direction vector points towards -z
    G4double ux = sin_theta_dir * std::cos(phi_dir);
    G4double uy = sin_theta_dir * std::sin(phi_dir);
    G4double uz = -cos_theta_dir; // Always negative to point downwards

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

    // --- 3. Randomize Energy ---
    // A range of energies ensures muons can pass through the detector.
    // G4double energy = 2.0 * GeV + G4UniformRand() * 8.0 * GeV;
    // --- TEST: Forcing a fixed low energy to test the shower hypothesis ---
    const G4double energy = 500.0 * MeV;
    // --- END TEST --- 
    fParticleGun->SetParticleEnergy(energy);

    // --- Create the primary particle for this event ---
    fParticleGun->GeneratePrimaryVertex(anEvent);
}