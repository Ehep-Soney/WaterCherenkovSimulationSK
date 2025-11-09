#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
PrimaryGeneratorAction();
~PrimaryGeneratorAction() override;

void GeneratePrimaries(G4Event* event) override;

private:
G4ParticleGun* fParticleGun = nullptr;
};
