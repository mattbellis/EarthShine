// PrimaryGeneratorAction.cc
#include "PrimaryGeneratorAction.hh"
#include <G4ParticleTable.hh>
#include <G4MuonMinus.hh>
#include <G4SystemOfUnits.hh>

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : fEnergy(100.*GeV)
{
    fGun = new G4ParticleGun(1);
    fGun->SetParticleDefinition(G4MuonMinus::Definition());
    fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));

    // UI: /mygen/energy <value> <unit>
    fMessenger = new G4GenericMessenger(this,"/mygen/","Beam control");
    fMessenger->DeclarePropertyWithUnit("energy","GeV",fEnergy,
        "Muon kinetic energy");
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
    // Set energy & starting position at z = +1 m above rock top
    fGun->SetParticleEnergy(fEnergy);
    fGun->SetParticlePosition(G4ThreeVector(0,0,1.*m));
    fGun->GeneratePrimaryVertex(evt);
}

