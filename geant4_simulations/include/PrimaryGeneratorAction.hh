// PrimaryGeneratorAction.hh
#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4GenericMessenger.hh>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    void GeneratePrimaries(G4Event*) override;

    void SetEnergy(G4double e) { fEnergy = e; }
    G4double GetEnergy() const { return fEnergy; }

  private:
    G4ParticleGun* fGun;
    G4double fEnergy;            // GeV
    G4GenericMessenger* fMessenger;
};

#endif

