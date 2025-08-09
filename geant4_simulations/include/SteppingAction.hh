// SteppingAction.hh
#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>
#include <G4ThreeVector.hh>
#include <G4GenericMessenger.hh>

class RunAction;

class SteppingAction : public G4UserSteppingAction {
  public:
    SteppingAction(RunAction* run);
    void UserSteppingAction(const G4Step*) override;

  private:
    RunAction* fRunAction;
    unsigned long long fStepID = 0;

    // --- new ---
    G4double            fKillE  = 1.0;//*CLHEP::GeV;   // default 1 GeV
    G4GenericMessenger* fMessenger;
};

#endif

