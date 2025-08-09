// SteppingAction.cc
#include "SteppingAction.hh"
#include "RunAction.hh"
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4SystemOfUnits.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>

#include "G4RunManager.hh"          // header you need


SteppingAction::SteppingAction(RunAction* run) : fRunAction(run) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const auto track = step->GetTrack();
    const auto pdg   = track->GetParticleDefinition()->GetPDGEncoding();

    auto event   = G4RunManager::GetRunManager()->GetCurrentEvent();
    G4int eventID  = event ? event->GetEventID() : -1;   // ‑1 = safety fallback
    G4int  trackID   = track->GetTrackID();          // 1 = primary muon
    G4int  stepInTrk = track->GetCurrentStepNumber();


    if(std::abs(pdg) == 13) { // only mu±

        if (stepInTrk % 100 == 0 || stepInTrk==1) {

        const auto pos = step->GetPostStepPoint()->GetPosition();
        const auto E   = step->GetPostStepPoint()->GetKineticEnergy();

        auto& out = fRunAction->GetStream();
        //out << track->GetEventID() << ','
        //out << track->GetParentID() << ','
        out << eventID << ','
            //<< fStepID++ << ','
            << trackID << ','
            << stepInTrk << ','
            << pos.x()/m << ','
            << pos.y()/m << ','
            << pos.z()/m << ','
            << E/GeV << '\n';

        // --- NEW: kill if below threshold ---
        if(E/GeV < fKillE) {
            track->SetTrackStatus(fStopAndKill);
            fStepID = 0;
        }
      }
    }
    else {
            // GEANT4 will track all the secondary particles produced so we should 
            // kill these otherwise it takes forever.
            track->SetTrackStatus(fStopAndKill);
    }
}

