// RunAction.cc
#include "RunAction.hh"
#include <G4ios.hh>

RunAction::RunAction() {}
RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
    //char outfilename[128];
    //sprintf(outfilename, "muonSteps_%d_GeV.csv", (int)fEnergy);

    fOut.open("muonSteps.csv");
    fOut << "event,track,step,x(m),y(m),z(m),E(GeV)\n";
}

void RunAction::EndOfRunAction(const G4Run*)
{
    fOut.close();
}

