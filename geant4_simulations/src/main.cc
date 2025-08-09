#include <G4RunManagerFactory.hh>
#include <FTFP_BERT.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

int main(int argc,char** argv)
{
    auto ui = (argc==1) ? new G4UIExecutive(argc,argv) : nullptr;

    // Run manager (single-threaded for simplicity)
    auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // User classes
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT);
    runManager->SetUserAction(new PrimaryGeneratorAction());

    auto runAction = new RunAction();
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(new SteppingAction(runAction));

    runManager->Initialize();

    // Visualization (optional)
    auto vis = new G4VisExecutive;
    vis->Initialize();

    auto UImanager = G4UImanager::GetUIpointer();
    if(ui) {
        // Interactive
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    } else {
        // Batch
        G4String macro = argv[1];
        UImanager->ApplyCommand("/control/execute "+macro);
    }

    delete vis;
    delete runManager;
}

