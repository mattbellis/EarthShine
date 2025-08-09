// DetectorConstruction.cc
#include "DetectorConstruction.hh"
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>


#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"


DetectorConstruction::DetectorConstruction()
 : fDepth(1000.*m), fConstructed(false)
{
    // UI command: /mydet/depth <value> <unit>
    fMessenger = new G4GenericMessenger(this,"/mydet/","Detector control");
    //fMessenger->DeclarePropertyWithUnit("depth","m",fDepth,
    //    "Rock thickness (1 m – 1e5 m)");

    // New
    fMessenger->DeclareMethodWithUnit( "depth","m", &DetectorConstruction::SetDepth,
    "Rock thickness (1 m – 1e5 m)");
}

void DetectorConstruction::SetDepth(G4double d)
{
    fDepth = d;
    fConstructed = false;                       // mark dirty
    //G4RunManager::GetRunManager()->GeometryHasChanged();
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
    if(fConstructed) return G4TransportationManager::GetTransportationManager()
                          ->GetNavigatorForTracking()->GetWorldVolume();

    auto nist = G4NistManager::Instance();

    auto granite = new G4Material("Granite", 2.65*g/cm3, 3);  // bulk density
    granite->AddMaterial(nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"), 0.75*perCent); // SiO2
    granite->AddMaterial(nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE"), 0.21*perCent); // Al2O3
    granite->AddMaterial(nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE"), 0.04*perCent); // K2O
    //auto granite = new G4Material("Granite", 2.65*g/cm3, 6);  // bulk density
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"), 0.72*perCent); // SiO2
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE"), 0.14*perCent); // Al2O3
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE"), 0.04*perCent); // K2O
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_SODIUM_OXIDE"),    0.04*perCent); // Na2O
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_CALCIUM_OXIDE"),   0.02*perCent); // CaO
    //granite->AddMaterial(nist->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE"), 0.01*perCent); // MgO

    auto rock = granite;

    //auto rock = nist->FindOrBuildMaterial("G4_GRANITE");
    // This works
    //auto rock = nist->FindOrBuildMaterial("G4_Si");

    // World size: add generous margin
    const G4double worldSize = fDepth + 20.*m;

    auto solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
    auto logicWorld = new G4LogicalVolume(solidWorld,
                                          nist->FindOrBuildMaterial("G4_AIR"),
                                          "World");

    auto physWorld  = new G4PVPlacement(nullptr, {}, logicWorld,
                                        "World", nullptr, false, 0);

    // Rock volume (a box in the center, thickness = fDepth along +Z)
    auto solidRock = new G4Box("Rock", 0.5*worldSize, 0.5*worldSize, 0.5*fDepth);
    auto logicRock = new G4LogicalVolume(solidRock, rock, "Rock");
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0,0,-0.25*worldSize), // put top at z=0
                      logicRock, "Rock", logicWorld, false, 0);

    fConstructed = true;
    return physWorld;
}

