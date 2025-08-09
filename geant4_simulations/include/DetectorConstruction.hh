// DetectorConstruction.hh
#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4ThreeVector.hh>
#include <G4LogicalVolume.hh>
#include <G4GlobalMagFieldMessenger.hh>
#include <G4GenericMessenger.hh>
#include <G4TransportationManager.hh>
#include "G4RunManagerFactory.hh"

/*
// https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/material.html
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"


G4String name, symbol;             // a=mass of a mole;
G4double a, z, density;            // z=mean number of protons;
G4int iz, n;                       // iz=nb of protons  in an isotope;
                                   // n=nb of nucleons in an isotope;
G4int ncomponents, natoms;
G4double abundance, fractionmass;
G4double temperature, pressure;

G4UnitDefinition::BuildUnitsTable();

// define Elements
a = 1.01*g/mole;
G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
///////////////////////////////////////////////////
*/

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;
    G4VPhysicalVolume* Construct() override;

    void SetDepth(G4double d); // { fDepth = d; fConstructed = false; }
    G4double GetDepth() const { return fDepth; }

  private:
    G4double fDepth;          // rock thickness (m)
    G4bool   fConstructed;    // dirty-flag
    G4GenericMessenger* fMessenger;
};

#endif

