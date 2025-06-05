//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/include/DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace B1
{

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    //G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    //G4LogicalVolume* fScoringVolume = nullptr;
    void DefineMaterials();
    void ConstructLaboratory();
    void SensitiveDete1();
    //Logical Volumes
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* gas_log;
    G4LogicalVolume* gase_log;
    G4LogicalVolume* pressure_log;
    G4LogicalVolume* ele_log;
    G4LogicalVolume* ppt_log;
    G4LogicalVolume* pmt1_log;
    G4LogicalVolume* pmt2_log;
    G4LogicalVolume* pc1_log;
    G4LogicalVolume* pc2_log;
    //Physical Volumes
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* gas_phys;
    G4VPhysicalVolume* gase_phys;
    G4VPhysicalVolume* pressure_phys;
    G4VPhysicalVolume* ele_phys;
    G4VPhysicalVolume* ppt_phys;
    G4VPhysicalVolume* pmt_phys;
    G4VPhysicalVolume* pc_phys;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
