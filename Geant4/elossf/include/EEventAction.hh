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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
//#include "G4Accumulable.hh"
#include "globals.hh"
#include "G4THitsMap.hh"

//#include <TFile.h>
//#include <TTree.h>

namespace B1
{

class RunAction;
class SteppingAction;
class PrimaryGeneratorAction;

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
    EventAction(SteppingAction* stepAction, RunAction* runAction, PrimaryGeneratorAction* generatorAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    virtual void RecordEvent(const G4Event* event);
    G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);

    int counter;

  private:
    SteppingAction* fStepAction = nullptr;
    RunAction* fRunAction = nullptr;
    PrimaryGeneratorAction* fGeneratorAction = nullptr;
    std::vector<G4String> theCollName;
    std::vector<G4int> theCollID;
    std::vector<G4THitsMap<G4double>*> theRunMap;
    G4THitsMap<G4double> mapSum;
    
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


