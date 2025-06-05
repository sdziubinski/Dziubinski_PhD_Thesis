
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

// User header files
#include "EEventAction.hh"
#include "ERunAction.hh"
#include "ESteppingAction.hh"
#include "EPrimaryGeneratorAction.hh"

// Managers and executives
#include "G4RunManager.hh"
#include "G4AccumulableManager.hh"
#include "G4SDManager.hh"

// Geant4 tools
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"

// C++ tools
#include <iostream>

// Action classes
class PrimaryGeneratorAction;
class SteppingAction;
class RunAction;

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(SteppingAction* stepAction, RunAction* runAction, PrimaryGeneratorAction* generatorAction)
: G4UserEventAction(),
  fStepAction(stepAction),
  fGeneratorAction(generatorAction),
  fRunAction(runAction)
{}

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void EventAction::BeginOfEventAction(const G4Event*)
{
  // Register multifunctional detectors
  G4SDManager* SDman = G4SDManager::GetSDMpointer(); 
  G4String detName = "IonPro1";
  G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detName));


  if ( mfd )
    {
      //--- Loop over the registered primitive scorers.
      for (G4int icol = 0; icol < mfd->GetNumberOfPrimitives(); icol++)
      {
        // Get Primitive Scorer object.
        G4VPrimitiveScorer* scorer = mfd->GetPrimitive(icol);
        // collection name and collectionID for HitsCollection,
        // where type of HitsCollection is G4THitsMap in case of primitive scorer.
        // The collection name is given by <MFD name>/<Primitive Scorer name>.
        G4String collectionName = scorer->GetName();
        //G4cout << "Scorer Name: " << scorer->GetName() << G4endl;
        G4String fullCollectionName = detName+ "/" +collectionName;
        G4int collectionID = SDman->GetCollectionID(fullCollectionName);
        //
        if ( collectionID >= 0 )
        {
          // Store obtained HitsCollection information into data members.
          // And, creates new G4THitsMap for accumulating quantities during RUN.
          theCollName.push_back(fullCollectionName);
          theCollID.push_back(collectionID);
          theRunMap.push_back(new G4THitsMap<G4double>(detName,collectionName));
          //G4cout << "** collection found wwwwwwwwwwwwwwwwwwwwwwwwwwwwww " << G4endl;
        }
        else
        {
          G4cout << "** collection " << fullCollectionName << " not found. "<< G4endl;
        } 
      }
    }


  counter += 1; 
  G4cout << "New Particle #" <<counter<< G4endl;
  fStepAction->Populate_Sig();    // Initialize a new sig vector at the start of each event
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* aEvent)
{
  // Dump photon counts from sig into photonN.dat 
  fStepAction->Dump_Stepping();  

  // Gets the energy deposited in the effective Xe volume
  RecordEvent(aEvent);

  //--- Clear HitsMap for RUN
  G4int Nmap = theRunMap.size();
  for ( G4int i = 0; i < Nmap; i++) // Loops over the number of different sensitive detectors
  {
    if(theRunMap[i] ) theRunMap[i]->clear();
  }
  theCollName.clear();
  theCollID.clear();
  theRunMap.clear();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


void EventAction::RecordEvent(const G4Event* aEvent)
{
  // Get worker thread ID
  int thr=G4Threading::G4GetThreadId();

  // Check for hits in this event
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  if (!HCE) return;

  //=======================================================
  // Sum up HitsMap of this Event  into HitsMap of this RUN
  //=======================================================
  G4int Ncol = theCollID.size();
  for ( G4int i = 0; i < Ncol ; i++ )
  {  // Loop over HitsCollection
    G4THitsMap<G4double>* EvtMap=0;
    if ( theCollID[i] >= 0 )
    {           // Collection is attached to HCE
      EvtMap = (G4THitsMap<G4double>*)(HCE->GetHC(theCollID[i])); 
    }
    else
    {
      G4cout <<" Error EvtMap Not Found "<< i << G4endl;
    }
    if ( EvtMap )  
    {
      std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *theRunMap[i] += *EvtMap;
      if ( EvtMap->entries() == 0)
      {
        G4cout << "No hits in this event :("<< G4endl;
      }
      else if ( EvtMap->entries() >= 1)
      {
        for(; itr != EvtMap->GetMap()->end(); itr++) 
        {
          // Open a file to dump energy deposit data into (specific to the worker thread)
          std::ofstream Pfile;
          string filename = "energy" + to_string(thr) + ".dat";
          Pfile.open(filename,std::ios::in | std::ios::app);
          G4double putin = *(itr->second)/MeV;

          // Report and record deposited energy
          G4cout << "Deposited Energy = " << putin << " MeV" << G4endl;
          Pfile << putin << G4endl;
          Pfile.close();
        }
      }
    }
  }
}

G4THitsMap<G4double>* EventAction::GetHitsMap(const G4String& fullName)
{
  G4int Ncol = theCollName.size();
  for ( G4int i = 0; i < Ncol; i++)
  {
    if ( theCollName[i] == fullName )
    {
      return theRunMap[i];
    }
  }
  return NULL;
}


}
