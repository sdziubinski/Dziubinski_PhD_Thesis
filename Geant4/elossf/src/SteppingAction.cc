
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "ESteppingAction.hh"
#include "EEventAction.hh"
#include "EDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "Randomize.hh"

#include <numeric>
#include <iostream>
#include <string>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  double pde = 0.28; // photon detection efficiency
  double r_rand;
  G4String post_volume;
  int pmt_num;
  
  // Update photon count for a PMT if an optical photon enters the PMT volume from the gas volume
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton" &&
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "gas_Vol")
  {
    post_volume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

    if (post_volume.substr(0,2) == "pc")
    {
      r_rand = G4UniformRand();
      if (r_rand<pde) 
      {
        
        pmt_num = std::stoi(post_volume.substr(2,post_volume.length()));
        sig.at(pmt_num-1)=sig.at(pmt_num-1)+1; 
      }
    }
  }
}
  
void SteppingAction::Populate_Sig()
{
  for (int i = 0; i < 120; i++) 
  {
    sig.push_back(0);
  }
}

void SteppingAction::Dump_Stepping()
{
  // Get worker thread ID
  int thr=G4Threading::G4GetThreadId();

  // Create worker thread specific file name
  std::string filename = "photon" + std::to_string(thr) + ".dat";
  
  // Sum photon counts for all 120 PMTs
  int sum = std::accumulate(sig.begin(), sig.end(), 0);
  G4cout<<"Sum: "<<sum<<G4endl;

  // Open and write to photonN.dat file
  std::ofstream Pfile;
  Pfile.open(filename,std::ios::in | std::ios::app);

  if (sum > 0)
  {
      for (int pos = 0; pos<120; pos++)
      {
          Pfile << sig.at(pos) << " ";
      }
      Pfile << G4endl;
      Pfile.close();
  }
  sig.clear();
}

}
