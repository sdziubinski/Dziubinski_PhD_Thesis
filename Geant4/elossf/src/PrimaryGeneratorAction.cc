
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

// User header files
#include "EPrimaryGeneratorAction.hh"
#include "ERunAction.hh"

// Managers and executives
#include "G4RunManager.hh"

// Geant4 tools
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Threading.hh"
#include "G4IonTable.hh"

// C++ tools
#include <random>
#include <cmath>
#include <iostream>

using namespace std;

//G4Mutex myMutex=G4MUTEX_INITIALIZER;

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(int Z, int A)
: G4VUserPrimaryGeneratorAction()
{
  fParticleGun  = new G4ParticleGun(1);
  Z_pg = Z;
  A_pg = A;
  x0 = 0;
  y0 = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  //

  // Get worker thread ID
  int thr=G4Threading::G4GetThreadId();
  G4cout << "Thread ID: " << thr <<G4endl;

  // Define primary particle
  G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4double energy = A_pg*80*MeV;
  G4double excite_energy = 0.*keV;
  G4ParticleDefinition* particle = ionTable->GetIon(Z_pg,A_pg,excite_energy);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(energy);
  
  // Set starting position -> random integer within ELOSS effective area
  x0 = std::round(G4UniformRand() * max_x1 - 295);
  y0 = std::round(G4UniformRand() * max_y1 - 145);
  G4double z0 = -3000.0*mm;
  
  G4cout << "x0=" << x0 << ", y0=" << y0 <<G4endl;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0*mm,y0*mm,z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));

  // Generate particle gun for the event
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Define file name for starting position data
  string filename = "pos" + to_string(thr) + ".dat";

  // Write to the worker thread's own position data file
  //myMutex.lock();
  ofstream out;
  out.open(filename,ofstream::app);
  out<<x0<<" "<<y0<<G4endl;
  out.close();
  //myMutex.unlock();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


