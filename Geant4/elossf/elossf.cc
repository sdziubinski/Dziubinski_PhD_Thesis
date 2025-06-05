//
//******************************************************************************
/// \author Sean Dziubinski
/// \date 03-04-2025
/// \file elossf.cc
/// \brief Main program of the elossf - a faster version of eloss_r8520_m
//******************************************************************************
//

// User header files
#include "EDetectorConstruction.hh"
#include "EActionInitialization.hh"

// Managers and executives
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

// Physics
#include "LBE.hh"
#include "G4OpticalPhysics.hh"

// Additional tools
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"
#include <iostream>
#include <sys/time.h>
#include <time.h>
#include <format>
#include <string>

using namespace std;
using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Declare input parameters
  int Z, A;
  string n;

  // Look for input parameters
  if (argc == 4)
  {
    Z = atoi(argv[1]);
    A = atoi(argv[2]);
    n = argv[3];
    G4cout<<G4endl
    <<"**************************************************************"<<G4endl
    <<"Input Parameters: "<<"Z="<<Z<<", A="<<A<<", n="<<n<<G4endl
    <<G4endl;
  }
  else
  {
    G4cout<<"Missing Input Parameters!!"<<G4endl
    <<"Need Z, A, and n!!!"<<G4endl;
    return 1;
  }

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new DetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physicsList = new LBE(0);
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->SetVerboseLevel(0);
  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);
  
  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(Z, A));

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Run simulation in batch mode 
  string sfmt = "/run/beamOn ";
  sfmt +=n;
  G4cout<<sfmt<<G4endl;
  G4cout<<"Running Batch Mode..."<<G4endl;
  UImanager->ApplyCommand("/run/initialize");
  UImanager->ApplyCommand(sfmt); // Choose number of events here!!

  delete runManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
