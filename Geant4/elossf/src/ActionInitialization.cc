//
//******************************************************************************
/// \author Sean Dziubinski
/// \date 10-21-2024
/// \file ActionInitialization.cc
/// \brief ELOSSAI of eloss based off of eloss_r8520_m and example B1
//******************************************************************************
//

#include "EActionInitialization.hh"
#include "EPrimaryGeneratorAction.hh"
#include "ERunAction.hh"
#include "EEventAction.hh"
#include "ESteppingAction.hh"

//#include "G4Mutex.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(int Z, int A)
{
  Z_pg = Z;
  A_pg = A;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
  auto runAction = new RunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  auto generatorAction = new PrimaryGeneratorAction(Z_pg, A_pg);
  SetUserAction(generatorAction);
  auto runAction = new RunAction;
  SetUserAction(runAction);
  auto stepAction = new SteppingAction;
  SetUserAction(stepAction);
  SetUserAction(new EventAction(stepAction, runAction, generatorAction)); // EventAction needs access to SteppingAction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
