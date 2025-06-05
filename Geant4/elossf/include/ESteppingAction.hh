
/// \file B1/include/SteppingAction.hh
/// \brief Definition of the B1::SteppingAction class

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;

//using namespace std;
namespace B1
{

/// Stepping action class

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    ~SteppingAction();

    // method from the base class
    void UserSteppingAction(const G4Step*);
    void Populate_Sig();
    void Dump_Stepping();
    std::vector<int> Getsig();

  private:
    //void Root_load(G4int fNpho, G4String fileNames);
    G4String datai_file;
    std::vector<G4String> photo_files;
    std::vector<int>sig;
    std::vector<G4double> tempoo[120];
    double tg;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
