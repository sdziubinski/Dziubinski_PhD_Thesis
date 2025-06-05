
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "EDetectorConstruction.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDManager.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
 :  G4VUserDetectorConstruction(),experimentalHall_log(0),gas_log(0),gase_log(0),pressure_log(0),ele_log(0),ppt_log(0),pmt1_log(0),pmt2_log(0),
    pc1_log(0),pc2_log(0), 
    experimentalHall_phys(0),gas_phys(0),gase_phys(0),pressure_phys(0),ele_phys(0),ppt_phys(0),pmt_phys(0),pc_phys(0)
{;}

DetectorConstruction::~DetectorConstruction()
{;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  
  DefineMaterials();
  ConstructLaboratory();
  SensitiveDete1();
  //
  //always return the physical World
  //
  return experimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4Element *H  = new G4Element("Hydrogen","H",1.,1.0079*g/mole);
  G4Element *Xe = new G4Element("Xenon","Xe", 54., 131.293*g/mole);
  G4Element *N  = new G4Element("Nitrogen","N",7.,14.007*g/mole);
  G4Element *O  = new G4Element("Oxygen","O",8.,15.999*g/mole);
  G4Element* Al = new G4Element("Aluminum","Al",13.,26.98*g/mole);
  G4Element *C  = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* I  = new G4Element("Iodine","I",53.,126.90447*g/mole);
  G4Element* Cs = new G4Element("Cesium","Cs",55.,132.90543*g/mole);

  //-------------------------------- Vacuum ------------------------------------
  G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3,2,kStateGas);
  Vacuum->AddElement(N, 0.755);
  Vacuum->AddElement(O, 0.245);

  //-------------------------------- Kevlar ------------------------------------
  G4Material *Kevlar = new G4Material("Kevlar",1.44*g/cm3,4,kStateSolid);
  Kevlar->AddElement(C, 14);
  Kevlar->AddElement(H, 10);
  Kevlar->AddElement(O, 2);
  Kevlar->AddElement(N, 2);
  
  //----------------------------- Polypropylene --------------------------------
  G4Material *Polypropylene = new G4Material("Polypropylene",0.94*g/cm3, 2,kStateSolid);
  Polypropylene->AddElement(C, 2);
  Polypropylene->AddElement(H, 4);
 
  //---------------------------- CsI Photocathode ------------------------------
  G4Material* CsI = new G4Material("CsI",4.51*g/cm3,2,kStateSolid);
  CsI-> AddElement(Cs, 1);
  CsI-> AddElement(I, 1);  
 
  //----------------------------- Gaseous Xenon --------------------------------
  //G4Material *GXe = new G4Material("GXe", 0.0062*g/cm3, 1, kStateGas); //, 298*kelvin , 1.0*atmosphere); /800 Torr 0.0062*g/cm3
  G4Material *GXe = new G4Material("GXe", 0.005887*g/cm3, 1, kStateGas);
  GXe->AddElement(Xe, 1);
  
  const G4int iNbEntries = 3;
  G4double pdGXePhotonMomentum[iNbEntries]   = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdGXeScintillation[iNbEntries]    = {0.1,     1.0,     0.1};
  G4double pdGXeRefractiveIndex[iNbEntries]  = {1.00,    1.00,    1.00};
  G4double pdGXeAbsorbtionLength[iNbEntries] = {100*m,   100*m,   100*m};
  G4double pdGXeScatteringLength[iNbEntries] = {100*m,   100*m,   100*m};
 
  G4MaterialPropertiesTable *pGXePropertiesTable = new G4MaterialPropertiesTable(); 
 
  pGXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
  pGXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
  pGXePropertiesTable->AddProperty("RINDEX", pdGXePhotonMomentum, pdGXeRefractiveIndex, iNbEntries);
  pGXePropertiesTable->AddProperty("ABSLENGTH", pdGXePhotonMomentum, pdGXeAbsorbtionLength, iNbEntries);
  pGXePropertiesTable->AddProperty("RAYLEIGH", pdGXePhotonMomentum, pdGXeScatteringLength, iNbEntries);

  pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 13.0/(keV)); //VALUS ASSUMED 13 ph/KeV for electroncs and 28 ph/KeV for alpha, 14 ph/KeV conservative
  //pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 26.0/(keV)); //Henriques 2024
  pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
  pGXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 4.1*ns);
  pGXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 21.*ns);
  pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

  GXe->SetMaterialPropertiesTable(pGXePropertiesTable); 

  // --------------------------- Aluminum - Al ---------------------------------
  G4Material* Aluminum = new G4Material("Aluminum",2.7*g/cm3,1,kStateSolid);
  Aluminum->AddElement(Al, 1.0);

  // ------------------------------------ Aluminum - Al ------------------------------------
  //G4Material* Oxygen_gas = new G4Material("Oxygen_gas",0.004413*g/cm3,1,kStateGas); //, 298*kelvin, 1.0*atmosphere);
  //Oxygen_gas->AddElement(O, 1.0);
  //G4double pdO2PhotonMomentum[iNbEntries]   = {6.91*eV, 6.98*eV, 7.05*eV};
  //G4double pdO2eAbsorbtionLength[iNbEntries] = {2.2*cm,   2.2*cm, 2.2*m};

  //G4MaterialPropertiesTable *pO2PropertiesTable = new G4MaterialPropertiesTable(); 
  //pO2PropertiesTable->AddProperty("ABSLENGTH", pdO2PhotonMomentum, pdO2eAbsorbtionLength, iNbEntries);
  //Oxygen_gas->SetMaterialPropertiesTable(pO2PropertiesTable); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructLaboratory()
{
  G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
  G4Material *GXe = G4Material::GetMaterial("GXe");
  //G4Material *Oxygen_gas = G4Material::GetMaterial("Oxygen_gas");
  G4Material *Aluminum = G4Material::GetMaterial("Aluminum");
  G4Material *Kevlar = G4Material::GetMaterial("Kevlar");
  G4Material *Polypropylene = G4Material::GetMaterial("Polypropylene"); 
  G4Material *CsI = G4Material::GetMaterial("CsI"); 
  //G4Material *Filling_gas = new G4Material("Filling_gas",0.005887*g/cm3,2);
  //Filling_gas->AddMaterial(GXe,99.999*perCent);
  //Filling_gas->AddMaterial(Oxygen_gas,0.001*perCent);
  
  G4double gap_pmt = -71.25*mm;//28.0*mm;
  G4double gap_pc = -71.25*mm;//28.0*mm;

  // Experimental hall (world volume) --------------
  G4double expHall_x = 800.0*mm;
  G4double expHall_y = 500.0*mm;
  G4double expHall_z = 20.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),experimentalHall_log,"expHall",0,false,0);

  //GXe volume
  G4double v_x = 400.0*mm;
  G4double v_y = 230.0*mm;
  G4double v_z = 120.0*mm;
  G4Box* gas_box = new G4Box("GXe Volume",v_x,v_y,v_z);
  gas_log = new G4LogicalVolume(gas_box,GXe,"gas_log",0,0,0);
  gas_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),gas_log,"gas_Vol",experimentalHall_log,false,0);

  //GXe volume effective
  G4double ve_x = 300.0*mm;
  G4double ve_y = 150.0*mm;
  G4double ve_z = ((17.5/2)+(172.5/2))*mm;
  G4Box* gase_box = new G4Box("GXe effective Volume",ve_x,ve_y,ve_z);
  gase_log = new G4LogicalVolume(gase_box,GXe,"gase_log",0,0,0);
  gase_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0),gase_log,"gase_Vol",gas_log,false,0);

  //pressure windows
  G4double p_z = 55.0*um; // assuming 15 mg/cm^2 thickness
  G4Box* pressure_box = new G4Box("Pressure Window",v_x,v_y,p_z);
  pressure_log = new G4LogicalVolume(pressure_box,Kevlar,"pressure_log",0,0,0);
  pressure_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,p_z+v_z),pressure_log,"pressure_Vol1",experimentalHall_log,false,0);
  pressure_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-p_z-v_z),pressure_log,"pressure_Vol2",experimentalHall_log,false,1);

  //electrode foils Al
  G4double ele_z = 0.075*um; // assuming 150 nm thickess evaporated Al
  G4Box* ele_box = new G4Box("electrode Volume",ve_x,ve_y,ele_z);
  ele_log = new G4LogicalVolume(ele_box,Aluminum,"ele_log",0,0,0);
  //ele_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ve_z+ele_z),ele_log,"cathode_Vol",gas_log,false,0);
  //ele_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-ve_z-ele_z),ele_log,"anode_Vol",gas_log,false,1);

  //electrode foils Polypropilene
  G4double ppt_z = 0.4*um; // assuming 70 um/cm^2 thickness
  G4Box* ppt_box = new G4Box("polypropylene Volume",ve_x,ve_y,ppt_z);
  ppt_log = new G4LogicalVolume(ppt_box,Polypropylene,"ppt_log",0,0,0);
  ppt_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0),ppt_log,"ppt_cathode_Vol",gase_log,false,0);
  ppt_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,47.5*mm-0.4*um),ppt_log,"ppt_anode_Vol",gase_log,false,1);
  ppt_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-47.5*mm+0.4*um),ppt_log,"ppt_anode_Vol",gase_log,false,2);
  ppt_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,2*47.5*mm-2*0.4*um),ppt_log,"ppt_anode_Vol",gase_log,false,3);
  ppt_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-2*47.5*mm+2*0.4*um),ppt_log,"ppt_anode_Vol",gase_log,false,4);
  
  //pmt R8520
  G4double pmt_lunga = 15*mm;//28.0*mm;
  G4double pmt_corta = 7.5*mm;
  G4double pmt_z = 15*mm;//28.0*mm;
  G4double gap = 60*mm;//28.0*mm;
  G4double gap_z = 47.5*2*mm;//28.0*mm;
  G4Box* pmt_box1 = new G4Box("pmt1",pmt_lunga,pmt_corta,pmt_z);
  G4Box* pmt_box2 = new G4Box("pmt2",pmt_corta,pmt_lunga,pmt_z);
  pmt1_log = new G4LogicalVolume(pmt_box1,GXe,"pmt1_log",0,0,0);
  pmt2_log = new G4LogicalVolume(pmt_box2,GXe,"pmt2_log",0,0,0);
  //fist sector
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt1",gas_log,false,0);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt2",gas_log,false,1);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt3",gas_log,false,2);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt4",gas_log,false,3);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt5",gas_log,false,4);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt6",gas_log,false,5);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt7",gas_log,false,6);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt8",gas_log,false,7);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt9",gas_log,false,8);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,-165.0*mm,0.0+gap_pmt),pmt1_log,"pmt10",gas_log,false,9);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+0*gap,0.0+gap_pmt),pmt2_log,"pmt11",gas_log,false,10);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+1*gap,0.0+gap_pmt),pmt2_log,"pmt12",gas_log,false,11);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+2*gap,0.0+gap_pmt),pmt2_log,"pmt13",gas_log,false,12);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+3*gap,0.0+gap_pmt),pmt2_log,"pmt14",gas_log,false,13);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+4*gap,0.0+gap_pmt),pmt2_log,"pmt15",gas_log,false,14);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt16",gas_log,false,15);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt17",gas_log,false,16);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt18",gas_log,false,17);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt19",gas_log,false,18);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt20",gas_log,false,19);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt21",gas_log,false,20);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt22",gas_log,false,21);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt23",gas_log,false,22);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt24",gas_log,false,23);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,165.0*mm,0.0+gap_pmt),pmt1_log,"pmt25",gas_log,false,24);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+4*gap,0.0+gap_pmt),pmt2_log,"pmt26",gas_log,false,25);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+3*gap,0.0+gap_pmt),pmt2_log,"pmt27",gas_log,false,26);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+2*gap,0.0+gap_pmt),pmt2_log,"pmt28",gas_log,false,27);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+1*gap,0.0+gap_pmt),pmt2_log,"pmt29",gas_log,false,28);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+0*gap,0.0+gap_pmt),pmt2_log,"pmt30",gas_log,false,29);  
  //second sector
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt31",gas_log,false,30);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt32",gas_log,false,31);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt33",gas_log,false,32);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt34",gas_log,false,33);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt35",gas_log,false,34);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt36",gas_log,false,35);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt37",gas_log,false,36);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt38",gas_log,false,37);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt39",gas_log,false,38);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,-165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt40",gas_log,false,39);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+0*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt41",gas_log,false,40);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+1*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt42",gas_log,false,41);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+2*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt43",gas_log,false,42);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+3*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt44",gas_log,false,43);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+4*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt45",gas_log,false,44);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt46",gas_log,false,45);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt47",gas_log,false,46);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt48",gas_log,false,47);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt49",gas_log,false,48);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt50",gas_log,false,49);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt51",gas_log,false,50);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt52",gas_log,false,51);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt53",gas_log,false,52);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt54",gas_log,false,53);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,165.0*mm,gap_z/2+gap_pmt),pmt1_log,"pmt55",gas_log,false,54);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+4*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt56",gas_log,false,55);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+3*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt57",gas_log,false,56);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+2*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt58",gas_log,false,57);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+1*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt59",gas_log,false,58);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+0*gap+gap/2,gap_z/2+gap_pmt),pmt2_log,"pmt60",gas_log,false,59);  
  //third sector
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt61",gas_log,false,60);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt62",gas_log,false,61);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt63",gas_log,false,62);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt64",gas_log,false,63);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt65",gas_log,false,64);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt66",gas_log,false,65);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt67",gas_log,false,66);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt68",gas_log,false,67);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt69",gas_log,false,68);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,-165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt70",gas_log,false,69);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+0*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt71",gas_log,false,70);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+1*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt72",gas_log,false,71);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+2*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt73",gas_log,false,72);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+3*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt74",gas_log,false,73);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+4*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt75",gas_log,false,74);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt76",gas_log,false,75);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt77",gas_log,false,76);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt78",gas_log,false,77);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt79",gas_log,false,78);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt80",gas_log,false,79);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt81",gas_log,false,80);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt82",gas_log,false,81);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt83",gas_log,false,82);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt84",gas_log,false,83);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,165.0*mm,2*gap_z/2+gap_pmt),pmt1_log,"pmt85",gas_log,false,84);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+4*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt86",gas_log,false,85);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+3*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt87",gas_log,false,86);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+2*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt88",gas_log,false,87);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+1*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt89",gas_log,false,88);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+0*gap,2*gap_z/2+gap_pmt),pmt2_log,"pmt90",gas_log,false,89);  
  //forth sector
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt91",gas_log,false,90);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt92",gas_log,false,91);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt93",gas_log,false,92);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt94",gas_log,false,93);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt95",gas_log,false,94);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt96",gas_log,false,95);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt97",gas_log,false,96);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt98",gas_log,false,97);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt99",gas_log,false,98);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,-165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt100",gas_log,false,99);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+0*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt101",gas_log,false,100);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+1*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt102",gas_log,false,101);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+2*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt103",gas_log,false,102);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+3*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt104",gas_log,false,103);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(315.0*mm,-ve_y+pmt_lunga+4*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt105",gas_log,false,104);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt106",gas_log,false,105);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt107",gas_log,false,106);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt108",gas_log,false,107);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt109",gas_log,false,108);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt110",gas_log,false,109);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt111",gas_log,false,110);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt112",gas_log,false,111);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt113",gas_log,false,112);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt114",gas_log,false,113);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,165.0*mm,3*gap_z/2+gap_pmt),pmt1_log,"pmt115",gas_log,false,114);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+4*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt116",gas_log,false,115);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+3*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt117",gas_log,false,116);
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+2*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt118",gas_log,false,117);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+1*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt119",gas_log,false,118);  
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(-315.0*mm,-ve_y+pmt_lunga+0*gap+gap/2,3*gap_z/2+gap_pmt),pmt2_log,"pmt120",gas_log,false,119);  

  //Photocathode
  G4double pc_lunga = 10.25*mm;//28.0*mm;
  G4double pc_corta = 75*nm;
  G4double pc_z = 10.25*mm;//28.0*mm;
  G4double posizione_165m = -165.0*mm+pmt_corta+pc_corta;
  G4double posizione_165p = 165.0*mm-pmt_corta-pc_corta;
  G4double posizione_315m = -315.0*mm+pmt_corta+pc_corta;
  G4double posizione_315p = 315.0*mm-pmt_corta-pc_corta;
  G4Box* pc_box1 = new G4Box("pc1",pc_lunga,pc_corta,pc_z);
  G4Box* pc_box2 = new G4Box("pc2",pc_corta,pc_lunga,pc_z);
  pc1_log = new G4LogicalVolume(pc_box1,CsI,"pc1_log",0,0,0);
  pc2_log = new G4LogicalVolume(pc_box2,CsI,"pc2_log",0,0,0);
  //fist sector
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc1",gas_log,false,0);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc2",gas_log,false,1);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc3",gas_log,false,2);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc4",gas_log,false,3);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc5",gas_log,false,4);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc6",gas_log,false,5);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc7",gas_log,false,6);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc8",gas_log,false,7);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc9",gas_log,false,8);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,posizione_165m,0.0+gap_pc),pc1_log,"pc10",gas_log,false,9);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+0*gap,0.0+gap_pc),pc2_log,"pc11",gas_log,false,10);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+1*gap,0.0+gap_pc),pc2_log,"pc12",gas_log,false,11);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+2*gap,0.0+gap_pc),pc2_log,"pc13",gas_log,false,12);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+3*gap,0.0+gap_pc),pc2_log,"pc14",gas_log,false,13);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+4*gap,0.0+gap_pc),pc2_log,"pc15",gas_log,false,14);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc16",gas_log,false,15);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc17",gas_log,false,16);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc18",gas_log,false,17);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc19",gas_log,false,18);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc20",gas_log,false,19);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc21",gas_log,false,20);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc22",gas_log,false,21);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc23",gas_log,false,22);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc24",gas_log,false,23);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,posizione_165p,0.0+gap_pc),pc1_log,"pc25",gas_log,false,24);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+4*gap,0.0+gap_pc),pc2_log,"pc26",gas_log,false,25);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+3*gap,0.0+gap_pc),pc2_log,"pc27",gas_log,false,26);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+2*gap,0.0+gap_pc),pc2_log,"pc28",gas_log,false,27);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+1*gap,0.0+gap_pc),pc2_log,"pc29",gas_log,false,28);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+0*gap,0.0+gap_pc),pc2_log,"pc30",gas_log,false,29);  
  //second sector
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc31",gas_log,false,30);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc32",gas_log,false,31);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc33",gas_log,false,32);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc34",gas_log,false,33);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc35",gas_log,false,34);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc36",gas_log,false,35);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc37",gas_log,false,36);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc38",gas_log,false,37);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc39",gas_log,false,38);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,posizione_165m,gap_z/2+gap_pc),pc1_log,"pc40",gas_log,false,39);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+0*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc41",gas_log,false,40);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+1*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc42",gas_log,false,41);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+2*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc43",gas_log,false,42);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+3*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc44",gas_log,false,43);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+4*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc45",gas_log,false,44);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc46",gas_log,false,45);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc47",gas_log,false,46);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc48",gas_log,false,47);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc49",gas_log,false,48);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc50",gas_log,false,49);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc51",gas_log,false,50);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc52",gas_log,false,51);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc53",gas_log,false,52);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc54",gas_log,false,53);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,posizione_165p,gap_z/2+gap_pc),pc1_log,"pc55",gas_log,false,54);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+4*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc56",gas_log,false,55);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+3*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc57",gas_log,false,56);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+2*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc58",gas_log,false,57);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+1*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc59",gas_log,false,58);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+0*gap+gap/2,gap_z/2+gap_pc),pc2_log,"pc60",gas_log,false,59);  
  //third sector
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc61",gas_log,false,60);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc62",gas_log,false,61);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc63",gas_log,false,62);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc64",gas_log,false,63);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc65",gas_log,false,64);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc66",gas_log,false,65);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc67",gas_log,false,66);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc68",gas_log,false,67);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc69",gas_log,false,68);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,posizione_165m,2*gap_z/2+gap_pc),pc1_log,"pc70",gas_log,false,69);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+0*gap,2*gap_z/2+gap_pc),pc2_log,"pc71",gas_log,false,70);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+1*gap,2*gap_z/2+gap_pc),pc2_log,"pc72",gas_log,false,71);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+2*gap,2*gap_z/2+gap_pc),pc2_log,"pc73",gas_log,false,72);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+3*gap,2*gap_z/2+gap_pc),pc2_log,"pc74",gas_log,false,73);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+4*gap,2*gap_z/2+gap_pc),pc2_log,"pc75",gas_log,false,74);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc76",gas_log,false,75);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc77",gas_log,false,76);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc78",gas_log,false,77);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc79",gas_log,false,78);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc80",gas_log,false,79);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc81",gas_log,false,80);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc82",gas_log,false,81);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc83",gas_log,false,82);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc84",gas_log,false,83);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap,posizione_165p,2*gap_z/2+gap_pc),pc1_log,"pc85",gas_log,false,84);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+4*gap,2*gap_z/2+gap_pc),pc2_log,"pc86",gas_log,false,85);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+3*gap,2*gap_z/2+gap_pc),pc2_log,"pc87",gas_log,false,86);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+2*gap,2*gap_z/2+gap_pc),pc2_log,"pc88",gas_log,false,87);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+1*gap,2*gap_z/2+gap_pc),pc2_log,"pc89",gas_log,false,88);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+0*gap,2*gap_z/2+gap_pc),pc2_log,"pc90",gas_log,false,89);  
  //forth sector
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc91",gas_log,false,90);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc92",gas_log,false,91);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc93",gas_log,false,92);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc94",gas_log,false,93);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc95",gas_log,false,94);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc96",gas_log,false,95);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc97",gas_log,false,96);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc98",gas_log,false,97);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc99",gas_log,false,98);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,posizione_165m,3*gap_z/2+gap_pc),pc1_log,"pc100",gas_log,false,99);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+0*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc101",gas_log,false,100);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+1*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc102",gas_log,false,101);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+2*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc103",gas_log,false,102);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+3*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc104",gas_log,false,103);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315p,-ve_y+pmt_lunga+4*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc105",gas_log,false,104);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+9*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc106",gas_log,false,105);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+8*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc107",gas_log,false,106);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+7*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc108",gas_log,false,107);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+6*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc109",gas_log,false,108);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+5*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc110",gas_log,false,109);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+4*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc111",gas_log,false,110);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+3*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc112",gas_log,false,111);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+2*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc113",gas_log,false,112);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+1*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc114",gas_log,false,113);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(-ve_x+pmt_lunga+0*gap+gap/2,posizione_165p,3*gap_z/2+gap_pc),pc1_log,"pc115",gas_log,false,114);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+4*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc116",gas_log,false,115);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+3*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc117",gas_log,false,116);
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+2*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc118",gas_log,false,117);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+1*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc119",gas_log,false,118);  
  pc_phys = new G4PVPlacement(0,G4ThreeVector(posizione_315m,-ve_y+pmt_lunga+0*gap+gap/2,3*gap_z/2+gap_pc),pc2_log,"pc120",gas_log,false,119);  
  G4VisAttributes * pc_VisAtt = new G4VisAttributes(G4Colour(0,1.0,0));
  pc_VisAtt->SetForceSolid(true);
  pc_VisAtt->SetForceAuxEdgeVisible(true);
  pc1_log->SetVisAttributes(pc_VisAtt);
  pc2_log->SetVisAttributes(pc_VisAtt);

  // ------------------------------- Surfaces -----------------------------------------------
  // ------------------------- Optical Propertiers ------------------------------------------
  const G4int iNbEntries_Al = 3;
  G4OpticalSurface* mirror = new G4OpticalSurface("mirror",glisur,polished,dielectric_metal);
  G4double AlPM[iNbEntries_Al]  = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double AlReflectivity[iNbEntries_Al] = {0.8, 0.8, 0.8};
  G4double AlEfficiency[iNbEntries_Al] = {0.0, 0.0, 0.0};
  G4MaterialPropertiesTable *AlPropertiesTable = new G4MaterialPropertiesTable();
  AlPropertiesTable->AddProperty("REFLECTIVITY", AlPM, AlReflectivity, iNbEntries_Al);
  AlPropertiesTable->AddProperty("EFFICIENCY", AlPM, AlEfficiency, iNbEntries_Al);
  mirror->SetMaterialPropertiesTable(AlPropertiesTable); 
  new G4LogicalSkinSurface("mirror_surf",ppt_log,mirror);

  //fScoringVolume = gase_log;
  
}

void DetectorConstruction::ConstructSDandField()
{
  SensitiveDete1();
}

void DetectorConstruction::SensitiveDete1()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->SetVerboseLevel(1);

  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("IonPro1");
  SDman->AddNewDetector(det);
  gase_log->SetSensitiveDetector(det);

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("nproton");

  G4String filterName,particleName;
  det->RegisterPrimitive(primitive);
}

}



