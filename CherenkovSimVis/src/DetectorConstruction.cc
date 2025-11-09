// DetectorConstruction.cc - FIXED VERSION
//
// Key fixes:
// 1. Added photocathode volume inside PMT glass
// 2. Defined optical border surface with EFFICIENCY property using QE curve
// 3. Attached sensitive detector to photocathode, not glass
// 4. Added proper optical photon detection at photocathode surface

#include "DetectorConstruction.hh"
#include "PMTObject.hh"
#include "PMTSensitiveDetector.hh"
#include "SimplePMTResponse.hh"
#include "G4SDManager.hh"

// Geant4 includes
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fWorldLogical(nullptr),
fTankLogical(nullptr),
fInnerDetectorLogical(nullptr),
fPMTLogical(nullptr),
fPhotocathodeLogical(nullptr),
fWorldPhysical(nullptr),
fTankPhysical(nullptr),
fInnerDetectorPhysical(nullptr),
fGdDopedWater(nullptr),
fSteel(nullptr),
fPMTGlass(nullptr),
fPhotocathodeMaterial(nullptr),
fWaterSteelSurface(nullptr),
fPhotocathodeOpticalSurface(nullptr),
fVisWorld(nullptr),
fVisTank(nullptr),
fVisInnerDetector(nullptr),
fVisPMT(nullptr),
fVisPhotocathode(nullptr),
fTankOuterRadius(16.9 * m),
fTankHeight(36.2 * m),
fInnerRadius(16.0 * m),
fInnerHeight(33.8 * m),
fNumPMTs(11146),
fPMTRadius(0.254 * m),
fPMTGlassThickness(0.4 * cm),
fPMTPlacementOffset(1.0 * cm),
fGdMassPercent(0.2 * perCent)
{
fPMTModel = std::make_unique<CherenkovSim::PMT::PMT20inch>();
fPMTRadius = fPMTModel->GetRadius_m() * m;
fPMTGlassThickness = fPMTModel->GetGlassThickness_m() * m;
}

DetectorConstruction::~DetectorConstruction() = default;

G4VPhysicalVolume* DetectorConstruction::Construct()
{
DefineMaterials();
BuildGeometry();
DefineOpticalSurfaces();
return fWorldPhysical;
}
//------------define materials///////////
void DetectorConstruction::DefineMaterials()
{
auto* nist = G4NistManager::Instance();

fSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
fPMTGlass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");
if (!fPMTGlass)
    fPMTGlass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");

// Add optical properties to PMT Glass
const G4int n_glass = 2;
G4double energy_glass[n_glass] = {2.0*eV, 3.5*eV}; // A simple range for visible light
G4double rindex_glass[n_glass] = {1.47, 1.47};     // Refractive index for Pyrex
G4double abs_glass[n_glass] = {50.0*m, 50.0*m};   // Make it very transparent

auto* mpt_glass = new G4MaterialPropertiesTable();
mpt_glass->AddProperty("RINDEX", energy_glass, rindex_glass, n_glass);
mpt_glass->AddProperty("ABSLENGTH", energy_glass, abs_glass, n_glass);
fPMTGlass->SetMaterialPropertiesTable(mpt_glass);

auto* water = nist->FindOrBuildMaterial("G4_WATER");

// Gadolinium element
auto* elGd = new G4Element("Gadolinium", "Gd", 64., 157.25 * g / mole);

// Gd-doped water
const G4double fracGd = fGdMassPercent;
const G4double fracWater = 1.0 - fracGd;
const G4double density = 1.00 * g / cm3;
fGdDopedWater = new G4Material("GdDopedWater", density, 2);
fGdDopedWater->AddMaterial(water, fracWater);
fGdDopedWater->AddElement(elGd, fracGd);

 // Water optical properties
const G4int nEntries = 5;
G4double energy[nEntries] = {2.00 * eV, 2.40 * eV, 2.75 * eV, 3.10 * eV, 3.50 * eV};
G4double rindex[nEntries] = {1.349, 1.345, 1.342, 1.339, 1.336};

// Set a single, constant absorption length for our measurement.
// 80m is a good, realistic value for clean water.
const G4double fixedAbsLength = 210.0 * m;
G4double absLength[nEntries];
for (int i = 0; i < nEntries; ++i) {
    absLength[i] = fixedAbsLength;
}

auto* mpt = new G4MaterialPropertiesTable();
// 1. Adding Refractive Index 
mpt->AddProperty("RINDEX", energy, rindex, nEntries);
// 2. Adding Absorption Length
mpt->AddProperty("ABSLENGTH", energy, absLength, nEntries);

// 3. ADDING THE NEW RAYLEIGH SCATTERING LENGTHS HERE
const G4int nRayleigh = 5;
G4double energy_rayleigh[nRayleigh]   = {2.00*eV, 2.40*eV, 2.75*eV, 3.10*eV, 3.50*eV};
G4double scattering_length[nRayleigh] = {210.0*m,  210.0*m,  210.0*m,  210.0*m,  210.0*m};
mpt->AddProperty("RAYLEIGH", energy_rayleigh, scattering_length, nRayleigh);

//// 4. Now, assign the completed table (with all 3 properties) to the water material
fGdDopedWater->SetMaterialPropertiesTable(mpt);

// **FIX 1: Create photocathode material (vacuum or low-density gas)**
fPhotocathodeMaterial = new G4Material("Photocathode", 1.0, 1.01*g/mole, 1e-25*g/cm3, kStateGas, 300*kelvin, 1e-19*pascal);

// Photocathode needs RINDEX for optical surface to work
auto* photocathodeMPT = new G4MaterialPropertiesTable();
G4double photocathodeRindex[nEntries] = {1.47, 1.47, 1.47, 1.47, 1.47}; // change from 1.0 ti 1.47
photocathodeMPT->AddProperty("RINDEX", energy, photocathodeRindex, nEntries);
fPhotocathodeMaterial->SetMaterialPropertiesTable(photocathodeMPT);

G4cout << "[DetectorConstruction] Materials defined with photocathode material" << G4endl;

}
//-------------------build geom------------------
void DetectorConstruction::BuildGeometry()
{
auto* nist = G4NistManager::Instance();

// World
const G4double worldSize = std::max({fTankOuterRadius * 1.2, fTankHeight * 1.2, 50.0 * m});
auto* solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
fWorldLogical = new G4LogicalVolume(solidWorld, nist->FindOrBuildMaterial("G4_AIR"), "WorldLogical"); //hmmm
// --- GET THE MATERIAL AND ADD PROPERTIES TO IT ---

G4Material* world_air = nist->FindOrBuildMaterial("G4_AIR");

// Create a new properties table for Air
auto* mpt_air = new G4MaterialPropertiesTable();

const G4int n_air = 2;
G4double energy_air[n_air] = {2.0*eV, 3.5*eV};
G4double rindex_air[n_air] = {1.0, 1.0}; // Refractive index of air is ~1.0
mpt_air->AddProperty("RINDEX", energy_air, rindex_air, n_air);

world_air->SetMaterialPropertiesTable(mpt_air);

// Now, make sure the logical volume uses this material pointer
fWorldLogical = new G4LogicalVolume(solidWorld, world_air, "WorldLogical"); // Re-assign with the modified material
fWorldPhysical = new G4PVPlacement(nullptr, {}, fWorldLogical, "WorldPhysical", nullptr, false, 0, true);

fVisWorld = new G4VisAttributes(G4Colour(1, 1, 1, 0));
fVisWorld->SetVisibility(false);
fWorldLogical->SetVisAttributes(fVisWorld);

// Steel Tank
auto* solidTank = new G4Tubs("Tank", fInnerRadius, fTankOuterRadius, fTankHeight / 2.0, 0, 360 * deg);
fTankLogical = new G4LogicalVolume(solidTank, fSteel, "TankLogical");
fTankPhysical = new G4PVPlacement(nullptr, {}, fTankLogical, "TankPhysical", fWorldLogical, false, 0, true);

fVisTank = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.2));
fVisTank->SetForceSolid(true);
fTankLogical->SetVisAttributes(fVisTank);

// Inner Detector (Gd-doped water)
auto* solidInner = new G4Tubs("InnerDetector", 0, fInnerRadius, fInnerHeight / 2.0, 0, 360 * deg);
fInnerDetectorLogical = new G4LogicalVolume(solidInner, fGdDopedWater, "InnerDetectorLogical");
fInnerDetectorPhysical = new G4PVPlacement(nullptr, {}, fInnerDetectorLogical,
                                           "InnerDetectorPhysical", fWorldLogical, false, 0, true);

fVisInnerDetector = new G4VisAttributes(G4Colour(0.0, 0.4, 1.0, 0.08));
fVisInnerDetector->SetForceSolid(true);
fInnerDetectorLogical->SetVisAttributes(fVisInnerDetector);

// **FIX 2: PMT glass shell (outer)**
const G4double pmtR = fPMTRadius;
const G4double glassT = fPMTGlassThickness;

// Glass shell: from pmtR-glassT to pmtR (hollow sphere cap)
auto* glassShell = new G4Sphere("PMTGlass", pmtR - glassT, pmtR, 0, 360 * deg, 0, 90 * deg);
fPMTLogical = new G4LogicalVolume(glassShell, fPMTGlass, "PMTLogical");

fVisPMT = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.3));
fVisPMT->SetForceSolid(true);
fPMTLogical->SetVisAttributes(fVisPMT);

// **FIX 3: Photocathode volume (thin hemisphere inside glass)**
const G4double photocathodeThickness = 0.1 * mm; // very thin
const G4double photocathodeOuterR = pmtR - glassT;
const G4double photocathodeInnerR = pmtR - glassT - photocathodeThickness;


auto* photocathodeSolid = new G4Sphere("Photocathode", 
                                       photocathodeInnerR, 
                                       photocathodeOuterR, 
                                       0, 360 * deg, 0, 90 * deg);
fPhotocathodeLogical = new G4LogicalVolume(photocathodeSolid, fPhotocathodeMaterial, "PhotocathodeLogical");

fVisPhotocathode = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5)); // red for photocathode
fVisPhotocathode->SetForceSolid(true);
fPhotocathodeLogical->SetVisAttributes(fVisPhotocathode);

// Place photocathode inside PMT (will be placed for each PMT)
// We'll do this in PlacePMTs()

// PMT black frame
const G4double frameThickness = 5.0 * mm;
const G4double frameRadius = fPMTRadius * 1.05;
auto* solidFrame = new G4Tubs("PMTFrame", 0, frameRadius, frameThickness / 2.0, 0, 360 * deg);
auto* blackMat = nist->FindOrBuildMaterial("G4_GRAPHITE");
fPMTFrameLogical = new G4LogicalVolume(solidFrame, blackMat, "PMTFrameLogical");
auto* visFrame = new G4VisAttributes(G4Colour(0, 0, 0));
visFrame->SetForceSolid(true);
fPMTFrameLogical->SetVisAttributes(visFrame);

PlacePMTs();

// **FIX 4: Attach Sensitive Detector to PHOTOCATHODE, not glass**
auto* sdManager = G4SDManager::GetSDMpointer();
auto* pmtSD = new PMTSensitiveDetector("PMTSD", "PMTHitsCollection");

// Use PMT model QE curve
auto wl_nm = fPMTModel->GetQEWavelength_nm();
auto qe = fPMTModel->GetQE();
std::vector<double> wl_nm_d(wl_nm.begin(), wl_nm.end());
std::vector<double> qe_d(qe.begin(), qe.end());

// Use SimplePMTResponse with QE for manual detection
SimplePMTResponse response(wl_nm_d, qe_d, 2.0 * ns);
pmtSD->SetPMTResponse(response);
pmtSD->SetVerboseLevel(0); // silenced for speed if want debug do 1

sdManager->AddNewDetector(pmtSD);
fPhotocathodeLogical->SetSensitiveDetector(pmtSD); // ATTACH TO PHOTOCATHODE!

G4cout << "[DetectorConstruction] PMT SensitiveDetector attached to PHOTOCATHODE" << G4endl;
G4cout << "[DetectorConstruction] QE curve loaded with " << wl_nm.size() << " points" << G4endl;

}

void DetectorConstruction::PlacePMTs()
{
const G4double R_inner = fInnerRadius - fPMTPlacementOffset;
const G4double halfHeight = fInnerHeight / 2.0;

const G4int nPhiBarrel = 180;
const G4int nRingBarrel = 35;
const G4double dz = fInnerHeight / (nRingBarrel - 1);

const G4int nCapRings = 17;
const G4double ringSpacing = R_inner / (nCapRings - 1);

G4int copyID = 0;

// **CRITICAL FIX: Place photocathode INSIDE the PMT glass as daughter volume**
// This ensures photons pass through glass first, then hit photocathode
new G4PVPlacement(nullptr,                    // no rotation relative to parent
                  G4ThreeVector(0, 0, 0),     // at center of PMT glass
                  fPhotocathodeLogical,        // the photocathode
                  "Photocathode",              // name
                  fPMTLogical,                 // PARENT is PMT glass!
                  false,                       // no boolean
                  0,                           // copy number (single instance per PMT)
                  false);                      // no overlap check

G4cout << "[PlacePMTs] Photocathode placed INSIDE PMT glass as daughter volume" << G4endl;

// BARREL
for (G4int ir = 0; ir < nRingBarrel; ++ir) {
    const G4double z = -halfHeight + ir * dz;

    for (G4int iphi = 0; iphi < nPhiBarrel; ++iphi) {
        const G4double phi = (2.0 * M_PI * iphi) / nPhiBarrel;

        const G4double x = R_inner * std::cos(phi);
        const G4double y = R_inner * std::sin(phi);
        const G4ThreeVector pos(x, y, z);

        G4ThreeVector centerAtSameZ(0.0, 0.0, z);
        G4ThreeVector dir = (centerAtSameZ - pos);
        if (dir.mag() <= 0.) dir = G4ThreeVector(0.,0.,-1.);
        dir = dir.unit();

        G4ThreeVector zAxis(0., 0., 1.);
        G4ThreeVector rotAxis = zAxis.cross(dir);
        G4double cosA = zAxis.dot(dir);
        if (cosA > 1.0) cosA = 1.0;
        if (cosA < -1.0) cosA = -1.0;
        G4double rotAngle = std::acos(cosA);

        G4RotationMatrix* rot = new G4RotationMatrix();

        if (rotAxis.mag2() > 1e-12) {
            rotAxis = rotAxis.unit();
            rot->rotate(rotAngle, rotAxis);
            rot->rotateY(180.0*deg);
        } else {
            if (cosA < 0.0) {
                rot->rotateX(180.0 * deg);
            }
        }

        // Place PMT glass (photocathode is already inside as daughter)
        new G4PVPlacement(rot, pos, fPMTLogical, "PMT_Barrel",
                          fInnerDetectorLogical, false, copyID, false);
        
        copyID++;
    }
}

// TOP CAP
const G4double zTop = +halfHeight - fPMTRadius;
for (G4int ir = 0; ir < nCapRings; ++ir) {
    const G4double r = ir * ringSpacing;
    const G4int nPhi = std::max(6, static_cast<int>(2.0 * M_PI * r / (2.0 * fPMTRadius)));

    for (G4int iphi = 0; iphi < nPhi; ++iphi) {
        const G4double phi = (2.0 * M_PI * iphi) / nPhi;
        const G4ThreeVector pos(r * std::cos(phi), r * std::sin(phi), zTop);

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateX(180 * deg);

        new G4PVPlacement(rot, pos, fPMTLogical, "PMT_Top",
                          fInnerDetectorLogical, false, copyID, false);
        copyID++;
    }
}

// BOTTOM CAP
const G4double zBot = -halfHeight + fPMTRadius;
for (G4int ir = 0; ir < nCapRings; ++ir) {
    const G4double r = ir * ringSpacing;
    const G4int nPhi = std::max(6, static_cast<int>(2.0 * M_PI * r / (2.0 * fPMTRadius)));

    for (G4int iphi = 0; iphi < nPhi; ++iphi) {
        const G4double phi = (2.0 * M_PI * iphi) / nPhi;
        const G4ThreeVector pos(r * std::cos(phi), r * std::sin(phi), zBot);

        G4RotationMatrix* rot = new G4RotationMatrix();

        new G4PVPlacement(rot, pos, fPMTLogical, "PMT_Bottom",
                          fInnerDetectorLogical, false, copyID, false);
        copyID++;
    }
}

G4cout << "-----------------------------------------------------" << G4endl;
G4cout << " PMTs placed: " << copyID << " (photocathode inside each PMT)" << G4endl;
G4cout << "-----------------------------------------------------" << G4endl;

}

void DetectorConstruction::DefineOpticalSurfaces()
{
// Get QE data from PMT model
auto wl_nm = fPMTModel->GetQEWavelength_nm();
auto qe = fPMTModel->GetQE();


// Convert wavelength to energy (eV)
// IMPORTANT: Wavelengths are in DECREASING order → energies in INCREASING order
const size_t nQE = wl_nm.size();
std::vector<G4double> energy_qe(nQE);
std::vector<G4double> efficiency(nQE);

// Fill arrays in REVERSE order so energies are increasing
for (size_t i = 0; i < nQE; ++i) {
    size_t idx = nQE - 1 - i; // reverse index
    // E = h*c/lambda (higher wavelength = lower energy)
    energy_qe[i] = (h_Planck * c_light) / (wl_nm[idx] * nm);
    efficiency[i] = qe[idx]; // corresponding QE
}

// **FIX 5: Create optical surface with EFFICIENCY property**
fPhotocathodeOpticalSurface = new G4OpticalSurface("PhotocathodeOpticalSurface");
fPhotocathodeOpticalSurface->SetType(dielectric_metal); // photocathode absorbs
fPhotocathodeOpticalSurface->SetModel(unified);
fPhotocathodeOpticalSurface->SetFinish(polished);

auto* photocathodeMPT = new G4MaterialPropertiesTable();
photocathodeMPT->AddProperty("EFFICIENCY", energy_qe.data(), efficiency.data(), nQE);
photocathodeMPT->AddProperty("REFLECTIVITY", energy_qe.data(), std::vector<G4double>(nQE, 0.0).data(), nQE);
fPhotocathodeOpticalSurface->SetMaterialPropertiesTable(photocathodeMPT);

// Apply to photocathode logical volume
new G4LogicalSkinSurface("PhotocathodeSkinSurface", fPhotocathodeLogical, fPhotocathodeOpticalSurface);

G4cout << "[DetectorConstruction] Photocathode optical surface defined with QE curve" << G4endl;

// Tank reflection (steel inner wall)
const G4int n = 5;
G4double energy[n] = {2.00 * eV, 2.40 * eV, 2.75 * eV, 3.10 * eV, 3.50 * eV};

fWaterSteelSurface = new G4OpticalSurface("WaterSteelSurface");
fWaterSteelSurface->SetType(dielectric_metal);
fWaterSteelSurface->SetModel(unified);
fWaterSteelSurface->SetFinish(ground);

G4double reflectivity[n] = {0.20, 0.20, 0.18, 0.15, 0.10};
auto* steelMPT = new G4MaterialPropertiesTable();
steelMPT->AddProperty("REFLECTIVITY", energy, reflectivity, n);
fWaterSteelSurface->SetMaterialPropertiesTable(steelMPT);

new G4LogicalSkinSurface("TankSkinSurface", fTankLogical, fWaterSteelSurface);

}