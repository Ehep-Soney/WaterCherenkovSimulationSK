// DetectorConstruction.hh - FIXED VERSION
// Added photocathode volume and optical surface members
#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>
#include <memory>

#include "PMTObject.hh"

class G4OpticalSurface;
class G4MaterialPropertiesTable;
class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;

/**
 * @class DetectorConstruction
 * @brief Defines a Super-K–like Cherenkov detector geometry with Gd-doped water and PMTs.
 * 
 * FIXED VERSION: Now includes photocathode volumes for proper optical photon detection
 */
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    
    // Mandatory user method
    G4VPhysicalVolume* Construct() override;

private:
    // --- Helper methods ---
    void DefineMaterials();
    void BuildGeometry();
    void PlacePMTs();
    void DefineOpticalSurfaces();
    
    // --- Volumes ---
    G4LogicalVolume* fWorldLogical = nullptr;
    G4LogicalVolume* fTankLogical = nullptr;
    G4LogicalVolume* fInnerDetectorLogical = nullptr;
    G4LogicalVolume* fPMTLogical = nullptr;
    G4LogicalVolume* fPhotocathodeLogical = nullptr;  // NEW: photocathode volume
    G4LogicalVolume* fPMTFrameLogical = nullptr;
    
    G4VPhysicalVolume* fWorldPhysical = nullptr;
    G4VPhysicalVolume* fTankPhysical = nullptr;
    G4VPhysicalVolume* fInnerDetectorPhysical = nullptr;
    
    // --- Materials ---
    G4Material* fGdDopedWater = nullptr;
    G4Material* fSteel = nullptr;
    G4Material* fPMTGlass = nullptr;
    G4Material* fPhotocathodeMaterial = nullptr;  // NEW: photocathode material
    
    // --- Optical surfaces ---
    G4OpticalSurface* fWaterSteelSurface = nullptr;
    G4OpticalSurface* fPhotocathodeOpticalSurface = nullptr;  // NEW: photocathode optical surface with QE
    G4OpticalSurface* fWaterWorldSurface;
    
    // --- Visualization attributes ---
    G4VisAttributes* fVisWorld = nullptr;
    G4VisAttributes* fVisTank = nullptr;
    G4VisAttributes* fVisInnerDetector = nullptr;
    G4VisAttributes* fVisPMT = nullptr;
    G4VisAttributes* fVisPhotocathode = nullptr;  // NEW: photocathode visualization
    
    // --- Geometry parameters (Super-K–like) ---
    G4double fTankOuterRadius = 16.9 * m;     // Outer tank radius
    G4double fTankHeight = 36.2 * m;          // Total height
    G4double fInnerRadius = 16.0 * m;         // Inner detector radius
    G4double fInnerHeight = 33.8 * m;         // Inner detector height
    G4int    fNumPMTs = 11146;                // Total number of PMTs
    
    // --- PMT-specific parameters ---
    G4double fPMTRadius = 0.254 * m;          // ~20-inch PMT
    G4double fPMTGlassThickness = 0.4 * cm;
    G4double fPMTPlacementOffset = 1.0 * cm;  // offset inward from wall
    G4double fGdMassPercent = 0.2 * perCent;  // Gadolinium doping fraction
    
    // --- PMT model (quantum efficiency, dimensions, etc.) ---
    std::unique_ptr<CherenkovSim::PMT::PMT20inch> fPMTModel;
    
    // --- Placement bookkeeping ---
    std::vector<G4ThreeVector> fPMTPositions;
};