// make_sk_plot.C
// This is the final script. It analyzes the "light bulb" data and generates
// a plot of Average Charge vs. Photon Path Length.
// VERSION: Updated to use a physically-correct exponential fit.

#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h" // Need to include TF1 header for the fit function
#include "iostream"
#include <cmath>
#include <vector>

void make_sk_plot2() {
    // ================================================================
    // ===                 CONFIGURE YOUR RUNS                  ===
    // ================================================================
    // filenames and their corresponding z_source positions
    std::vector<std::string> filenames = {"hits_z0.root", "hits_z4.root", "hits_z8.root", "hits_z10.root", "hits_z12.root", "hits_z14.root", "hits_z15.root", "hits_z16.root"};
    std::vector<double> z_positions_m = {0.0, 4.0, 8.0, 10.0, 12.0, 14.0, 15.0, 16.0};
    // ================================================================

    if (filenames.size() != z_positions_m.size()) {
        std::cerr << "Error: Number of files and Z positions do not match!" << std::endl;
        return;
    }

    // --- CREATE THE HISTOGRAM (A TProfile) ---
    // TProfile is perfect for this: for each bin on the X-axis (path length),
    // it automatically calculates the average of the Y-values (charge).
    TProfile* profile = new TProfile("profile", "Effective Observed Charge vs. Path Length;Path Length l (m);Average Charge Q (PEs)",
                                     50, 0, 25); // 50 bins from 0 to 25 meters

    // --- LOOP OVER ALL FILES ---
    for (size_t i = 0; i < filenames.size(); ++i) {
        std::string filename = filenames[i];
        double z_source = z_positions_m[i];

        TFile* inputFile = TFile::Open(filename.c_str());
        if (!inputFile || inputFile->IsZombie()) continue;
        TTree* tree = (TTree*)inputFile->Get("hits");
        if (!tree) continue;
        
        std::cout << "Processing file: " << filename << " (source at z=" << z_source << "m)" << std::endl;

        // Connect variables to the TTree
        std::vector<double>* pmtX = nullptr;
        std::vector<double>* pmtY = nullptr;
        std::vector<double>* pmtZ = nullptr;
        std::vector<int>*    numPEs = nullptr;
        tree->SetBranchAddress("pmtX", &pmtX);
        tree->SetBranchAddress("pmtY", &pmtY);
        tree->SetBranchAddress("pmtZ", &pmtZ);
        tree->SetBranchAddress("numPEs", &numPEs);

        // --- Loop over events and hits in this file ---
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            for (size_t hit = 0; hit < pmtX->size(); ++hit) {
                
                // THE CRUCIAL CALCULATION: TRUE 3D PHOTON PATH LENGTH
                double dx = pmtX->at(hit);
                double dy = pmtY->at(hit);
                double dz = pmtZ->at(hit) - z_source;
                double pathLength = std::sqrt(dx*dx + dy*dy + dz*dz);

                int charge = numPEs->at(hit);

                // Fill the profile. ROOT does the averaging for us.
                profile->Fill(pathLength, charge);
            }
        }
        inputFile->Close();
    }
    std::cout << "All files processed." << std::endl;

    // =========================================================================
    // ===                 PLOTTING & FITTING (UPDATED SECTION)              ===
    // =========================================================================
    TCanvas* c1 = new TCanvas("c1", "SK Attenuation Plot", 800, 600);
    gStyle->SetOptStat(0);   // Turn off the default statistics box
    gStyle->SetOptFit(1111); // Display the fit results box on the plot

    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(1.0);
    profile->SetLineColor(kBlack); // Ensure data points are black
    profile->Draw("E1");           // Draw with error bars

    // 2. Define the EXPONENTIAL fitting function based on the Beer-Lambert law.
    //    The physics model for attenuation is Q(l) = Q_0 * exp(-l / lambda)
    //    In ROOT's notation:
    //    - "x" is the path length l
    //    - "[0]" is the initial charge Q_0
    //    - "[1]" is the attenuation length lambda
    TF1 *expFit = new TF1("expFit", "[0]*exp(-x/[1])", 7.0, 22.0); // Fit over a clean range of the data

    // 3. Give the fit some reasonable starting guesses for the parameters.
    //    This helps the algorithm find the best fit.
    expFit->SetParameter(0, 1060); // Initial guess for Q_0, based on the y-value of your data at the start
    expFit->SetParameter(1, 80.0);   // Initial guess for the attenuation length (should be a large, positive number for pure water)
    
    // 4. (Optional but good practice) Label the parameters in the fit box.
    expFit->SetParName(0, "Q_0 (PEs)");
    expFit->SetParName(1, "Lambda (m)");
    
    // 5. Change the line color to make it stand out.
    expFit->SetLineColor(kRed);

    // 6. Perform the fit. The "R" tells it to use the range we specified in the TF1 constructor.
    //    The "+" tells ROOT to add this function to its internal list, which prevents it from being deleted.
    profile->Fit(expFit, "R+");

    // 7. Save the final plot.
    c1->SaveAs("final_sk_plot_exponential_fit2.png");
    std::cout << "Plot saved as final_sk_plot_exponential_fit2.png" << std::endl;
}
