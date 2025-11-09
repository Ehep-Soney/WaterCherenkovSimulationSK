// src/EventAction.cc - NEW ROOT-ENABLED VERSION
#include "EventAction.hh"
#include "DataManager.hh" // <-- Include our new DataManager
#include "PMTHit.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// The constructor is where we get the pointer to our DataManager instance.
EventAction::EventAction() 
  : G4UserEventAction(), 
    fLastEventPE(0),
    fDataManager(DataManager::GetInstance()) // Get the singleton instance
{}

// At the beginning of an event, we reset our counters AND clear our data vectors.
void EventAction::BeginOfEventAction(const G4Event*) {
    fLastEventPE = 0;
    
    // Clear the vectors to ensure no data from the previous event remains
    fPmtID.clear();
    fPmtX.clear();
    fPmtY.clear();
    fPmtZ.clear();
    fPmtNumPEs.clear();
}

// At the end of an event, we collect data and then pass it to the DataManager.
void EventAction::EndOfEventAction(const G4Event* event) {
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    if (!HCE) return;

    G4SDManager* sdman = G4SDManager::GetSDMpointer();
    G4int hcID = sdman->GetCollectionID("PMTHitsCollection");
    if (hcID < 0) return;

    auto* hits = static_cast<PMTHitsCollection*>(HCE->GetHC(hcID));
    if (!hits) return;

    // Loop through all hits to collect data and calculate the total PE for the summary.
    for (size_t i = 0; i < hits->entries(); ++i) {
        PMTHit* hit = (*hits)[i];
        if (!hit) continue;

        // Add this hit's PE count to the event total
        fLastEventPE += hit->GetTotalPE();

        // --- NEW DATA COLLECTION LOGIC ---
        // Instead of writing to a file, we store the data in our member vectors.
        fPmtID.push_back(hit->GetPMTID());
        fPmtX.push_back(hit->GetPosition().x() / m);
        fPmtY.push_back(hit->GetPosition().y() / m);
        fPmtZ.push_back(hit->GetPosition().z() / m);
        fPmtNumPEs.push_back(hit->GetTotalPE());
    }

    // --- NEW DATA SAVING STEP ---
    // After collecting all the hits, make ONE call to the DataManager to fill the TTree.
    if (fDataManager) {
        fDataManager->FillEventData(event->GetEventID(), fPmtID, fPmtX, fPmtY, fPmtZ, fPmtNumPEs);
    }

    // Visualization code remains the same
    G4VVisManager* visman = G4VVisManager::GetConcreteInstance();
    if (visman) {
        for (size_t i = 0; i < hits->entries(); ++i) {
            PMTHit* h = (*hits)[i];
            int npe = h->GetTotalPE();
            if (npe <= 0) continue;

            G4double screenSize = std::min(20.0, 2.0 + std::sqrt((double)npe) * 2.0);
            G4Circle circle(h->GetPosition());
            circle.SetScreenSize(screenSize);
            circle.SetFillStyle(G4Circle::filled);

            G4Colour colour;
            if (npe < 3) colour = G4Colour(0.2, 0.4, 1.0);
            else if (npe < 8) colour = G4Colour(1.0, 1.0, 0.0);
            else colour = G4Colour(1.0, 0.2, 0.0);

            circle.SetVisAttributes(G4VisAttributes(colour));
            visman->Draw(circle);
        }
    }

    // Print the summary for the event
    G4cout << "Event " << event->GetEventID() << " total PEs = " << fLastEventPE << G4endl;
}