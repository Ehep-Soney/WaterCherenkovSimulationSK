// EventAction.cc
#include "EventAction.hh"
#include "PMTHit.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"

EventAction::EventAction() : G4UserEventAction(), fLastEventPE(0) {}

void EventAction::BeginOfEventAction(const G4Event*) {
fLastEventPE = 0;
}

void EventAction::EndOfEventAction(const G4Event* event) {
fLastEventPE = 0;

G4HCofThisEvent* HCE = event->GetHCofThisEvent();
if (!HCE) return;

G4SDManager* sdman = G4SDManager::GetSDMpointer();
// Ensure the collection name matches the one used in PMTSensitiveDetector
G4int hcID = sdman->GetCollectionID("PMTHitsCollection");
if (hcID < 0) return;

PMTHitsCollection* hits = (PMTHitsCollection*)(HCE->GetHC(hcID));
if (!hits) return;

// Sum PE and (optionally) draw visualization markers
for (size_t i = 0; i < hits->entries(); ++i) {
    PMTHit* h = (*hits)[i];
    if (!h) continue;
    fLastEventPE += h->GetTotalPE();
}

// Visualization: draw circles at PMT positions scaled by PE
G4VVisManager* visman = G4VVisManager::GetConcreteInstance();
if (visman) {
    for (size_t i = 0; i < hits->entries(); ++i) {
        PMTHit* h = (*hits)[i];
        if (!h) continue;
        int npe = h->GetTotalPE();
        if (npe <= 0) continue;

        // screen size proportional to sqrt(npe) (tunable)
        G4double screenSize = std::min(20.0, 2.0 + std::sqrt((double)npe) * 2.0);

        G4Circle circle(h->GetPosition());
        circle.SetScreenSize(screenSize);
        circle.SetFillStyle(G4Circle::filled);

        // color map: small npe -> blue, medium -> yellow, large -> red
        G4Colour colour;
        if (npe < 3) colour = G4Colour(0.2, 0.4, 1.0);        // blue
        else if (npe < 8) colour = G4Colour(1.0, 1.0, 0.0);   // yellow
        else colour = G4Colour(1.0, 0.2, 0.0);               // red

        circle.SetVisAttributes(G4VisAttributes(colour));
        visman->Draw(circle);
    }
}

// Optional: print summary
G4cout << "Event " << event->GetEventID() << " total PEs = " << fLastEventPE << G4endl;

}