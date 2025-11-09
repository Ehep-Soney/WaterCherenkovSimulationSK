// PMTHit.hh
#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <iostream>

class PMTHit : public G4VHit {
public:
    PMTHit();
    explicit PMTHit(int pmtID);
    ~PMTHit() override = default;

    // basic setters
    void SetPMTID(int id) { fPMTID = id; }
    void SetPosition(const G4ThreeVector& pos) { fPosition = pos; }

    // add 1 photoelectron at time t (ns)
    void AddPE(double time_ns) {
        ++fTotalPE;
        fTimes.push_back(time_ns);
    }

    // getters
    int GetPMTID() const { return fPMTID; }
    const G4ThreeVector& GetPosition() const { return fPosition; }
    int GetTotalPE() const { return fTotalPE; }
    const std::vector<double>& GetTimes() const { return fTimes; }

    void Print() const;

private:
    int fPMTID;
    int fTotalPE;
    G4ThreeVector fPosition;
    std::vector<double> fTimes;
};

typedef G4THitsCollection<PMTHit> PMTHitsCollection;