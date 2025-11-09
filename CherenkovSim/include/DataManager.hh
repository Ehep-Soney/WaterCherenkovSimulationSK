// include/DataManager.hh
#pragma once

#include <string>
#include <vector>

// Forward declarations for ROOT classes. This is a C++ trick to speed up
// compilation by not needing to include the full ROOT headers here.
class TFile;
class TTree;

// DataManager is a "singleton" class. This is a design pattern that ensures
// we only ever have one single instance of it, which is perfect for managing
// a single output file for the whole application.
class DataManager {
public:
    // This is the static method to get the one and only instance of the class.
    static DataManager* GetInstance();

    // These are the public functions we will call from other parts of our code.
    void OpenFile(const std::string& filename);
    void CloseFile();
    void FillEventData(int eventID, const std::vector<int>& pmtIDs, 
                       const std::vector<double>& pmtXs, const std::vector<double>& pmtYs, 
                       const std::vector<double>& pmtZs, const std::vector<int>& pmtPEs);
    void ClearEventData();

private:
    // The constructor and destructor are private to ensure no one else can
    // create or destroy an instance. Only GetInstance() can.
    DataManager();
    ~DataManager();

    // --- Member Variables ---

    // Pointers to our ROOT objects.
    TFile* fRootFile = nullptr;
    TTree* fTree = nullptr;

    // These variables will hold the data for a single event.
    // They will be connected to the "branches" of our TTree.
    int fEventID;
    std::vector<int>    fPmtID;
    std::vector<double> fPmtX;
    std::vector<double> fPmtY;
    std::vector<double> fPmtZ;
    std::vector<int>    fPmtNumPEs;
};