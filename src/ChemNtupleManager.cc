//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file ChemNtupleManager.cc
/// \brief Implementation of the ChemNtupleManager class

#include "ChemNtupleManager.hh"

#include "G4UnitsTable.hh"
#include "G4Filesystem.hh"

ChemNtupleManager* ChemNtupleManager::fInstance = nullptr; //Mathieu
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ChemNtupleManager* ChemNtupleManager::Instance(){ //Mathieu

    if (fInstance == nullptr) {

        static ChemNtupleManager parparser;
        fInstance = &parparser;

    }

    return fInstance;

}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::Book()
{
   





    // open output file
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4bool fileOpen = analysisManager->OpenFile(fFileName.c_str());
    if (!fileOpen) {
        G4cout << "\n---> HistoManager::book(): cannot open " << fFileName
               << G4endl;
        return;
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::Save()
{  
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::FillNtupleIColumn(G4int icol, G4int ival)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(icol,ival);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::FillNtupleFColumn(G4int icol, G4float fval)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleFColumn(icol,fval);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::FillNtupleDColumn(G4int id, G4int icol, G4double dval)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(id,icol,dval);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::AddNtupleRow()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::AddNtupleRow(G4int id)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->AddNtupleRow(id);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemNtupleManager::CreateNtuples()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetDefaultFileType("root");
    // Create directories
    const G4String directoryName = "ntuple";
    //if (analysisManager->GetNtupleDirectoryName() != directoryName) { //Mathieu
        analysisManager->SetNtupleDirectoryName(directoryName);
    //}
    // create ntuple

    if (analysisManager->GetFirstNtupleId() != 1) analysisManager->SetFirstNtupleId(1);
    
    // DBScan

    analysisManager->CreateNtuple("ntuple_2","DB_chemical_stage");
    analysisManager->CreateNtupleDColumn(1,"strand");
    analysisManager->CreateNtupleDColumn(1,"copyNumber");
    //analysisManager->CreateNtupleDColumn(1,"xp"); //Mathieu
    //analysisManager->CreateNtupleDColumn(1,"yp");
    //analysisManager->CreateNtupleDColumn(1,"zp");
    analysisManager->CreateNtupleDColumn(1,"time");
    analysisManager->CreateNtupleDColumn(1,"base");
    analysisManager->CreateNtupleDColumn(1,"eventNumber"); // Mathieu
    analysisManager->CreateNtupleDColumn(1,"voxelCopyNumber"); // Mathieu

    analysisManager->FinishNtuple(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
