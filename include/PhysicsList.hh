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
//
/// \file PhysicsList.hh
/// \brief Definition of the PhysicsList class

#ifndef PhysicsList_h
#define PhysicsList_h 1
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <memory>

#include "G4VPhysicsConstructor.hh"
#include "PhysicsMessenger.hh"
#include "G4PeriodicBoundaryProcess.hh" // ajout de ce include afin de retrouver la class G4PeriodicBoundaryProcess (Mathieu)


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4PeriodicBoundaryProcess; // Appel de G4PeriodicBoundaryProcess (Mathieu)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:
    PhysicsList(const G4String &name = "Periodic", bool per_x = true, bool per_y = true, bool per_z = false, bool ref_walls = false, bool _activate_pbc = false);
    // Modification du constructeur de class en ajoutant les variables name, per_x, per_y, per_z, ref_walls, _activate_pbc (Mathieu)
    ~PhysicsList() override = default;
    
    void ConstructParticle() override;
    void ConstructProcess() override;

    void RegisterPhysicsList(const G4String& name);

    void apply_pbc_electrons(auto aParticleIterator, auto pbc_process); // declaration de fonction "publique" nommee apply_pbc_electrons (Mathieu)
private:
    std::unique_ptr<G4VPhysicsConstructor>    fDNAPhysicsList{nullptr};
    std::unique_ptr<PhysicsMessenger>         fPhysMsg{nullptr};
    G4String                                  fPhysDNAName{""};

    // declaration des 5 variables prives reflecting_walls, periodic_x, periodic_y, periodic_z, activate_pbc (Mathieu)

    bool reflecting_walls;
    bool periodic_x, periodic_y, periodic_z;
    bool activate_pbc;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif