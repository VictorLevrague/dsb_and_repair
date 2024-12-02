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
//redistributed by https://github.com/amentumspace/g4pbc
#pragma once

#include "G4VPhysicsConstructor.hh"

class G4PeriodicBoundaryProcess;

class G4PeriodicBoundaryPhysics : public G4VPhysicsConstructor {

public:

  explicit G4PeriodicBoundaryPhysics(const G4String &name = "Periodic", bool per_x = true,
                            bool per_y = true, bool per_z = false, bool ref_walls = false,
                            bool _activate_pbc = false);

  ~G4PeriodicBoundaryPhysics() override;

  void apply_pbc_electrons(auto aParticleIterator, auto pbc_process);

protected:

   void ConstructParticle() override;

   void ConstructProcess() override;

private:
  bool reflecting_walls;
  bool periodic_x, periodic_y, periodic_z;
  bool activate_pbc;

};
