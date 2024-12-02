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

#include "G4LogicalVolume.hh"

/*a derived class of logical volume that is used to indicate whether
periodic boundary conditions should be applied to this volume*/

class G4LogicalVolumePeriodic : public G4LogicalVolume {
public:
  G4LogicalVolumePeriodic(G4VSolid *pSolid,
                          G4Material *pMaterial,
                          const G4String &name) :
      G4LogicalVolume(pSolid, pMaterial, name) {
    //
  };

  ~G4LogicalVolumePeriodic() override = default;

  G4bool IsExtended() const override { return true; }
  // Return true if it is not a base-class object.


};

