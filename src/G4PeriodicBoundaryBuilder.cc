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
#include "G4PeriodicBoundaryBuilder.hh"
#include "G4LogicalVolumePeriodic.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

G4LogicalVolume *G4PeriodicBoundaryBuilder::Construct(G4LogicalVolume *logical_world) {
  G4cout << "Construction de la geometrie aux bords periodique de G4PeriodicBoundary" << G4endl; // Ajout pour savoir si le construct marche bien (Mathieu)
  auto *world = (G4Box *) logical_world->GetSolid();

  double buffer = 0.01 * micrometer;

  double periodic_world_hx = world->GetXHalfLength();
  double periodic_world_hy = world->GetYHalfLength();
  double periodic_world_hz = world->GetZHalfLength();

  /*reset the size of the world volume to be slightly larger to accommodate the
  the cyclic world volume without sharing a surface*/
  world->SetXHalfLength(world->GetXHalfLength() + buffer);
  world->SetYHalfLength(world->GetYHalfLength() + buffer);
  world->SetZHalfLength(world->GetZHalfLength() + buffer);

  auto *periodic_world = new G4Box("cyclic", periodic_world_hx, periodic_world_hy,
                                    periodic_world_hz);

  logical_periodic = new G4LogicalVolumePeriodic(periodic_world,
                                                 logical_world->GetMaterial(), "logical_periodic");

  logical_periodic->SetVisAttributes(G4Color::Magenta());

  new G4PVPlacement(nullptr, G4ThreeVector(), logical_periodic, "physical_cyclic",
                    logical_world, false, 0, true); //check for overlaps

  return logical_periodic;
}
