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
#include "G4PeriodicBoundaryPhysics.hh"
#include "G4PeriodicBoundaryProcess.hh"

#include "G4PhysicsConstructorFactory.hh"
#include "G4ProcessManager.hh"

#include "globals.hh"

#include "PhysRunAction.hh" // remplacement de RunAction.hh avec PhysRunAction.hh (Mathieu)

G4_DECLARE_PHYSCONSTR_FACTORY(G4PeriodicBoundaryPhysics);

G4PeriodicBoundaryPhysics::G4PeriodicBoundaryPhysics(const G4String &name,
                                                     bool per_x, bool per_y, bool per_z, bool ref_walls, bool _activate_pbc)
    : G4VPhysicsConstructor(name) {

  verboseLevel = 0;
  periodic_x = per_x;
  periodic_y = per_y;
  periodic_z = per_z;
  reflecting_walls = ref_walls;
  activate_pbc = _activate_pbc;
}

G4PeriodicBoundaryPhysics::~G4PeriodicBoundaryPhysics() = default;

void G4PeriodicBoundaryPhysics::ConstructParticle() {}

void G4PeriodicBoundaryPhysics::ConstructProcess() {

  if (verboseLevel > 0)
    G4cout << "Constructing cyclic boundary physics process" << G4endl;

  auto *pbc = new G4PeriodicBoundaryProcess("Cyclic",
                                                                 fNotDefined, periodic_x, periodic_y, periodic_z,
                                                                 reflecting_walls);

  if (verboseLevel > 0) pbc->SetVerboseLevel(verboseLevel);

  auto aParticleIterator = GetParticleIterator();

  if (activate_pbc)
  {
  G4cout << "Activating pbc" << G4endl;
  apply_pbc_electrons(aParticleIterator, pbc);
  }

}

void G4PeriodicBoundaryPhysics::apply_pbc_electrons(auto aParticleIterator, auto pbc_process)
{ aParticleIterator->reset();

  G4ProcessManager *pManager = nullptr;

  while ((*aParticleIterator)())
  {

    G4ParticleDefinition *particle = aParticleIterator->value();

    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();

    if (!pManager)
    {
      std::ostringstream o;
      o << "Particle " << particleName << "without a Process Manager";
      G4Exception("G4PeriodicBoundaryPhysics::ConstructProcess()", "",
                  FatalException, o.str().c_str());
      return;
    }

    if (pbc_process->IsApplicable(*particle) and particleName == "e-")
    {
      // if (verboseLevel > 0)
        G4cout << "Adding pbc to " << particleName << G4endl;
      pManager->AddDiscreteProcess(pbc_process);
    }
  }
}
