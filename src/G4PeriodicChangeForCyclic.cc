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
#include "G4DynamicParticle.hh"
#include "G4ParticleChangeForPeriodic.hh"
#include "G4Step.hh"
#include "G4Track.hh"


G4ParticleChangeForPeriodic::G4ParticleChangeForPeriodic() : G4VParticleChange() {}

G4ParticleChangeForPeriodic::~G4ParticleChangeForPeriodic() = default;

G4Step *G4ParticleChangeForPeriodic::UpdateStepForPostStep(G4Step *pStep) {
  G4StepPoint *pPostStepPoint = pStep->GetPostStepPoint();

  pPostStepPoint->SetMomentumDirection(proposedMomentumDirection);
  pPostStepPoint->SetPolarization(proposedPolarization);
  pPostStepPoint->SetPosition(proposedPosition);

  if (isParentWeightProposed) {
    pPostStepPoint->SetWeight(theParentWeight);
  }

  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);

  return pStep;
}

void G4ParticleChangeForPeriodic::AddSecondary(G4DynamicParticle *aParticle) {
  auto *aTrack = new G4Track(aParticle, currentTrack->GetGlobalTime(),
                                currentTrack->GetPosition());

  aTrack->SetTouchableHandle(currentTrack->GetTouchableHandle());

  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChangeForPeriodic::DumpInfo() const {
  G4VParticleChange::DumpInfo();
  G4int oldprc = G4cout.precision(3);

  G4cout << "        Momentum Direction: "
         << std::setw(20) << proposedMomentumDirection
         << G4endl;
  G4cout << "        Polarization: "
         << std::setw(20) << proposedPolarization
         << G4endl;
  G4cout << "        Position: "
         << std::setw(20) << proposedPosition
         << G4endl;
  G4cout.precision(oldprc);
}
