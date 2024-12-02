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
#ifndef G4PeriodicBoundaryProcess_h
#define G4PeriodicBoundaryProcess_h 1

#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4OpticalPhoton.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoTau.hh"
#include "G4ParticleChangeForPeriodic.hh"
#include "G4TransportationManager.hh"
#include "G4VDiscreteProcess.hh"

#include "globals.hh"

enum G4PeriodicBoundaryProcessStatus {
  Undefined,
  Reflection,
  Cycling,
  StepTooSmall,
  NotAtBoundary
};

class G4PeriodicBoundaryProcess : public G4VDiscreteProcess {

public:

  explicit G4PeriodicBoundaryProcess(const G4String &processName = "CycBoundary",
                            G4ProcessType type = fNotDefined, bool per_x = true,
                            bool per_y = true, bool per_z = true, bool ref_walls = false);

  ~G4PeriodicBoundaryProcess() override = default;

  G4PeriodicBoundaryProcess(const G4PeriodicBoundaryProcess &right) = delete;

  G4PeriodicBoundaryProcess &operator=(const G4PeriodicBoundaryProcess &right) = delete;

public:

  G4bool IsApplicable(const G4ParticleDefinition &) override;

  G4double GetMeanFreePath(const G4Track &,
                           G4double,
                           G4ForceCondition *condition) override;
  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'Forced' condition for the DoIt to be invoked at
  // every step. However, only at a boundary will any action be
  // taken.

  G4PeriodicBoundaryProcessStatus GetStatus() const;

  G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &) override;

protected:

  G4ParticleChangeForPeriodic fParticleChange;

private:

  void BoundaryProcessVerbose() const;

  G4PeriodicBoundaryProcessStatus theStatus = Undefined;
  G4ThreeVector OldPosition;
  G4ThreeVector NewPosition;
  G4ThreeVector OldMomentum;
  G4ThreeVector NewMomentum;
  G4ThreeVector OldPolarization;
  G4ThreeVector NewPolarization;
  G4ThreeVector theGlobalNormal;
  G4double kCarTolerance;

  G4bool reflecting_walls;

  bool periodic_x;
  bool periodic_y;
  bool periodic_z;

};

inline G4bool G4PeriodicBoundaryProcess::IsApplicable(const G4ParticleDefinition &
aParticleType) {

  bool applicable = true;

  // do not apply to neutrinos as can lead to very long simulation times if normal 
  // to periodic boundary
  // also do not apply to optical photons as boundary process triggered after PBC

  if (&aParticleType == G4AntiNeutrinoE::AntiNeutrinoE())
    applicable = false;
  else if (&aParticleType == G4NeutrinoE::NeutrinoE())
    applicable = false;
  else if (&aParticleType == G4AntiNeutrinoMu::AntiNeutrinoMu())
    applicable = false;
  else if (&aParticleType == G4NeutrinoMu::NeutrinoMu())
    applicable = false;
  else if (&aParticleType == G4AntiNeutrinoTau::AntiNeutrinoTau())
    applicable = false;
  else if (&aParticleType == G4NeutrinoTau::NeutrinoTau())
    applicable = false;
  else if (&aParticleType == G4OpticalPhoton::OpticalPhoton())
    applicable = false;

  return applicable;

}

inline G4PeriodicBoundaryProcessStatus G4PeriodicBoundaryProcess::GetStatus() const {
  return theStatus;
}

#endif
