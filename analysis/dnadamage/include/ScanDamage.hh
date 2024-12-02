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
/// \file ScanDamage.hh
/// \brief Definition of the ScanDamage class

#ifndef ScanDamage_h
#define ScanDamage_h

#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "Damage.hh"
#include <filesystem>
namespace fs = std::filesystem;
class TFile;
using ullint = unsigned long long int;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

struct VoxelData
{
    VoxelData(std::string voxelName,
	    int chromo, int domain, ullint firstBpCN,
	    double x, double y, double z,
	    double xx, double xy, double xz,
	    double yx, double yy, double yz,
	    double zx, double zy, double zz) //Floriane
  
    {
        fVoxelName = voxelName;//Mathieu
        fChromosome = chromo;
        fDomain = domain;
        fFirstBpCopyNum = firstBpCN;

        fx = x;//Mathieu
        fy = y;
        fz = z;

        // We need the inverse (Mathieu)
    //
        /*double det = xx * (yy * zz - zy * yz) - xy * (yx * zz - yz * zx) + xz * (yx * zy - yy * zx);
        double invdet = 1 / det;
        fxx = (yy * zz - zy * yz) * invdet;
        fxy = (xz * zy - xy * zz) * invdet;
        fxz = (xy * yz - xz * yy) * invdet;
        fyx = (yz * zx - yx * zz) * invdet;
        fyy = (xx * zz - xz * zx) * invdet;
        fyz = (yx * xz - xx * yz) * invdet;
        fzx = (yx * zy - zx * yy) * invdet;
        fzy = (zx * xy - xx * zy) * invdet;
        fzz = (xx * yy - yx * xy) * invdet;*/

        fxx = xx;
        fxy = xy;
        fxz = xz;
        fyx = yx;
        fyy = yy;
        fyz = yz;
        fzx = zx;
        fzy = zy;
        fzz = zz;



    }

    ~VoxelData() {}

    std::string fVoxelName{""};//Mathieu
    int fChromosome{0};
    int fDomain{0};
    ullint fFirstBpCopyNum{0};

    double fx{0.};//Mathieu
    double fy{0.};
    double fz{0.};

    double fxx{0.};//Mathieu
    double fxy{0.};
    double fxz{0.};
    double fyx{0.};
    double fyy{0.};
    double fyz{0.};
    double fzx{0.};
    double fzy{0.};
    double fzz{0.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

typedef std::vector<std::vector<ullint> > Table;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class Trier {
public:
    bool operator()(const std::vector<ullint>& a, const std::vector<ullint>& b)
    {
        bool bb = false;
        if(a[0] < b[0]) bb = true;
        return bb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ScanDamage
{
public:
    ScanDamage();
    ~ScanDamage() = default;  
    std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > ExtractDamage();
    void SetThresholdEnergy(double e) {fThresholdEnergy = e;}
    void SetPathToOutputs(std::string path){fPathToOutputs = path;} //Mathieu
    void SetProbabilityForIndirectSBSelection(double p) {fProbabilityForIndirectSB = p;}
    double GetThresholdEnergy() {return fThresholdEnergy;}
    double GetProbabilityForIndirectSBSelection() {return fProbabilityForIndirectSB;}
    std::map<int, Table> GetMergedSBData() {return fMergedTables;}
    double GetEdepSumInNucleus() {return fEdepSumInNucleus;} //eV
    double GetTotalNbBpPlacedInGeo() {return fTotalNbBpPlacedInGeo;}
    double GetTotalNbHistonePlacedInGeo() {return fTotalNbHistonePlacedInGeo;}
    double GetNucleusVolume() {return fNucleusVolume;} 
    double GetNucleusMassDensity() {return fNucleusMassDensity;}
    double GetNucleusMass() {return fNucleusMass;}
    std::map<int,ullint> GetChromosomeBpSizesMap() {return fChromosomeBpMap;}
    void SkipScanningIndirectDamage() {fSkipScanningIndirectDamage = true;}
    bool SkippedScanningIndirectDamage() {return fSkipScanningIndirectDamage;}
private:
    void ScanDamageFromPhys();
    void ScanDamageFromChem();
    void RetrieveVoxelBp();
    void FillVoxelData();
    void AnaPhysRootFile(const std::string fileName);
    void AnaChemRootFile(fs::directory_entry entry);
    void AnaPhysRootTree1(TFile*);
    void AnaPhysRootTree2(TFile*);
    void SortPhysTableWithSelection();
    void SortChemTableWithSelection();
    void ReadCellandVoxelDefFilePaths();
    void MergeDamageFromPhysChem();
    std::tuple<unsigned int, unsigned int> GetEventNberAndVoxelNberFromChemRoot(const std::string fileNam);
    double fThresholdEnergy{17.5};//eV
    std::string fPathToOutputs{""};//Mathieu
    std::string fCellDefFilePath{""};
    std::set<std::string> fVoxelDefFilesList;
    std::map<std::string, int> fBpPerVoxel;
    std::vector<VoxelData> fVoxels;
    std::map<int, Table> fphysTables, fphysSlectedTables, fchemTables, fchemSlectedTables,fMergedTables;
    std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > fDamage;
    double fEdepSumInNucleus{0}; //eV
    int corruptedFiles = 0;
    double fTotalNbBpPlacedInGeo{0};
    double fTotalNbHistonePlacedInGeo{0};
    double fNucleusVolume{0};
    double fNucleusMassDensity{0};
    double fNucleusMass{0};
    double fProbabilityForIndirectSB{0.4};
    std::map<int,ullint> fChromosomeBpMap; //Store number of Bp in each Chomosomes;
    bool fSkipScanningIndirectDamage{false};

    double fXCell{0};//Mathieu
    double fYCell{0};//Mathieu
    double fZCell{0};//Mathieu

    //---coordinates of DNA
  std::map<std::string, std::map<int, std::vector<double> > > f_ph1_coord;   //Floriane
  std::map<std::string, std::map<int, std::vector<double> > > f_ph2_coord;   //Floriane
  std::map<std::string, std::map<int, std::vector<double> > > f_deox1_coord; //Floriane
  std::map<std::string, std::map<int, std::vector<double> > > f_deox2_coord; //Floriane
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


