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
/// \file ScanDamage.cc
/// \brief Implementation of the ScanDamage class

#include "ScanDamage.hh"

#include "TSystemDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
TRandom gRandomGen; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScanDamage::ScanDamage(): fSkipScanningIndirectDamage(false)
{
    fThresholdEnergy = 17.5; //eV
    fEdepSumInNucleus = 0;
    fProbabilityForIndirectSB = 0.40; // 40%
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > ScanDamage::ExtractDamage(){

    // Read imp.info
    ReadCellandVoxelDefFilePaths();

    // clearing vectors
    fMergedTables.clear();
    fDamage.clear();

    // Go over voxel type file 
    // Stores the number of base pair of each voxel type in the map fBpPerVoxel
    RetrieveVoxelBp();

    // Go over the nucleus voxel file
    FillVoxelData();

    //----> Calls 
    //      (1) AnaPhysRootFile
    //         (1.1) AnaPhysRootTree1 
    //         (1.2) AnaPhysRootTree2
    //      (2) SortPhysTableWithSelection
    ScanDamageFromPhys();

    //-----> Calls
    // (1) AnaChemRootFile
    // (2) SortChemTableWithSelection
    if (!fSkipScanningIndirectDamage) ScanDamageFromChem();

    // clearing coordinates vectors
    f_ph1_coord.clear(); // netoyage du vecteur dans lequel les coordonnées de phosphate et deoxyriboses sont enregistrés (Mathieu)
    f_ph2_coord.clear();
    f_deox1_coord.clear();
    f_deox2_coord.clear();

    fVoxels.clear();

    // Put the list of damages from phys and chem into one table fMergedTables 
    MergeDamageFromPhysChem();



    for (const auto& [pChrom,table] : fMergedTables) {
        unsigned int evt;
        unsigned int strand;
        ullint cpyNb;
        unsigned int isBase;
        double time;
        double edep;
        double origin;

        unsigned long int domain;
        double x,y,z;

        Damage::DamageType pType;
        Damage::DamageCause pOrigin;
        std::map<unsigned int,std::vector<Damage> > perChromoDamage;
        for (const auto& v : table) {
            evt = v.at(0);
            strand = v.at(1);
            cpyNb = v.at(2);
            isBase = v.at(3);
            time = v.at(4);
            edep = v.at(5);
            origin = v.at(6);

            domain = v.at(7);
            x = v.at(8)/1000. - fXCell;
            y = v.at(9)/1000. - fYCell;
            z = v.at(10)/1000. - fZCell;

            std::cout << "-----------------" << std::endl;
            std::cout << "Event: " << evt 
            << "\nstrand: " << strand << ", copy number (Mbp): " << cpyNb/1e6 << ", is base: " << isBase << ", origin: " << origin
            << "\ntime: " << time << ", edep: " << edep 
            << "\n Chr:" << pChrom << ", TAD: " << domain
            << "\nx: " << x << ", y: " << y << ", z:" << z << std::endl;

            if(isBase==0)
                pType=Damage::DamageType::fBackbone;
            else
                pType=Damage::DamageType::fBase;
            if(origin==0)
                pOrigin=Damage::DamageCause::fDirect;
            else
                pOrigin=Damage::DamageCause::fIndirect;
            if (perChromoDamage.find(evt) == perChromoDamage.end()) {
                std::vector<Damage> dm{Damage(pType,pChrom,evt,strand,cpyNb,domain,Position(x,y,z),
                                            pOrigin,Damage::DamageChromatin::fUnspecified)};
                perChromoDamage.insert({evt,dm});
            } else perChromoDamage[evt].push_back(Damage(pType,pChrom,evt,strand,cpyNb,domain,Position(x,y,z),
                                            pOrigin,Damage::DamageChromatin::fUnspecified));
        }
        fDamage.insert({pChrom,perChromoDamage});
    }
    return fDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::RetrieveVoxelBp()
{
    fBpPerVoxel.clear();

    f_ph1_coord.clear(); // Mathieu
    f_ph2_coord.clear();
    f_deox1_coord.clear();
    f_deox2_coord.clear();
    
    for (const auto & entry : fVoxelDefFilesList) {
        std::ifstream file(entry);
        if(!file.good() )
        {
            std::cerr<<"**** Fatal Error *****"<<std::endl;
            std::cerr<<"ScanDamage::RetrieveVoxelBp: No file named "<<entry<<std::endl;
            std::cerr<<"*************** *****"<<std::endl;
            exit(EXIT_FAILURE);
        }
        std::string voxelName = "noName";
        std::string line;
        bool foundName = false;
        bool foundNumOfBp = false;

        std::map<int, std::vector<double>> ph1_coord; // Mathieu
        std::map<int, std::vector<double>> ph2_coord;
        std::map<int, std::vector<double>> deox1_coord;
        std::map<int, std::vector<double>> deox2_coord;

        while(std::getline(file, line))
            // && !foundNumOfBp) // Supprimer cette condition (Mathieu)
        {
            std::istringstream iss(line);
            std::string flag;
            iss >> flag;
            std::string charac;
            iss >> charac;
            // Look for the name of the voxel
            if(flag=="_Name")
            {
                voxelName = charac;
                foundName = true;
            }
            // Look for the flag "_Number"
            // And the characteristic "voxelBasePair"
            if(flag=="_Number" && charac=="voxelBasePair")
            {
                int numOfBp;
                iss >> numOfBp;
                if(!foundName)
                {
                    std::cerr<<"*** Fatal Error ***"<<std::endl;
                    std::cerr<<"ScanDamage::RetrieveVoxelBp: The number of bp was found before the name "
                    <<"of the voxel... This is an unexpected case."<<std::endl;
                    std::cerr<<"******"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                else
                {
                    fBpPerVoxel[voxelName] = numOfBp;
                    std::cout<<voxelName<<" has "<<numOfBp<<" bp"<<std::endl;
                    foundNumOfBp = true;
                }
            }

            if (flag == "_pl") // Mathieu
            {
                std::string mat;
                int strand;
                int copyNumber;

                double x,y,z;

                std::vector<double> coord;
                iss >> mat;
                iss >> strand;
                iss >> copyNumber;

                iss >> x;
                iss >> y;
                iss >> z;

                coord.push_back(x); // Recuperation des coordonnes dans l'objet coord
                coord.push_back(y);
                coord.push_back(z);
                //std::cout << "reading coordinate " << x << " " << y << " " << z << " " << copyNumber << " " << voxelName << std::endl;
                //getchar();

                if (charac == "phosphate1"){

                    ph1_coord.insert({copyNumber,coord});
                    
                    }

                if (charac == "phosphate2"){

                    ph2_coord.insert({copyNumber,coord});
                    
                    }
                
                if (charac == "deoxyribose1"){

                    deox1_coord.insert({copyNumber,coord});
                    
                    }

                if (charac == "deoxyribose2"){

                    deox2_coord.insert({copyNumber,coord});
                    
                    } // Ajout des coordonees dans les objets ph1_coord, ect... (Mathieu)
            }
        }
        file.close();

        f_ph1_coord.insert({voxelName,ph1_coord}); // Association des noms des voxels avec les coordonnes (Mathieu)
        f_ph2_coord.insert({voxelName,ph2_coord});
        f_deox1_coord.insert({voxelName,deox1_coord});
        f_deox2_coord.insert({voxelName,deox2_coord});
    }
    
    if (fBpPerVoxel.size() == 0) {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"ScanDamage::RetrieveVoxelBp: No Bp found in voxel definition files. \n Or make"
        <<" sure that file imp.info exists in working directory!"<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::FillVoxelData()
{   
    std::ifstream file(fCellDefFilePath);
    if(!file.good() )
    {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"FillVoxelData: No file named "<<fCellDefFilePath<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
    ullint bpCount = 0;
    unsigned int voxelCount = 0;
    int chromo_previous = 0;
    // Read the file line by line
    std::string line;
    while(std::getline(file, line) )
    {
        std::istringstream iss(line);
        std::string flag;
        iss >> flag;

        // Read cell geometry (Mathieu)
        if (flag == "_Type")
        {
            std::string type;
            iss >> type;
            if (type == "Ellipsoid")
            {
                iss >> fXCell;
                iss >> fYCell;
                iss >> fZCell;
                //std::cout << "Checking cell dimension " << fXCell << " " << fYCell << " " << fZCell << std::endl;
            }
        }
        // If the flag correspond to the placement of a voxel (Mathieu)
        if(flag == "_pl")
        {
            double x,y,z;
            double xx,xy,xz,yx,yy,yz,zx,zy,zz;

            std::string voxelName;
            iss >> voxelName;
            int chromo;
            iss >> chromo;
            int domain;
            iss >> domain;

            //--->Read coordinates a
            iss >> x; //in nm
            iss >> y; //in nm
            iss >> z; //in nm

            //--->Read rot. matrix
            iss >> xx;
            iss >> yx;
            iss >> zx;
            iss >> xy;
            iss >> yy;
            iss >> zy;
            iss >> xz;
            iss >> yz;
            iss >> zz;

            //std::cout << "Checking reading of voxel data ... " << std::endl;
            //std::cout << x << " " << y << " " << z << std::endl;
            //std::cout << xx << " " << xy << " " << xz << std::endl;
            //std::cout << yx << " " << yy << " " << yz << std::endl;
            //std::cout << zx << " " << zy << " " << zz << std::endl;

            //getchar();
            
            // If we change of chromosome then reset the number of bp.
            // Each chromosome starts at 0 bp.
            if(chromo != chromo_previous)
            {
                bpCount = 0;
                chromo_previous = chromo;
            }
            // Fill the data structure
            fVoxels.push_back( VoxelData(voxelName, chromo, domain, bpCount, x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz) ); // ajout dans fVoxels le voxel ainsi que ces coordonnées 
            int numBpinthisvoxel = int(fBpPerVoxel[voxelName]);
            bpCount += numBpinthisvoxel;// bpcount à vérifier mais bp dans chromosome 
            if (fChromosomeBpMap.find(chromo) == fChromosomeBpMap.end()) {
                fChromosomeBpMap.insert({chromo,numBpinthisvoxel});
            } else {
                fChromosomeBpMap[chromo] += numBpinthisvoxel;
            }
            voxelCount++;
        }
    }
    file.close();
    

    if (fVoxels.size() == 0) {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"ScanDamage::FillVoxelData: NofVoxels info found in files "<<fCellDefFilePath<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout<<"Num of voxels: "<< fVoxels.size()<<" placed in Cell Nucleus."<<std::endl;
    std::cout<<"=====Choromosome sizes====="<<std::endl;
    std::cout<<"Chromosome ID\t Number of Bp"<<std::endl;
    for (auto const& [chrom, nBp] : fChromosomeBpMap) {
         std::cout<<chrom<< "\t"<<nBp<<std::endl;
    }
    std::cout<<"==========================="<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ScanDamageFromPhys()
{
    std::cout<<"===== Start Scanning Damages From Phys =====\n";
    fEdepSumInNucleus = 0;
    fphysTables.clear();
    fphysSlectedTables.clear();
    fs::path currentP{fPathToOutputs+"/phys_output"};

    std::cout << fPathToOutputs+"/phys_output" << std::endl;

    fs::file_status s = fs::file_status{};
    auto isExist = fs::status_known(s) ? fs::exists(s) : fs::exists(currentP);

    std::cout << isExist << std::endl;
    //getchar();

    if (isExist) {
        bool isFoundRootFiles = false;
        for (const auto entry : fs::directory_iterator(currentP)) {
            if (entry.path().extension() == ".root") {
                std::cout <<"ScanDamageFromPhys(): Processing file: "<< entry.path().filename()<< std::endl;
                AnaPhysRootFile(entry.path());
                if (!isFoundRootFiles) isFoundRootFiles=true;
            }
        }
        if (!isFoundRootFiles) {
            std::cout<<"=====>> No root files found in folder \"phys_output\"!!! Skip Scanning Damages From Phys =====\n";
        }
        if (fphysTables.size() > 0) {
            SortPhysTableWithSelection();
        }
    } else {
        std::cout<<"=====>> Cannot find folder \"phys_output\"!!! Skip Scanning Damages From Phys =====\n";
    }

    std::cout<<"===== End Scanning Damages From Phys =====\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ScanDamageFromChem()
{
    std::cout<<"===== Start Scanning Damages From Chem =====\n";
    fs::path currentP{fPathToOutputs+"/chem_output"};
    fs::file_status s = fs::file_status{};
    auto isExist = fs::status_known(s) ? fs::exists(s) : fs::exists(currentP);
    if (isExist) {
        bool isFoundRootFiles = false;
        for (const auto entry : fs::directory_iterator(currentP)) {
            if (entry.path().extension() == ".root") {
                AnaChemRootFile(entry);
                if (!isFoundRootFiles) isFoundRootFiles=true;
            }
        }
        if (!isFoundRootFiles) {
            std::cout<<"=====>> No root files found in folder \"chem_ouput\"!!! Skip Scanning Damages From Chem =====\n";
            fSkipScanningIndirectDamage = true;
        }
        if (fchemTables.size() > 0) {
            SortChemTableWithSelection();
        }
    } else {
        std::cout<<"=====>> Cannot find folder \"chem_ouput\"!!! Skip Scanning Damages From Chem =====\n";
        fSkipScanningIndirectDamage = true;
    }
    
    std::cout<<"===== End Scanning Damages From Chem =====\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootFile(const std::string fileName)
{
    TFile* f = new TFile(fileName.c_str());
    if(f->IsZombie() ){
        // File is corrupted
        std::cerr<<"*********** Warning *************"<<std::endl;
        std::cerr<<"The file "<<fileName<<" seems to be corrupted..."<<std::endl;
        std::cerr<<"We will skip it."<<std::endl;
        std::cerr<<"**********************************"<<std::endl;
        return;
    }
    AnaPhysRootTree1(f);
    AnaPhysRootTree2(f);
    f->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaChemRootFile(fs::directory_entry entry)
{
    /*std::string fnameWithoutExtension = entry.path().stem().string();
    auto [eventNumber, voxelNumber] = GetEventNberAndVoxelNberFromChemRoot(fnameWithoutExtension);
    
    // Voxel data retrieval
    VoxelData& voxelData = fVoxels.at(voxelNumber);

    int chromo = voxelData.fChromosome;
    ullint firstBpNum = voxelData.fFirstBpCopyNum; // Mathieu
    */
    // *******************
    // Analyse of the ntuple to detect SB and associate them with a bpNumCorrected
    // *******************

    // Load the file, the directory and the ntuple
    //TFile f(entry.path().c_str());
    TFile* f = new TFile(entry.path().c_str()); // Mathieu
    if(f->IsZombie() ){
        corruptedFiles++;
        // File is corrupted
        std::cerr<<"*********** Warning *************"<<std::endl;
        std::cerr<<"The file "<<entry.path().string()<<" seems to be corrupted..."<<std::endl;
        std::cerr<<"We will skip it."<<std::endl;
        std::cerr<<"Number of corrupted files: "<< corruptedFiles<<std::endl;
        std::cerr<<"**********************************"<<std::endl;
    } else {
        //TDirectoryFile *d = dynamic_cast<TDirectoryFile*> (f.Get("ntuple") ); //(Mathieu)
        TTree* chemTree = (TTree*) f->Get("ntuple_2"); //Mathieu remplacement de d par f (file)

        if( (int) chemTree->GetEntries() >0)
        {
            double strand;
            double copyNumber;
            //double xp; //Mathieu
            //double yp;
            //double zp;
            double time;
            double base;
            double dnaMol;
            double eventNumber;
            double voxelNumber;
            chemTree->SetBranchAddress("strand", &strand);
            chemTree->SetBranchAddress("copyNumber", &copyNumber);
            //chemTree->SetBranchAddress("xp", &xp); // zp, yp et zp sont les coordonnées de la trace microdosimétrique (Mathieu)
            //chemTree->SetBranchAddress("yp", &yp);
            //chemTree->SetBranchAddress("zp", &zp);
            chemTree->SetBranchAddress("time", &time);
            chemTree->SetBranchAddress("base", &dnaMol);

            chemTree->SetBranchAddress("eventNumber", &eventNumber);
            chemTree->SetBranchAddress("voxelCopyNumber", &voxelNumber);

            unsigned int entryNumber = (int) chemTree->GetEntries();
            for (unsigned int e=0;e<entryNumber;e++)
            {
                chemTree->GetEntry(e);

                VoxelData& voxelData = fVoxels.at(size_t(voxelNumber));

                int chromo = voxelData.fChromosome;
                int domain = voxelData.fDomain;
                ullint firstBpNum = voxelData.fFirstBpCopyNum;

                ullint cpNumCorrected = firstBpNum+int(copyNumber);
                if (dnaMol == 8 || dnaMol == 9){base = 0;}
                else {base = 1;}


                double xDna, yDna, zDna;

                if (int(strand) == 1 and int(dnaMol) == 8) {
                xDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][0]; // xDna, yDna et zDna sont tous les trois les coordonnées des molécules touchés dans le voxel (Mathieu)
                yDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][1];
                zDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][2];

                std::cout << "Coordinates deoxyribose1" << "xDna = " << xDna << " yDna = " << yDna << " zDna = " << zDna << std::endl;
                }

                if (int(strand) == 1 and int(dnaMol) == 9) {
                xDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][0];
                yDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][1];
                zDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][2];

                std::cout << "Coordinates phosphate1" << "xDna = " << xDna << " yDna = " << yDna << " zDna = " << zDna << std::endl;

                }

                if (int(strand) == 2 and int(dnaMol) == 8) {
                xDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][0];
                yDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][1];
                zDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][2];

                std::cout << "Coordinates deoxyribose2" << "xDna = " << xDna << " yDna = " << yDna << " zDna = " << zDna << std::endl;
                }

                if (int(strand) == 2 and int(dnaMol) == 9) {
                xDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][0];
                yDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][1];
                zDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][2];

                std::cout << "Coordinates phosphate2" << "xDna = " << xDna << " yDna = " << yDna << " zDna = " << zDna << std::endl;
                }

                double x1,y1,z1;

                x1 = voxelData.fxx * xDna + voxelData.fxy * yDna + voxelData.fxz * zDna; // Apply rotation (Mathieu)
                y1 = voxelData.fyx * xDna + voxelData.fyy * yDna + voxelData.fyz * zDna; 
                z1 = voxelData.fzx * xDna + voxelData.fzy * yDna + voxelData.fzz * zDna;

                x1 += fXCell + voxelData.fx; // Apply translation (Mathieu)
                y1 += fYCell + voxelData.fy; // Application de la translation (Mathieu)
                z1 += fZCell + voxelData.fz;
                




                            

                std::vector<ullint> newLine{ 
                    (ullint)eventNumber,
                    (ullint)strand,
                    cpNumCorrected,
                    (ullint)base,
                    (ullint)(time*1000000),
                    (ullint)domain,
                    (ullint)(x1*1000),
                    (ullint)(y1*1000),
                    (ullint)(z1*1000)}; //Mathieu *1000 pour garder la precision des chiffres apres la virgule

                auto itr = fchemTables.find(chromo);
                if ( itr == fchemTables.end()) {
                    Table tbforThisChro{newLine};
                    fchemTables.insert({chromo,tbforThisChro});
                } else (itr->second).push_back(newLine);
            }
        }
    }//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::tuple<unsigned int, unsigned int> ScanDamage::GetEventNberAndVoxelNberFromChemRoot(
    const std::string fileNameWithoutExtension)
{
    unsigned int evnN, volxelN;
    auto fristPos = fileNameWithoutExtension.find_first_of("_");
    auto secondPos = fileNameWithoutExtension.substr(fristPos+1).find_first_of("_");
    auto lastPos = fileNameWithoutExtension.find_last_of("_");
    evnN = std::stoul(fileNameWithoutExtension.substr(fristPos+1,secondPos));
    volxelN = std::stoul(fileNameWithoutExtension.substr(lastPos+1));
    return  {evnN, volxelN};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootTree1(TFile* f)
{
    TDirectoryFile* d = dynamic_cast<TDirectoryFile*> (f->Get("ntuple") );
    TTree* tPhys = dynamic_cast<TTree*> (d->Get("ntuple_1") );

    if( tPhys->GetEntries()  > 0)
        {
        double flagParticle;
        //double flagParentID; //Mathieu
        double flagProcess;
        //double x;
        //double y;
        //double z;
        double edep;
        double eventNumber;
        double volumeName;
        double copyNumber;
        double lastMetVoxelCopyNum;

        tPhys->SetBranchAddress("flagParticle", &flagParticle);
        //tPhys->SetBranchAddress("flagParentID", &flagParentID); // Mathieu
        tPhys->SetBranchAddress("flagProcess", &flagProcess);
        //tPhys->SetBranchAddress("x", &x);
        //tPhys->SetBranchAddress("y", &y);
        //tPhys->SetBranchAddress("z", &z);
        tPhys->SetBranchAddress("edep", &edep);
        tPhys->SetBranchAddress("eventNumber", &eventNumber);
        tPhys->SetBranchAddress("volumeName", &volumeName);
        tPhys->SetBranchAddress("copyNumber", &copyNumber);
        tPhys->SetBranchAddress("lastMetVoxelCopyNum", &lastMetVoxelCopyNum);

        unsigned int entryNumber =  tPhys->GetEntries() ;

        // Loop on all the "lines" (ie entry) of the ntuple
        for(unsigned int e=0; e<entryNumber; e++)
        {
            // Set all the variables to the values corresponding to the entry number
            tPhys->GetEntry(e);
            // Check if the process is an ionisation
            // Only ionisation should trigger the removal of a DNA molecule from the chemical step
            if(flagProcess == 13 // e-_DNAIonisation
                    || flagProcess == 113 // e-_DNAPTBIonisation
                    || flagProcess == 18 // proton_DNAIonisation
                    || flagProcess == 21 // hydrogen_DNAIonisation
                    || flagProcess == 24 // alpha_DNAIonisation
                    || flagProcess == 27 // alpha+_DNAIonisation
                    || flagProcess == 31 // helium_DNAIonisation
                    || flagProcess == 12 // e-_DNAExcitation
                    || flagProcess == 112 // e-_DNAPTBExcitation
                    || flagProcess == 15 // e-_DNAVibExcitation
                    || flagProcess == 17 // proton_DNAExcitation
                    || flagProcess == 20 // hydrogen_DNAExcitation
                    || flagProcess == 23 // alpha_DNAExcitation
                    || flagProcess == 26 // alpha+_DNAExcitation
                    || flagProcess == 30 // helium_DNAExcitation
                    || flagProcess == 33 // GenericIon Ionization (Mathieu)
                    ) {
                // Check the interaction happened in a dna molecule or its hydration shell
                if(volumeName == 1 // d1
                        || volumeName == 11 // p1
                        || volumeName == 2 // d2
                        || volumeName == 22 // p2
                        || volumeName == 7 // d1_w
                        || volumeName == 71 // p1_w
                        || volumeName == 8 // d2_w
                        || volumeName == 81 // p2_w
                        )
                {
                    // *************

                    // Retrieve the voxel copy number
                    double voxelCopyNumber  = lastMetVoxelCopyNum;
                    if (voxelCopyNumber >= 0 && voxelCopyNumber<fVoxels.size()){
                        // Chromosome, domain and firstNucleotideNum
                        const VoxelData& voxelData = fVoxels.at(size_t(voxelCopyNumber) );
                        int chromo = voxelData.fChromosome;
                        int domain = voxelData.fDomain;
                        ullint firstBpCN = voxelData.fFirstBpCopyNum;
                        ullint cpNumCorrected = firstBpCN+int(copyNumber);

                        // Get the event number
                        double eventNum = eventNumber;

                        // Determine the strand
                        double strand (-1);
                        if(volumeName==1
                                || volumeName==11
                                || volumeName==7
                                || volumeName==71
                                || volumeName==6 // ade
                                || volumeName==9 // ade
                                || volumeName==4 // gua
                                || volumeName==10) // gua
                            strand = 1;
                        else if(volumeName==2
                                || volumeName==22
                                || volumeName==8
                                || volumeName==81
                                || volumeName==5 // thy
                                || volumeName==12 // thy
                                || volumeName==3 // cyto
                                || volumeName==13) // cyto
                            strand = 2;

                        // Retrieve DNA damage coordinate (Mathieu)

                        double xDna, yDna, zDna;
                        // ici, toutes ces conditions sont valable sur le backbone uniquement (Mathieu)
                        if (volumeName == 1 
                            || volumeName == 7)
                            {
                                if (copyNumber < f_deox1_coord[voxelData.fVoxelName].size()){
                                xDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][0];
                                yDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][1];
                                zDna = f_deox1_coord[voxelData.fVoxelName][int(copyNumber)][2];
                                }
                            }

                        if (volumeName == 2 
                            || volumeName == 8)
                            {
                                if (copyNumber < f_deox2_coord[voxelData.fVoxelName].size()){
                                xDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][0];
                                yDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][1];
                                zDna = f_deox2_coord[voxelData.fVoxelName][int(copyNumber)][2];
                                }
                            }

                        if (volumeName == 11 
                            || volumeName == 71)
                            {

                                if (copyNumber < f_ph1_coord[voxelData.fVoxelName].size()){
                                xDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][0];
                                yDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][1];
                                zDna = f_ph1_coord[voxelData.fVoxelName][int(copyNumber)][2];
                                }
                            }

                        if (volumeName == 22 
                            || volumeName == 81)
                            {
                                if (copyNumber < f_ph2_coord[voxelData.fVoxelName].size()){
                            xDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][0];
                            yDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][1];
                            zDna = f_ph2_coord[voxelData.fVoxelName][int(copyNumber)][2];
                                }
                            }

                        double x1,y1,z1;

                        x1 = voxelData.fxx * xDna + voxelData.fxy * yDna + voxelData.fxz * zDna; // Apply rotation (Mathieu)
                        y1 = voxelData.fyx * xDna + voxelData.fyy * yDna + voxelData.fyz * zDna; // Application de la rotation inverse (Mathieu)
                        z1 = voxelData.fzx * xDna + voxelData.fzy * yDna + voxelData.fzz * zDna;

                        x1 += fXCell + voxelData.fx; // Apply translation (Mathieu)
                        y1 += fYCell + voxelData.fy; // Application de la translation (Mathieu)
                        z1 += fZCell + voxelData.fz;

                        

                        std::cout << "Cheking coordinate " << "x1 = " << x1 << " y1 = " << y1 << " z1 = " << z1 << std::endl;



                        // Check if the chromo has already been registered
                        std::vector<ullint> newLine{ 
                            (ullint)eventNum,
                            (ullint)strand,
                            cpNumCorrected,
                            (ullint)volumeName,
                            (ullint)flagProcess,
                            (ullint)(edep*1000000),
                            (ullint)domain,
                            (ullint)(x1*1000),
                            (ullint)(y1*1000),
                            (ullint)(z1*1000)}; // Mathieu
                            // *1000000 in edep is taken back after. This is done 
                            // because the number is "unsigned long long int" instead of "double".

                        auto itr = fphysTables.find(chromo);
                        if ( itr == fphysTables.end()) {
                            Table tbforThisChro{newLine};
                            fphysTables.insert({chromo,tbforThisChro});
                        } else (itr->second).push_back(newLine);
                    }
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootTree2(TFile* f)
{
    TDirectoryFile* d2 = dynamic_cast<TDirectoryFile*> (f->Get("ntuple") );
    TTree* tPhys2 = dynamic_cast<TTree*> (d2->Get("ntuple_3") );

    if( int(tPhys2->GetEntries() ) > 0)
    {
        double edep;
        double eventNumber;

        tPhys2->SetBranchAddress("edep", &edep);
        tPhys2->SetBranchAddress("eventNumber", &eventNumber);

        unsigned int entryNumber = int( tPhys2->GetEntries() );

        // Loop on all the "lines" (ie entry) of the ntuple
        for(unsigned int e=0; e<entryNumber; e++)
        {
            tPhys2->GetEntry(e);
            fEdepSumInNucleus += edep;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::SortPhysTableWithSelection()
{
    for(const auto [chrom, physTable] : fphysTables)
    {
        // Final table
        Table physTableWithSelection;

        std::map<ullint,std::map<ullint,std::map<ullint, ullint > > > energyMap;

        std::map<ullint,std::map<ullint,std::map<ullint, std::vector <ullint> > > > infoMap;
        
        // Loop on all the lines of the table
        for(unsigned int line=0, eline=physTable.size(); line<eline; ++line)
        {
            ullint eventNum = physTable[line][0];
            ullint strand = physTable[line][1];
            ullint copyNumber = physTable[line][2];
            ullint volumeFlag = physTable[line][3];
            ullint processFlagr = physTable[line][4];
            ullint energy = physTable[line][5];

            ullint dom = physTable[line][6];
            ullint x = physTable[line][7];
            ullint y = physTable[line][8];
            ullint z = physTable[line][9];

            // Cumulate the energy value
            energyMap[eventNum][strand][copyNumber] += energy;

            std::vector<ullint> info{dom,x,y,z};

            infoMap[eventNum][strand][copyNumber] = info;

        }

        int notDuplicatedLine = 0;

        // Loop on all the events
        std::map<ullint, std::map<ullint, std::map<ullint, ullint> > >::iterator iit = energyMap.begin();
        std::map<ullint, std::map<ullint, std::map<ullint, ullint> > >::iterator iite = energyMap.end();

        std::map<ullint, std::map<ullint, std::map<ullint, std::vector <ullint> > > >::iterator iit2 = infoMap.begin();

        for(; iit!=iite;++iit,++iit2)
        {
            ullint eventNum = iit->first;

            // Loop on all the strands
            std::map<ullint, std::map<ullint, ullint> >::iterator itt = iit->second.begin();
            std::map<ullint, std::map<ullint, ullint> >::iterator itte = iit->second.end();

            std::map<ullint, std::map<ullint, std::vector<ullint> > >::iterator itt2 = iit2->second.begin();

            for(; itt!=itte;++itt,++itt2)
            {
                ullint strand = itt->first;
                // Loop on all the copy numbers
                std::map<ullint, ullint>::iterator ittt = itt->second.begin();
                std::map<ullint, ullint>::iterator ittte = itt->second.end();

                std::map<ullint, std::vector <ullint> >::iterator ittt2 = itt2->second.begin();

                for(; ittt!=ittte;++ittt,++ittt2)
                {
                    ullint copyNumber = ittt->first;
                    double currentE = double(ittt->second) / 1000000; // eV

                    std::vector<ullint> info = ittt2->second;

                    // Energy condition(s) are set here
                    bool fill = false;                  
                     // Threshold condition
                     if(currentE < fThresholdEnergy)
                          fill=false;
                     else
                          fill=true;

                    if(fill)
                    {
                        // Add a line
                        physTableWithSelection.push_back(std::vector<ullint>());
                        // Fill the line
                        physTableWithSelection[notDuplicatedLine].push_back(eventNum);
                        physTableWithSelection[notDuplicatedLine].push_back(strand);
                        physTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                        physTableWithSelection[notDuplicatedLine].push_back(0);
                        physTableWithSelection[notDuplicatedLine].push_back(0);
                        physTableWithSelection[notDuplicatedLine].push_back((ullint)(currentE*1000000) );
                        physTableWithSelection[notDuplicatedLine].push_back(0);

                        physTableWithSelection[notDuplicatedLine].push_back(info[0]);
                        physTableWithSelection[notDuplicatedLine].push_back(info[1]);
                        physTableWithSelection[notDuplicatedLine].push_back(info[2]);
                        physTableWithSelection[notDuplicatedLine].push_back(info[3]);

                        ++notDuplicatedLine;
                    }
                }
            }
        }
        // *******************************************
        // Print the "physTableWithSelection" table for the current chromosome
        // *******************************************
        std::cout << "### Phys SB for chromosome "<<chrom<<" : " << physTableWithSelection.size() << " ###" << std::endl;
        if (physTableWithSelection.size()>0) fphysSlectedTables.insert({chrom,physTableWithSelection});
    }
    fphysTables.clear();//Free memory
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::SortChemTableWithSelection()
{
    for (const auto& [chromo,chemTable] : fchemTables)
    {
      	// Final table
        Table chemTableWithSelection;
        
        int notDuplicatedLine = 0;

        for(unsigned int line=0, eline=chemTable.size(); line<eline; ++line)
        {
            ullint eventNum = chemTable[line][0];
            ullint strand = chemTable[line][1];
            ullint copyNumber = chemTable[line][2];
            ullint base = chemTable[line][3];
            ullint time = chemTable[line][4];

            ullint domain = chemTable[line][5];
            ullint x = chemTable[line][6];
            ullint y = chemTable[line][7];
            ullint z = chemTable[line][8];
        
     		// Random number between 0 and <1
			if (base==1) {
                // Add a line
                chemTableWithSelection.push_back(std::vector<ullint>());
                // Fill the line
                chemTableWithSelection[notDuplicatedLine].push_back(eventNum);
                chemTableWithSelection[notDuplicatedLine].push_back(strand);
                chemTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                chemTableWithSelection[notDuplicatedLine].push_back(base);
                chemTableWithSelection[notDuplicatedLine].push_back(time);
                chemTableWithSelection[notDuplicatedLine].push_back(0);
                chemTableWithSelection[notDuplicatedLine].push_back(1);  

                chemTableWithSelection[notDuplicatedLine].push_back(domain);  
                chemTableWithSelection[notDuplicatedLine].push_back(x);   
                chemTableWithSelection[notDuplicatedLine].push_back(y);    
                chemTableWithSelection[notDuplicatedLine].push_back(z);   
                ++notDuplicatedLine;
			}
			else { 
			    //double r = double(std::rand() ) / RAND_MAX;
                double r = gRandomGen.Rndm(); 
                //std::cout<<r<<std::endl;
			    if(r <= fProbabilityForIndirectSB){
                    // Add a line
                    chemTableWithSelection.push_back(std::vector<ullint>());
                    // Fill the line
                    chemTableWithSelection[notDuplicatedLine].push_back(eventNum);
                    chemTableWithSelection[notDuplicatedLine].push_back(strand);
                    chemTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                    chemTableWithSelection[notDuplicatedLine].push_back(base);
                    chemTableWithSelection[notDuplicatedLine].push_back(time);
                    chemTableWithSelection[notDuplicatedLine].push_back(0);
                    chemTableWithSelection[notDuplicatedLine].push_back(1); 

                    chemTableWithSelection[notDuplicatedLine].push_back(domain);  
                    chemTableWithSelection[notDuplicatedLine].push_back(x);   
                    chemTableWithSelection[notDuplicatedLine].push_back(y);    
                    chemTableWithSelection[notDuplicatedLine].push_back(z);                           
                    ++notDuplicatedLine;
                }
	        }
      	}   	
        if (chemTableWithSelection.size() > 0) fchemSlectedTables.insert({chromo,chemTableWithSelection});
	}
    fchemTables.clear();//free memory
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::MergeDamageFromPhysChem()
{
    for(int  chromo=0; chromo<46; chromo++)
    {
        // *****************************
        // Add one table after the other
        // *****************************

        // MergedTable to be built
        Table mergedTable;

        auto itr = fphysSlectedTables.find(chromo);
        if(itr != fphysSlectedTables.end())
        {
            Table physTable= itr->second;
            for(auto const& vec : physTable) mergedTable.push_back(vec);
        }
        itr = fchemSlectedTables.find(chromo);
        if(itr != fchemSlectedTables.end())
        {
            Table chemTable = itr->second;
            for(auto const& vec : chemTable) mergedTable.push_back(vec);
        }

        // *******************************************
        // Sort the merged table to put the event in the correct order
        // *******************************************
        Trier tri;
        Table::iterator it = mergedTable.begin();
        Table::iterator ite = mergedTable.end();
        std::sort(it, ite, tri);

        // *******************************************
        // Delete duplicate SB
        // *******************************************

        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>> removeDuplicateValueMap;

        // Put all the values of the mergedTable in the map created just above.
        // Loop on all the lines of the mergedTable
        for(unsigned int line=0; line<mergedTable.size(); ++line)
        {
            ullint eventNum = mergedTable[line][0];
            ullint strand = mergedTable[line][1];
            ullint copyNumber = mergedTable[line][2];
            ullint base = mergedTable[line][3];
            std::vector<ullint> lineV;
            // If more elements are presents, add them here
            int lineSize = mergedTable[line].size();
            if(lineSize > 4)
            {
                for(int i=4; i<lineSize; i++)
                {
                    lineV.push_back(mergedTable[line][i]);
                }
                removeDuplicateValueMap[eventNum][strand][copyNumber][base]=lineV;
            }
            else
            {
                lineV.push_back(0);
                removeDuplicateValueMap[eventNum][strand][copyNumber][base]=lineV;
            }
        }

        // *******************************************
        // Create the "mergedTableWithoutDuplicatedSB" table
        // *******************************************

        // At this point, we have created a map named "removeDuplicateValueMap" that organized all the mergedTable content
        // AND that does not contain any duplicate because of the override caracteristic of a map.
        // Indeed, doing map[2] = "hello" followed by map[2]  = "bye" will put map[2] value to "bye" because it overrided the first "hello".
        // This is a cheap way to remove duplicates by overriding them.

        // The next part is dedicated to the creation of the "mergedTableWithoutDuplicatedSB" table from the "removeDuplicateValueMap" map.
        // Only the event, strand and copynumber and  will be put in this final table.
        Table mergedTableWithoutDuplicatedSB;
        Table mergedTableWithoutDuplicatedSBandbases;
        int notDuplicatedLine = 0;
        int notDuplicatedLine2= 0;
        // Loop on all the events
        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>>::iterator 
            iit = removeDuplicateValueMap.begin();
        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>>::iterator 
            iite = removeDuplicateValueMap.end();
        for(; iit!=iite;++iit)
        {
            ullint eventNum = iit->first;
            // Loop on all the strands
            std::map<ullint, std::map<ullint, std::map<ullint, std::vector<ullint > > > >::iterator 
                itt = iit->second.begin();
            std::map<ullint, std::map<ullint, std::map<ullint, std::vector<ullint > > > >::iterator 
                itte = iit->second.end();
            for(; itt!=itte;++itt)
            {
                ullint strand = itt->first;
                // Loop on all the copy numbers
                std::map<ullint, std::map<ullint, std::vector<ullint > > >::iterator ittt = itt->second.begin();
                std::map<ullint, std::map<ullint, std::vector<ullint > > >::iterator ittte = itt->second.end();
                for(; ittt!=ittte;++ittt)
                {
                    ullint copyNumber = ittt->first;                    
		            // Loop on all the base flag
                    std::map<ullint, std::vector<ullint > >::iterator itttt = ittt->second.begin();
                    std::map<ullint, std::vector<ullint > >::iterator itttte = ittt->second.end();
                    for(; itttt!=itttte;++itttt)
                    {
                        ullint base = itttt->first;
                        // Fill the table
                        if (strand>0 && strand<3)
                        {

                            ullint time = removeDuplicateValueMap[eventNum][strand][copyNumber][base][0];
                            ullint energy = removeDuplicateValueMap[eventNum][strand][copyNumber][base][1];
                            ullint type = removeDuplicateValueMap[eventNum][strand][copyNumber][base][2]; // type physique = 0 et chem = 1 ou vis vers ca 
                            ullint domain = removeDuplicateValueMap[eventNum][strand][copyNumber][base][3];
                            ullint x = removeDuplicateValueMap[eventNum][strand][copyNumber][base][4];
                            ullint y = removeDuplicateValueMap[eventNum][strand][copyNumber][base][5];
                            ullint z = removeDuplicateValueMap[eventNum][strand][copyNumber][base][6];

                            // Add a line
                            mergedTableWithoutDuplicatedSB.push_back(std::vector<ullint>());                           
                            // Fill the line
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(eventNum);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(strand);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(copyNumber);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(base);

                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(time);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(energy);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(type);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(domain);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(x);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(y);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(z);
                            ++notDuplicatedLine;                           
                            if (base==0)
                            {
                                // Add a line
                                mergedTableWithoutDuplicatedSBandbases.push_back(std::vector<ullint>());                              
                                // Fill the line
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(eventNum);
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(strand);
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(copyNumber);                               
                                ++notDuplicatedLine2;
                            }
                        }
                    }
                }
            }
        }

        // *******************************************
        // Print the "mergedTableWithoutDuplicatedSB" table for the current chromosome
        // *******************************************
        if (mergedTableWithoutDuplicatedSB.size() > 0) {
            //PrintTable(fMergeFolder + "/chromo_"+std::to_string(chromo)+".dat", 
            //mergedTableWithoutDuplicatedSB, "eventNum, strand, copyNumber, isbase, 
            //time(ns*1000000), edep, phy:0 chem:1");
            fMergedTables.insert({chromo,mergedTableWithoutDuplicatedSB});
        }
        mergedTable.clear();
        mergedTableWithoutDuplicatedSBandbases.clear();
        mergedTableWithoutDuplicatedSB.clear();
    }
    //free memory:
    fphysSlectedTables.clear();
    fchemSlectedTables.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ReadCellandVoxelDefFilePaths()
{
    //fs::path thisP = fs::current_path();
    fs::path thisP = fPathToOutputs;//Mathieu
    for (const auto entry : fs::directory_iterator(thisP)){
        if (entry.path().filename() == "imp.info") {
            std::ifstream file(entry.path().c_str());
            if(!file.good() ){
                std::cerr<<"**** Fatal Error *****"<<std::endl;
                std::cerr<<"ScanDamage::ReadCellandVoxelDefFilePaths(): File corupted: "
                <<entry.path()<<std::endl;
                std::cerr<<"*************** *****"<<std::endl;
                exit(EXIT_FAILURE);
            }

            std::string line;
            while(std::getline(file, line) ){
                std::istringstream iss(line);
                std::string flag;
                iss >> flag;
                if ( flag == "_geovolxelpath") {
                    std::string voxname;
                    iss >> voxname;
                    fVoxelDefFilesList.insert(voxname);
                }
                if ( flag == "_geocellpath") {
                    std::string cellpname;
                    iss >> cellpname;
                    fCellDefFilePath = cellpname;
                }
                if ( flag == "_numberOfBasepairs") {
                    iss >> fTotalNbBpPlacedInGeo;
                }
                if ( flag == "_numberOfHistones") {
                    iss >> fTotalNbHistonePlacedInGeo;
                }
                if ( flag == "_nucleusVolume") {
                    iss >> fNucleusVolume;
                }
                if ( flag == "_nucleusMassDensity") {
                    iss >> fNucleusMassDensity;
                }
                if ( flag == "_nucleusMass") {
                    iss >> fNucleusMass;
                }
            }
            file.close();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....