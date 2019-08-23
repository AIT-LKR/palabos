#include <../src/core/globalDefs.h>
#include <../src/io/parallelIO.h>
#include <../externalLibraries/tinyxml/tinyxml.h>

#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>

#include "./tools.h"
#include "./tools.hh"
#include "./tools_files.h"
#include "./tools_files.hh"
#include "./plb_physicalFlowParam.h"
#include "./tools_plb.h"

typedef double T;

#ifndef TOOLS_PLB_HH
#define TOOLS_PLB_HH

using namespace std;
using namespace plb;


void renameAllExisting(string fileBaseName) {
    int existingFiles=0;
    while (doesFileExist(fileBaseName+(static_cast<ostringstream*>( &(ostringstream() << existingFiles) )->str()))) existingFiles++;
    while (existingFiles>0) {
        string sourceName = fileBaseName+(static_cast<ostringstream*>( &(ostringstream() << (existingFiles-1)) )->str());
        string destinationName = fileBaseName+(static_cast<ostringstream*>( &(ostringstream() << existingFiles) )->str());
        copyFile(sourceName, destinationName);
        existingFiles--;
    }
}


bool copyFile(string sourceName, string destinationName) {
    bool success;
    if (global::mpi().isMainProcessor()) {
        std::ifstream src(sourceName.c_str(), std::ios::binary);
        std::ofstream dst(destinationName.c_str(), std::ios::binary);
        if (dst << src.rdbuf())
            success = true;
        else
            success = false;
    }
    global::mpi().bCast(&success, 1);
    plbIOError(!success, std::string("Unsuccessful writing into file.")+destinationName);
    return success;
}

void rewriteFile(std::string data, std::string fileName) {
    bool success;
    if (global::mpi().isMainProcessor()) {
        plb::plb_ofstream ofile;
        ofile.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc );
        ofile << data;
        ofile.close();
        success = true;
    }
    global::mpi().bCast(&success, 1);
}

void appendToFile(std::string data, std::string fileName) {
    plb::plb_ofstream ofile;
    ofile.open(fileName.c_str(), std::ios_base::out | std::ios_base::app );
    ofile << data;
    ofile.close();
}


/* **************************************************************************** */

template<typename T>
PoiseuilleVelocity<T>::PoiseuilleVelocity(T uLB_, T uDev_, plint powerPoiseuilleVel_, char dir_,
        plint inletCentre_, plint inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
     inletCentreA(inletCentre_), inletCentreB(-1), inletRadius(inletRadius_),
     seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T>(0., uDev);
}

template<typename T>
PoiseuilleVelocity<T>::PoiseuilleVelocity(T uLB_, T uDev_, plint powerPoiseuilleVel_, char dir_,
        plint inletCentreA_, plint inletCentreB_, plint inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
      inletCentreA(inletCentreA_), inletCentreB(inletCentreB_), inletRadius(inletRadius_),
      seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T>(0., uDev);
}

template<typename T>
void PoiseuilleVelocity<T>::operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {
    T main_velocity = uLB;
    T side_velocity = 0.;
    T down_velocity = 0.;
    if (uDev) {
        main_velocity = (*main_distribution)(*generator);
        side_velocity = (*side_distribution)(*generator)/4;
        down_velocity = (*side_distribution)(*generator)/4;
    }
    plint r;

    if (dir == 'x') {
        if (inletCentreB == -1)
            r = distanceFromPoint(iY,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint(iY,iZ, inletCentreA,inletCentreB);
        u[0] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[1] = side_velocity;
        u[2] = down_velocity;
    }
    if (dir == 'y') {
        if (inletCentreB == -1)
            r = distanceFromPoint(iX,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint(iX,iZ, inletCentreA,inletCentreB);
        u[0] = side_velocity;
        u[1] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[2] = down_velocity;
    }
}

template<typename T>
void PoiseuilleVelocity<T>::operator()(plint iX, plint iY, Array<T,2>& u) const {
    T main_velocity = uLB;
    T side_velocity = 0.;
    if (uDev) {
        main_velocity = (*main_distribution)(*generator);
        side_velocity = (*side_distribution)(*generator)/2;
    }

    if( dir == 'x' ) {
        u[0] = poiseuilleVelocity(iY-inletCentreA, main_velocity, powerPoiseuilleVel, inletRadius);
        u[1] = side_velocity;
    }
    if( dir == 'y' ) {
        u[0] = side_velocity;
        u[1] = poiseuilleVelocity(iX-inletCentreA, main_velocity, powerPoiseuilleVel, inletRadius);
    }
}

template<typename T>
void writeLogFile(PhysicalFlowParam<T> const param,
                  std::string const& fileName) {
    stringstream logText;
    logText << endl
    << "--------------- PhysicalFlowParam: --------------"                  << "\n"
    << "-- dimensionless:"                                                  << "\n"
    << " Reynolds number               Re = " << param.getRe()              << "\n"
    << " System length                 lx = " << param.getlx()              << "\n"
    << "        width                  ly = " << param.getly()              << "\n"
    << "        height                 lz = " << param.getlz()              << "\n"
    << "-- physical:"                                                       << "\n"
    << " Velocity:      [m/s]        phU0 = " << param.getPhysicalU()       << "\n"
    << " Char. length   [m]          phL0 = " << param.getPhysicalLength()  << "\n"
    << " Char. time     [s]          phT0 = " << param.getPhysicalT()       << "\n"
    << " Viscosity      [m2/s]       phNu = " << param.getPhysicalNu()      << "\n"
    << " System length  [m]          phLx = " << param.getLx()              << "\n"
    << "        width   [m]          phLy = " << param.getLy()              << "\n"
    << "        height  [m]          phLz = " << param.getLz()              << "\n"
    << "-- numerical:"                                                      << "\n"
    << " Velocity                    lbU0 = " << param.getLatticeU()        << "\n"
    << " Grid spacing deltaX       deltaX = " << param.getDeltaX()          << "\n"
    << " Time step    deltaT       deltaT = " << param.getDeltaT()          << "\n"
    << " Viscosity                   lbNu = " << param.getLatticeNu()       << "\n"
    << " Lattice resolution"                                                << "\n"
    << "        nodes per meter         r = " << param.getNodesPerMeter()   << "\n"
    << "        nodes per char. length  x = " << param.getResolution()      << "\n"
    << " Extent of the system:         Nx = " << param.getNx()              << "\n"
    << "                               Ny = " << param.getNy()              << "\n"
    << "                               Nz = " << param.getNz()              << "\n"
    << " Relaxation frequency       omega = " << param.getOmega()           << "\n"
    << " Relaxation time              tau = " << param.getTau()             << "\n"
    << "-------------------------------------------------"                  << "\n";
    // print to screen
    pcout << logText.str();
    // save as file
    rewriteFile(logText.str(), fileName);
}


namespace plb {

template<typename T>
TopoWriter<T>::TopoWriter(plint valueRange_)
    : valueRange(valueRange_)
{ }

template<typename T>
void TopoWriter<T>::writePpm (
        std::string const& fName, bool colored,
        MultiScalarField2D<T>& field) const
{
    global::profiler().start("io");
    ScalarField2D<T> localField(field.getNx(), field.getNy());
    copySerializedBlock(field, localField);
    writePpmImplementation(fName, colored, localField);
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writePpm (
        std::string const& fName, bool colored,
        MultiScalarField3D<T>& field, int dir) const
{
    global::profiler().start("io");
    plint nx=field.getNx();
    plint ny=field.getNy();
    plint nz=field.getNz();

    if (dir == 0) {
        ScalarField2D<T> localFieldX(ny,nz);
        serializerToUnSerializer(
                field.getBlockSerializer(layer(nx/2,ny,nz), IndexOrdering::forward),
                localFieldX.getBlockUnSerializer(localFieldX.getBoundingBox(), IndexOrdering::forward) );
        writePpmImplementation(fName, colored, localFieldX);
    }
    if (dir == 1) {
        ScalarField2D<T> localFieldY(nx,nz);
        serializerToUnSerializer(
                field.getBlockSerializer(layer(nx,ny/2,nz), IndexOrdering::forward),
                localFieldY.getBlockUnSerializer(localFieldY.getBoundingBox(), IndexOrdering::forward) );
        writePpmImplementation(fName, colored, localFieldY);
    }
    if (dir == 2) {
        ScalarField2D<T> localFieldZ(nx,ny);
        serializerToUnSerializer(
                field.getBlockSerializer(layer(nx,ny,nz/2), IndexOrdering::forward),
                localFieldZ.getBlockUnSerializer(localFieldZ.getBoundingBox(), IndexOrdering::forward) );
        writePpmImplementation(fName, colored, localFieldZ);
    }
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writePpmImplementation (
        std::string const& fName, bool colored,
        ScalarField2D<T>& localField) const
{
    if (global::mpi().isMainProcessor()) {
        T maxVal = computeMax(localField);
        T minVal = computeMin(localField);
        if (maxVal > valueRange || minVal < 0.) {
            pcout << "ScalarField has values that exceed valueRange!" << endl
                  << "maxVal = " << maxVal << ", minVal = " << minVal << endl
                  << "No ppm written." << endl;
            return;
        }
        std::string fullName;
        if (colored)
            fullName = global::directories().getImageOutDir() + fName+".ppm";
        else
            fullName = global::directories().getImageOutDir() + fName+".pgm";
        std::ofstream fout(fullName.c_str());
        if (colored) fout << "P3\n";
        else         fout << "P2\n";
        fout << localField.getNx() << " " << localField.getNy() << "\n";
        fout << (valueRange-1) << "\n";

        for (plint iY=localField.getNy()-1; iY>=0; --iY) {
            for (plint iX=0; iX<localField.getNx(); ++iX) {
                // create variables to hold rgb value for PPM export:
                double red, green, blue;
                // red always represents a flag value of topology
                red = (double) (localField.get(iX,iY));

                switch( (int)(red) ) {
                    // however, some flags are too specific, and replaced with BB or bulk
                    case flagSedimenting:  green=blue=red= flagBB;   break;
                    // for these flags, colors are nice to spot them easily
                    case flagEroding:      green= 255.; blue=red= 0.; break;
                    case flagBuffer:       green= 0.; blue=   0.;    break;
                    case flagConstraintBB: green= 0.; blue= 255.;    break;
                    // BB and Bulk are not changed, rgb is all set to same value
                    default: green=blue=red;
                }

                if (colored) {
                    fout << (int) (red)   << " "
                         << (int) (green) << " "
                         << (int) (blue)  << "\n";
                } else {
                    fout << (int) (red)   << "\n"; // greyValue
                }
            }
        }
        fout.close();
    }
}

template<typename T>
void TopoWriter<T>::writeBinary (
        std::string const& fName,
        MultiScalarField3D<T>& field) const
{
    global::profiler().start("io");
    plint nx=field.getNx();
    plint ny=field.getNy();
    plint nz=field.getNz();
    ScalarField3D<T> localField(nx,ny,nz);
    copySerializedBlock(field, localField);
    writeBinaryImplementation(fName, localField);
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writeBinaryImplementation (
        std::string const& fName,
        ScalarField3D<T>& localField) const
{
    if (global::mpi().isMainProcessor()) {
        T maxVal = computeMax(localField);
        T minVal = computeMin(localField);
        if (maxVal > valueRange || minVal < 0.) {
            pcout << "ScalarField has values that exceed valueRange!" << endl
                  << "maxVal = " << maxVal << ", minVal = " << minVal << endl
                  << "No file written." << endl;
            return;
        }
        std::string fullName;
        fullName = global::directories().getImageOutDir() + fName+".binary";

        std::ofstream fout(fullName.c_str(), std::ios::out | std::ios::binary);

        for (plint iZ=0; iZ<localField.getNz(); ++iZ) {
            //for (plint iY=localField.getNy()-1; iY>=0; --iY) {
            for (plint iY=0; iY<localField.getNy(); ++iY) {
                for (plint iX=0; iX<localField.getNx(); ++iX) {
                    char value = localField.get(iX,iY,iZ);
                    fout.write( (char*) &value, sizeof(value) );

                    //T value = localField.get(iX,iY,iZ);
                    //fout.write( (char*) &value, sizeof(T) );
                }
            }
        }
        fout.close();
    }
}
}

template<typename T>
inRange<T>::inRange(T min_, T max_)
    : min(min_), max(max_)
{}
template<typename T>
bool inRange<T>::operator()(T value) const {
    if (value < min || value > max)
        return false;
    else
        return true;
}

template<typename T>
bool isAnyNearestNeighbourValue2D(ScalarField2D<T> field, plint iX, plint iY, T value) {
    for(int x=-1;x<=1;x+=2)
            if (field.get(iX+x,iY) == value) return true;
    for(int y=-1;y<=1;y+=2)
            if (field.get(iX,iY+y) == value) return true;
    return false;
}
template<typename T>
bool isAnyNearestNeighbourValue3D(ScalarField3D<T> field, plint iX, plint iY, plint iZ, T value) {
    for(int x=-1;x<=1;x+=2)
        if (field.get(iX+x,iY,iZ) == value) return true;
    for(int y=-1;y<=1;y+=2)
        if (field.get(iX,iY+y,iZ) == value) return true;
    for(int z=-1;z<=1;z+=2)
        if (field.get(iX,iY,iZ+z) == value) return true;
    return false;
}

template<typename T>
bool isAllNearestNeighbourValue2D(ScalarField2D<T> field, plint iX, plint iY, T value) {
    for(int x=-1;x<=1;x+=2)
        if (field.get(iX+x,iY) != value) return false;
    for(int y=-1;y<=1;y+=2)
        if (field.get(iX,iY+y) != value) return false;
    return true;
}
template<typename T>
bool isAllNearestNeighbourValue3D(ScalarField3D<T> field, plint iX, plint iY, plint iZ, T value) {
    for(int x=-1;x<=1;x+=2)
        if (field.get(iX+x,iY,iZ) != value) return false;
    for(int y=-1;y<=1;y+=2)
        if (field.get(iX,iY+y,iZ) != value) return false;
    for(int z=-1;z<=1;z+=2)
        if (field.get(iX,iY,iZ+z) != value) return false;
    return true;
}

template<typename T>
bool isAnyDirectNeighbourValue2D(ScalarField2D<T> field, plint iX, plint iY, T value) {
    for(int x=-1;x<=1;x++)
        for(int y=-1;y<=1;y++) {
            if (x==0 && y==0) continue;
            if (field.get(iX+x,iY+y) == value) return true;
        }
    return false;
}
template<typename T>
bool isAnyDirectNeighbourValue3D(ScalarField3D<T> field, plint iX, plint iY, plint iZ, T value) {
    for(int x=-1;x<=1;x++)
        for(int y=-1;y<=1;y++)
            for(int z=-1;z<=1;z++) {
                if (x==0 && y==0 && z==0) continue;
                if (field.get(iX+x,iY+y,iZ+z) == value) return true;
            }
    return false;
}

template<typename T>
bool isAllDirectNeighbourValue2D(ScalarField2D<T> field, plint iX, plint iY, T value) {
    for(int x=-1;x<=1;x++)
        for(int y=-1;y<=1;y++) {
            if (x==0 && y==0) continue;
            if (field.get(iX+x,iY+y) != value) return false;
        }
    return true;
}
template<typename T>
bool isAllDirectNeighbourValue3D(ScalarField3D<T> field, plint iX, plint iY, plint iZ, T value) {
    for(int x=-1;x<=1;x++)
        for(int y=-1;y<=1;y++)
            for(int z=-1;z<=1;z++) {
                if (x==0 && y==0 && z==0) continue;
                if (field.get(iX+x,iY+y,iZ+z) != value) return false;
            }
    return true;
}

template<typename T>
bool isAnyNextNeighbourValue3D(ScalarField3D<T> field, plint iX, plint iY, plint iZ, T value) {
    for(int x=-2;x<=2;x+=4)
        if (field.get(iX+x,iY,iZ) == value) return true;
    for(int y=-2;y<=2;y+=4)
        if (field.get(iX,iY+y,iZ) == value) return true;
    for(int z=-2;z<=2;z+=4)
        if (field.get(iX,iY,iZ+z) == value) return true;
    return false;
}


template<typename T>
void lessEqual(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<int>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_le_alpha_functional3D<T>(scalar), domain, field, result );
}


/* ******** A_le_alpha_functional3D ************************************* */

template<typename T>
A_le_alpha_functional3D<T>::A_le_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_le_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<int>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    A.get(iX,iY,iZ) <= alpha ? 1 : 0;
            }
        }
    }
}

template<typename T>
A_le_alpha_functional3D<T>* A_le_alpha_functional3D<T>::clone() const {
    return new A_le_alpha_functional3D<T>(*this);
}

template<typename T>
void A_le_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_le_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
void greaterEqual(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<int>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_ge_alpha_functional3D<T>(scalar), domain, field, result );
}


/* ******** A_ge_alpha_functional3D ************************************* */

template<typename T>
A_ge_alpha_functional3D<T>::A_ge_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_ge_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<int>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    A.get(iX,iY,iZ) >= alpha ? 1 : 0;
            }
        }
    }
}

template<typename T>
A_ge_alpha_functional3D<T>* A_ge_alpha_functional3D<T>::clone() const {
    return new A_ge_alpha_functional3D<T>(*this);
}

template<typename T>
void A_ge_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_ge_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

inline std::string createLongFileName(std::string name, plint number, plint width, std::string info, plint infoNumber) {
    std::stringstream fNameStream;
    if (infoNumber)
        fNameStream << name << std::setfill('0') << std::setw(width) << number << "_" << info << infoNumber;
    else
        fNameStream << name << std::setfill('0') << std::setw(width) << number << "_init";
    return fNameStream.str();
}


#endif
