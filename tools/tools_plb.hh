#ifndef TOOLS_PLB_HH
#define TOOLS_PLB_HH

#include <../src/core/globalDefs.h>
#include <../src/io/parallelIO.h>
#include <../externalLibraries/tinyxml/tinyxml.h>

#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>

#include "./tools.h"
#include "./tools.hh"
#include "./tools_geometric.h"
#include "./tools_geometric.hh"
#include "./tools_files.h"
#include "./tools_files.hh"
#include "./plb_physicalFlowParam.h"
#include "./tools_plb.h"


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

/* ******** MaskedReplaceAlphaFunctional3D ************************************** */
// based on: MaskedBoxScalarMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_SS<T,int>
// in: dataAnalysisFunctional3D.h:512

template<typename T1, typename T2>
MaskedReplaceAlphaFunctional3D<T1,T2>::MaskedReplaceAlphaFunctional3D(T1 alpha_, std::vector<int>& flags_)
      : alpha(alpha_),
        flags(flags_)
{ }

template<typename T1, typename T2>
MaskedReplaceAlphaFunctional3D<T1,T2>::MaskedReplaceAlphaFunctional3D(T1 alpha_, int flag_)
      : alpha(alpha_)
{ flags.push_back(flag_); }

template<typename T1, typename T2>
void MaskedReplaceAlphaFunctional3D<T1,T2>::process (
        Box3D domain,
        ScalarField3D<T1>& scalarField,
        ScalarField3D<T2>& mask )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (unsigned int i=0; i<flags.size(); i++) {
                    if (mask.get(iX+offset.x, iY+offset.y, iZ+offset.z)==flags[i]) {
                        scalarField.get(iX,iY,iZ) = alpha;
                    }
                }
            }
        }
    }
}

/* ******** MaskedBoxScalarListFunctional3D ******* */
template<typename T, class BoolMask>
void MaskedBoxScalarListFunctional3D<T, BoolMask>::process (Box3D domain,
                                                       ScalarField3D<T>& scalarField,
                                                       ScalarField3D<int>& mask ) {
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (condition(mask.get(iX+offset.x, iY+offset.y, iZ+offset.z))) {
                    statistics.gatherList(listScalarId, (T)scalarField.get(iX,iY,iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T, class BoolMask>
std::vector<T> MaskedBoxScalarListFunctional3D<T, BoolMask>::getListScalar() const {
    std::vector<T> list = this->getStatistics().getList(listScalarId);
    // The list is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return std::vector<T>(list.begin(), list.end());
    }
    return list;
}

/* ************* Class PseudomaskedSmoothen3D ******************* */

template<typename T>
void PseudomaskedSmoothen3D<T>::process(Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result)
{
    Dot3D offset = computeRelativeDisplacement(data, result);
    //T weightDiagonalNeighbour = .25;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T *res = &result.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                *res = data.get(iX,iY,iZ) * weightCentreCell ;
                if (*res == 0.)
                    continue; // skip this cell, because it's not part of fluid
                float n = weightCentreCell;
                for (int i = -1; i < 2; i++) {
                    plint nextX = iX + i;
                    for (int j = -1; j < 2; j++) {
                        plint nextY = iY + j;
                        for (int k = -1; k < 2; k++) {
                            plint nextZ = iZ + k;
                            if (!(i == 0 && j == 0 && k == 0)) {
                                T currentCell = data.get(nextX, nextY, nextZ);
                                if (currentCell != 0.) {
                                    n+= weightDirectNeighbour;
                                    *res += currentCell * weightDirectNeighbour;
                                }
                            }
                        }
                    }
                }
                *res /= (T) n;
            }
        }
    }
}

template<typename T>
PseudomaskedSmoothen3D<T>* PseudomaskedSmoothen3D<T>::clone() const
{
    return new PseudomaskedSmoothen3D<T>(*this);
}

template<typename T>
void PseudomaskedSmoothen3D<T>::getTypeOfModification (std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* *************** PseudomaskedSmoothen3D ******************************************* */

template<typename T>
void pseudomaskedSmoothen(MultiScalarField3D<T>& data, MultiScalarField3D<T>& result, Box3D domain, T weightDirectNeighbour)
{
    applyProcessingFunctional(new PseudomaskedSmoothen3D<T>(weightDirectNeighbour), domain, data, result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > pseudomaskedSmoothen(MultiScalarField3D<T>& data, Box3D domain, T weightDirectNeighbour)
{
    std::auto_ptr<MultiScalarField3D<T> > result =
        generateMultiScalarField<T>(data, domain);
    pseudomaskedSmoothen<T>(data, *result, domain, weightDirectNeighbour);
    return result;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > pseudomaskedSmoothen(MultiScalarField3D<T>& data, T weightDirectNeighbour)
{
    return pseudomaskedSmoothen<T>(data, data.getBoundingBox(), weightDirectNeighbour);
}

/* **************************************************************************** */

template<typename T1, typename T2>
PoiseuilleVelocity<T1,T2>::PoiseuilleVelocity(T2 uLB_, T2 uDev_, T1 powerPoiseuilleVel_, char dir_,
        T1 inletCentre_, T1 inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
     inletCentreA(inletCentre_), inletCentreB(-1), inletRadius(inletRadius_),
     seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T2>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T2>(0., uDev);
}

template<typename T1, typename T2>
PoiseuilleVelocity<T1,T2>::PoiseuilleVelocity(T2 uLB_, T2 uDev_, T1 powerPoiseuilleVel_, char dir_,
        T1 inletCentreA_, T1 inletCentreB_, T1 inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
      inletCentreA(inletCentreA_), inletCentreB(inletCentreB_), inletRadius(inletRadius_),
      seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T2>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T2>(0., uDev);
}

template<typename T1, typename T2>
void PoiseuilleVelocity<T1,T2>::operator()(T1 iX, T1 iY, T1 iZ, Array<T2,3>& u) const {
    T2 main_velocity = uLB;
    T2 side_velocity = 0.;
    T2 down_velocity = 0.;
    if (uDev) {
        main_velocity = (*main_distribution)(*generator);
        side_velocity = (*side_distribution)(*generator)/4;
        down_velocity = (*side_distribution)(*generator)/4;
    }
    T1 r;

    if (dir == 'x') {
        if (inletCentreB == -1)
            r = distanceFromPoint<T2>(iY,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint<T2>(iY,iZ, inletCentreA,inletCentreB);
        u[0] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[1] = side_velocity;
        u[2] = down_velocity;
    }
    if (dir == 'y') {
        if (inletCentreB == -1)
            r = distanceFromPoint<T2>(iX,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint<T2>(iX,iZ, inletCentreA,inletCentreB);
        u[0] = side_velocity;
        u[1] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[2] = down_velocity;
    }
}

template<typename T1, typename T2>
void PoiseuilleVelocity<T1,T2>::operator()(T1 iX, T1 iY, Array<T2,2>& u) const {
    T2 main_velocity = uLB;
    T2 side_velocity = 0.;
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


#endif
