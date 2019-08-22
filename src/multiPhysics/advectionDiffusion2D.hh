/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ADVECTION_DIFFUSION_2D_HH
#define ADVECTION_DIFFUSION_2D_HH

#include "multiPhysics/advectionDiffusion2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

#include <math.h>
#include <random>

namespace plb {

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>::LatticeToPassiveAdvDiff2D(T scaling_)
    : scaling(scaling_)
{ }

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>::process (
              Box2D domain,
              BlockLattice2D<T,FluidDescriptor>& fluid,
              BlockLattice2D<T,ScalarDescriptor>& scalar )
{
    Dot2D offset = computeRelativeDisplacement(fluid, scalar);
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T *u = scalar.get(iX+offset.x,iY+offset.y).getExternal(velOffset);
            Array<T,2> velocity;
            fluid.get(iX,iY).computeVelocity(velocity);
            velocity *= scaling;
            velocity.to_cArray(u);
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>*
    LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>::clone() const
{
    return new LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>(*this);
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
BlockDomain::DomainT LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void latticeToPassiveAdvDiff (
        MultiBlockLattice2D<T,FluidDescriptor>& fluid, MultiBlockLattice2D<T,ScalarDescriptor>& scalar, Box2D domain )
{
    applyProcessingFunctional (
            new LatticeToPassiveAdvDiff2D<T,FluidDescriptor,ScalarDescriptor>(), domain, fluid, scalar );
}


template<typename T, template<typename U> class Descriptor>
void BoxTemperatureFromScalarFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T temperature = scalarField.get(iX+offset.x,iY+offset.y);
            // pcout << temperature << " at " << iX << ", " << iY << std::endl;
            lattice.get(iX,iY).defineDensity(temperature);
            iniCellAtEquilibrium(lattice.get(iX,iY), temperature, Array<T,2>((T)0.,(T)0.));
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BoxTemperatureFromScalarFunctional2D<T,Descriptor>* BoxTemperatureFromScalarFunctional2D<T,Descriptor>::clone() const
{
    return new BoxTemperatureFromScalarFunctional2D<T,Descriptor>(*this);
}
//
template<typename T, template<typename U> class Descriptor>
void BoxTemperatureFromScalarFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    // modified[1] = modif::dataStructure;
}
//
template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT BoxTemperatureFromScalarFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor>
BoxTemperatureRateFunctional2D<T,Descriptor>::BoxTemperatureRateFunctional2D(T temperaturRate_)
    : temperatureRate(temperaturRate_)
{ }

template<typename T, template<typename U> class Descriptor>
void BoxTemperatureRateFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T temperature = lattice.get(iX,iY).computeDensity();
            T changedTemperature = temperature + temperatureRate;
            lattice.get(iX,iY).defineDensity(changedTemperature);
            iniCellAtEquilibrium(lattice.get(iX,iY),(T)changedTemperature, Array<T,2> ((T)0.,(T)0.));
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BoxTemperatureRateFunctional2D<T,Descriptor>* BoxTemperatureRateFunctional2D<T,Descriptor>::clone() const
{
    return new BoxTemperatureRateFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void BoxTemperatureRateFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT BoxTemperatureRateFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_2D_HH

