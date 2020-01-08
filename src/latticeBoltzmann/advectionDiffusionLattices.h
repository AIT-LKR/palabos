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

/** \file
 * Descriptors for advection-diffusion physics. In principle, thanks
 * to the fact that the Palabos code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef ADVECTION_DIFFUSION_LATTICES_H
#define ADVECTION_DIFFUSION_LATTICES_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/externalFields.h"
#include "latticeBoltzmann/roundOffPolicy.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include <vector>

namespace plb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors 
{
//===========================================================================//
//=================== AdvectionDiffusion Lattice Descriptors=================//
//===========================================================================//
    
    /// D2Q5 lattice
    template <typename T> struct D2Q5Constants
    {
        enum { d = 2, q = 5 };      ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
    };

    template <typename T> struct D2Q5DescriptorBase
        : public D2Q5Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D2Q5DescriptorBase<T> BaseDescriptor;
        enum { numPop=D2Q5Constants<T>::q };
    };
    
    /// AD D2Q5 lattice
    template <typename T> struct AdvectionDiffusionD2Q5Descriptor
        : public D2Q5DescriptorBase<T>, public Velocity2dDescriptorBase
    {
        static const char name[];
    };
    
    
    template <typename T> struct AdvectionDiffusionWithSourceD2Q5Descriptor 
    : public D2Q5DescriptorBase<T>, public VelocityAndScalar2dBase
    {
        static const char name[];
    };
    
    /// AD D2Q9 lattice
    template <typename T> struct AdvectionDiffusionD2Q9Descriptor
        : public D2Q9DescriptorBase<T>, public Velocity2dDescriptorBase
    {
        static const char name[];
    };
    
    
    template <typename T> struct AdvectionDiffusionWithSourceD2Q9Descriptor 
    : public D2Q9DescriptorBase<T>, public VelocityAndScalar2dBase
    {
        static const char name[];
    };
    
    /// D3Q7 lattice
    template <typename T> struct D3Q7Constants {
        enum { d = 3, q = 7 };     ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
    };

    template <typename T> struct D3Q7DescriptorBase
        : public D3Q7Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D3Q7DescriptorBase<T> BaseDescriptor;
        enum { numPop=D3Q7Constants<T>::q };
    };

    template <typename T> struct AdvectionDiffusionD3Q7Descriptor 
    : public D3Q7DescriptorBase<T>, public Velocity3dBase
    {
        static const char name[];
    };

    template <typename T> struct AdvectionDiffusionWithSourceD3Q7Descriptor 
    : public D3Q7DescriptorBase<T>, public VelocityAndScalar3dBase
    {
        static const char name[];
    };

    /// D3Q19 lattice
    template <typename T> struct AdvectionDiffusionD3Q19Descriptor 
    : public D3Q19DescriptorBase<T>, public Velocity3dBase
    {
        static const char name[];
    };

    template <typename T> struct AdvectionDiffusionWithSourceD3Q19Descriptor 
    : public D3Q19DescriptorBase<T>, public VelocityAndScalar3dBase
    {
        static const char name[];
    };

    /// D3Q27 lattice
    template <typename T> struct AdvectionDiffusionD3Q27Descriptor 
    : public D3Q27DescriptorBase<T>, public Velocity3dBase
    {
        static const char name[];
    };

    template <typename T> struct AdvectionDiffusionWithSourceD3Q27Descriptor 
    : public D3Q27DescriptorBase<T>, public VelocityAndScalar3dBase
    {
        static const char name[];
    };

}  // namespace descriptors

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_LATTICES_H

