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
 * The CombinedStatistics class -- generic implementation.
 */
#include "parallelism/mpiManager.h"
#include "parallelism/parallelStatistics.h"
#include <cmath>

namespace plb {

#ifdef PLB_MPI_PARALLEL

ParallelCombinedStatistics* ParallelCombinedStatistics::clone() const
{
    return new ParallelCombinedStatistics(*this);
}

void ParallelCombinedStatistics::reduceStatistics (
            std::vector<double>& averageObservables,
            std::vector<double>& sumWeights,
            std::vector<std::vector<double>>& listObservables,
            std::vector<double>& sumWeightsList,
            std::vector<double>& sumObservables,
            std::vector<double>& maxObservables,
            std::vector<plint>& intSumObservables ) const
{
    // Averages
    for (pluint iAverage=0; iAverage<averageObservables.size(); ++iAverage) {
        double globalAverage, globalWeight;
        global::mpi().reduce(averageObservables[iAverage]*sumWeights[iAverage], globalAverage, MPI_SUM);
        global::mpi().reduce(sumWeights[iAverage], globalWeight, MPI_SUM);
        if (global::mpi().isMainProcessor() && std::fabs(globalWeight) > 0.5) {
            globalAverage /= globalWeight;
        }
        global::mpi().bCast(&globalAverage, 1);
        averageObservables[iAverage] = globalAverage;
    }

    // List
    for (pluint iList=0; iList<listObservables.size(); ++iList) {
        int numProcesses = global::mpi().getSize();

        // Each process tells the root how many elements it holds
            // count has to be double, because gatherv_impl takes values from
            // sumWeightList[iList] which are of type double
        double countsDouble[numProcesses];
        int rank = global::mpi().getRank();
        global::mpi().gather_impl(&sumWeightsList[iList], 1, &countsDouble[0], 1, global::mpi().bossId());
        global::mpi().bCast(&countsDouble[0], numProcesses);
            // total number of to be collected values:
        double globalWeightList = 0;
        for (int i = 0; i < numProcesses; i++) globalWeightList += countsDouble[i];
        // For next steps, convert to int:
        int globalListLength = int(globalWeightList);
        int counts[numProcesses];
        for (int i = 0; i < numProcesses; i++) counts[i] = (int)countsDouble[i];
/*
        // Check if gatherv_impl worked correctly & numProceesses is correct:
        double test_globalWeightList = 0;
        global::mpi().reduce(sumWeightsList[iList], test_globalWeightList, MPI_SUM);
        if (global::mpi().isMainProcessor() && globalWeightList != test_globalWeightList ) exit(1);
*/

        // Displacements in the receive buffer for MPI_GATHERV
            // for now 'only' int, because gatherv_impl does not support long or unsigned int!!!!!
            // need to be changed for higher resolution?
        int disps[numProcesses];
        for (int i = 0; i < numProcesses; i++) disps[i] = (i > 0) ? (disps[i-1] + counts[i-1] ) : 0;
        // Collection of all list values
        std::vector<double> globalList(globalListLength);
        global::mpi().gatherv_impl(&listObservables[iList].front(), counts[rank],
                                   &globalList.front(), counts, disps, global::mpi().bossId());
        global::mpi().bCast(&globalList[0], globalListLength);
        listObservables[iList] = globalList;
    }

    // Sum
    for (pluint iSum=0; iSum<sumObservables.size(); ++iSum) {
        double globalSum;
        global::mpi().reduce(sumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        sumObservables[iSum] = globalSum;
    }

    // Max
    for (pluint iMax=0; iMax<maxObservables.size(); ++iMax) {
        double globalMax;
        global::mpi().reduce(maxObservables[iMax], globalMax, MPI_MAX);
        global::mpi().bCast(&globalMax, 1);
        maxObservables[iMax] = globalMax;
    }

    // Integer sum
    for (pluint iSum=0; iSum<intSumObservables.size(); ++iSum) {
        plint globalSum;
        global::mpi().reduce(intSumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        intSumObservables[iSum] = globalSum;
    }
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb
