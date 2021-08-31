/*
* ==============================================================================
*
*  Copyright (C) 2021 Marissa Gee
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------------
*
* File: CPolicyEvaluator.cpp
*
* Author: Marissa Gee
*
* Description:
*
* ==============================================================================
*/


#include "CPolicyEvaluator.hpp"

/** ------- Libraries ------------------------------------------------------- */
#include <string>
#include <iostream>
#include <cmath>
#include <Eigen/Sparse>

/** ------- Project-specific header files ----------------------------------- */
#include "MemoryAllocations.hpp"
#include "CFMM.hpp"
#include "WriteToFile.hpp"
#include "GlobalConfiguration.hpp"


/** ------- Namespaces ------------------------------------------------------ */
using namespace std;
using namespace memory;
using namespace Eigen;

/** ============================================================================
*   Constructors
==============================================================================*/
CPolicyEvaluator::CPolicyEvaluator(const shared_ptr<CGrid> aDamagedValue,
                                   const shared_ptr<CGrid> aTravelValue,
                                   const shared_ptr<array2D_t<bool>> aGoal,
                                   const shared_ptr<array2D_t<bool>> aDepots,
                                   const shared_ptr<array2D_t<double>> aGoalBoundaryValue,
                                   const shared_ptr<array2D_t<double>> aDepotBoundaryValue,
                                   const shared_ptr<array2D_t<bool>> aDomain){
  fDamagedValue = aDamagedValue;
  fWorkingValue = aTravelValue;
  fGoal = aGoal;
  fDepots = aDepots;
  fDomain = aDomain;
  fGoalBoundaryValue = aGoalBoundaryValue;
  fDepotBoundaryValue = aDepotBoundaryValue;

  fNx = fDamagedValue->getGridSizeX();
  fNy = fDamagedValue->getGridSizeY();
  fDx = fDamagedValue->getDx();
  fDy = fDamagedValue->getDy();
  fMatrixSize = 2*fNx*fNy;
};

/** ============================================================================
*   Main Compute Function
==============================================================================*/
void CPolicyEvaluator::evaluate() {
  /** Allocate vectors for FD matrix entries and rhs*/
  shared_ptr<vector<Triplet<double>>> entries = make_shared<vector<Triplet<double>>>();
  shared_ptr<VectorXd> rhs = make_shared<VectorXd>(VectorXd(fMatrixSize));
  VectorXd soln(fMatrixSize);

  initializeFDSystem(entries, rhs);
  SparseMatrix<double> matrix(fMatrixSize, fMatrixSize);
  matrix.setFromTriplets(entries->begin(), entries->end());

  /** Solve FD system via Sparse LU factorization */
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  solver.compute(matrix);

  if (solver.info() != Success) {
    cout << "LU Decomposition Failed" << endl;
    assert(solver.info() == Success);
  }
  soln = solver.solve(*rhs);

  /** Update value functions with FD solutions */
  for (int i = 1; i < fNx-1; i++) {
    for (int j = 1; j < fNy-1; j++) {
      fDamagedValue->setValue(i, j, soln[i*fNy+j]);
      fWorkingValue->setValue(i, j, soln[i*fNy+j+fNx*fNy]);
    }
  }
};

void CPolicyEvaluator::initializeFDSystem(shared_ptr<vector<Triplet<double>>> aEntries, shared_ptr<VectorXd> aRHS) {
  /** Initialize the component of the matrix corresponding to mode 2*/
  for (int k = 0; k < fNx*fNy; k++) {
    /** Convert indices from vector to matrix form */
    const int rI = floor(k/fNy);
    const int rJ = k % fNy;

    if ((*fDomain)[rI][rJ] == false) {
      /** Set value function to infinity outside the domain*/
      aEntries->push_back(Triplet<double>(k, k, 1));
      (*aRHS)[k] = INF;
    } else {
      if ((*fDepots)[rI][rJ] == false) {
        const int xDir = fDamagedValue->getXDirection(rI, rJ);
        const int yDir = fDamagedValue->getYDirection(rI, rJ);
        const double xControl = computeXControl(rI, rJ, "Damaged");
        const double yControl = computeYControl(rI, rJ, "Damaged");

        const double K = fDamagedValue->getCost(rI, rJ);
        const double f = fDamagedValue->getSpeed(rI, rJ);
        const double lambda = fDamagedValue->getBreakdownRate(rI, rJ);

        /** Compute upwind discretization and store in FD matrix */
        aEntries->push_back(Triplet<double>(k, k, -1*xDir*xControl/fDx
                                                    - yDir*yControl/fDy
                                                    + lambda/f));
        aEntries->push_back(Triplet<double>(k, k + yDir, yDir*yControl/fDy));
        aEntries->push_back(Triplet<double>(k, k + xDir*fNy, xDir*xControl/fDx));

        /** Term corresponding to coupling with s */
        aEntries->push_back(Triplet<double>(k, k + fNy*fNx, -lambda/f));

        (*aRHS)[k] = K/f;
      } else {
        /** Enforce the boundary condition at the depots */
        aEntries->push_back(Triplet<double>(k, k, 1));
        aEntries->push_back(Triplet<double>(k, k + fNy*fNx, -1));
        (*aRHS)[k] = (*fDepotBoundaryValue)[rI][rJ];
      }
    }
  }

  /** Initialize the component of the matrix corresonding to mode 1 */
  for (int k = 0; k < fNx*fNy; k++) {
    /** Convert indices from vector to matrix form */
    const int sI = floor(k/fNy);
    const int sJ = k % fNy;
    const int Mk = k + fNy*fNx;

    if ((*fDomain)[sI][sJ] == false) {
      /** Set value function to infinity outside the domain*/
      aEntries->push_back(Triplet<double>(Mk, Mk, 1));
      (*aRHS)[Mk] = INF;
    } else {
      if ((*fGoal)[sI][sJ] == false) {
        const int xDir = fWorkingValue->getXDirection(sI, sJ);
        const int yDir = fWorkingValue->getYDirection(sI, sJ);

        const double xControl = computeXControl(sI, sJ, "Working");
        const double yControl = computeYControl(sI, sJ, "Working");

        const double phi = fWorkingValue->getBreakdownRate(sI, sJ);
        const double f = fWorkingValue->getSpeed(sI, sJ);
        const double K = fWorkingValue->getCost(sI, sJ);

        /** Compute upwind discretization and store in FD matrix */
        aEntries->push_back(Triplet<double>(Mk, Mk, -1*xDir*xControl/fDx
                                                      - yDir*yControl/fDy
                                                      + phi/f));
        aEntries->push_back(Triplet<double>(Mk, Mk + yDir, yDir*yControl/fDy));
        aEntries->push_back(Triplet<double>(Mk, Mk + xDir*fNy, xDir*xControl/fDx));

        /** Term corresponding to coupling with r */
        aEntries->push_back(Triplet<double>(Mk, Mk-fNx*fNy, -phi/f));

        (*aRHS)[Mk] = K/f;
      } else {
        /** Enforce the boundary condition at the goal */
        aEntries->push_back(Triplet<double>(Mk, Mk, 1));
        (*aRHS)[Mk] = (*fGoalBoundaryValue)[sI][sJ];
      }
    }
  }
}

/** ============================================================================
*  Helper Functions
==============================================================================*/
double CPolicyEvaluator::computeXControl(const int aI, const int aJ, const string aMode) const {
  double xDerivative;
  double yDerivative;
  if (aMode == "Working") {
    xDerivative = fWorkingValue->computeXDerivative(aI, aJ);
    yDerivative = fWorkingValue->computeYDerivative(aI, aJ);
  } else if (aMode == "Damaged") {
    xDerivative = fDamagedValue->computeXDerivative(aI, aJ);
    yDerivative = fDamagedValue->computeYDerivative(aI, aJ);
  } else {
    cout << "Invalid mode for control" << endl;
    assert(false);
  }

  double xControl;
  if (xDerivative == 0) {
    xControl = 0;
  } else {
    xControl = xDerivative/pow((pow(xDerivative,2) + pow(yDerivative, 2)), 0.5);
  }

  return xControl;
}

double CPolicyEvaluator::computeYControl(const int aI, const int aJ, const string aMode) const {
  double xDerivative;
  double yDerivative;
  if (aMode == "Working") {
    xDerivative = fWorkingValue->computeXDerivative(aI, aJ);
    yDerivative = fWorkingValue->computeYDerivative(aI, aJ);
  } else if (aMode == "Damaged") {
    xDerivative = fDamagedValue->computeXDerivative(aI, aJ);
    yDerivative = fDamagedValue->computeYDerivative(aI, aJ);
  } else {
    cout << "Invalid mode for control" << endl;
    assert(false);
  }

  double yControl;
  if (yDerivative == 0) {
    yControl = 0;
  } else {
    yControl = yDerivative/pow((pow(xDerivative,2) + pow(yDerivative, 2)), 0.5);
  }

  return yControl;
}
