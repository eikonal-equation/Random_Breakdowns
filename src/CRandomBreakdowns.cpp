/*
* ==============================================================================
*
*  Copyright (C) 2020 Marissa Gee
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
* File: CRandomBreakdowns.cpp
*
* Author: Marissa Gee
*
* Description: This is a class for solving the Random Breakdowns problem by
               iteratively solving a system of coupled PDEs.
* It has three modes: "OnlyTotalBreakdowns", "OnlyPartialBreakdowns", and "Combination".
* Each mode that requires an iterative solver can be solved by value iterations
*       value-policy iterations.
* This is determined by appending "V" or "VP" to the mode.
* This class can handle synthetic or real-world speed, running cost, and
*       breakdown rate.
*
* ==============================================================================
*/

#include "CRandomBreakdowns.hpp"

/** ------- Libraries ------------------------------------------------------- */
#include <string>
#include <iostream>
#include <cmath>
#include <boost/progress.hpp>

/** ------- Project-specific header files ----------------------------------- */
#include "MemoryAllocations.hpp"
#include "CFMM.hpp"
#include "CPolicyEvaluator.hpp"
#include "WriteToFile.hpp"
#include "GlobalConfiguration.hpp"


/** ------- Namespaces ------------------------------------------------------ */
using namespace std;
using namespace memory;

/** ============================================================================
*   Constructors
==============================================================================*/
/** Function based constructor for OnlyTotalBreakdowns (non-iterative) solver. */
CRandomBreakdowns::CRandomBreakdowns(const int aNx, const int aNy,
                  const double aDx, const double aDy,
                  double aDepotList[][2], const int aNumDepots,
                  double aGoalList[][2], const int aNumGoals,
                  function<double(double,double)> aBorderFunction,
                  function<double(double,double)> aWorkingSpeedFunction,
                  function<double(double,double)> aWorkingRunningCostFunction,
                  function<double(double,double)> aWorkingTotalBreakdownRateFunction,
                  function<double(double,double)> aRepairSpeedFunction,
                  function<double(double,double)> aRepairRunningCostFunction,
                  shared_ptr<vector<double>> aTotalRepairCostArray,
                  function<double(double,double)> aFieldRepairCostFunction) {
  /** Initialize solver parameters and arrays. */
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;
  fMode = "OnlyTotalBreakdowns";

  fDomain = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fDepots = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fGoal = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fWorkingSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fWorkingRunningCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fRepairSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fRepairRunningCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fWorkingTotalBreakdownRate = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fTotalRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fFieldRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fGoalBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      if (aBorderFunction(x,y) >= 0 ) {
        // we are in the domain (includes the boundary)
        (*fDomain)[i][j] = true;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;

        (*fWorkingSpeed)[i][j] = aWorkingSpeedFunction(x,y);
        (*fWorkingRunningCost)[i][j] = aWorkingRunningCostFunction(x,y);
        (*fRepairSpeed)[i][j] = aRepairSpeedFunction(x,y);
        (*fRepairRunningCost)[i][j] = aRepairRunningCostFunction(x,y);
        (*fWorkingTotalBreakdownRate)[i][j] = aWorkingTotalBreakdownRateFunction(x,y);
        (*fFieldRepairCost)[i][j] = aFieldRepairCostFunction(x,y);
      } else {
        // we are not in the domain
        (*fDomain)[i][j] = false;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;

        (*fWorkingSpeed)[i][j] = 0;
        (*fWorkingRunningCost)[i][j] = LARGE_NUMBER;
        (*fRepairSpeed)[i][j] = 0;
        (*fRepairRunningCost)[i][j] = LARGE_NUMBER;
        (*fWorkingTotalBreakdownRate)[i][j] = 0;
        (*fFieldRepairCost)[i][j] = 0;
      }
    }
  }

  /** Record depot locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumDepots; ++i) {
    const double depotX = aDepotList[i][0];
    const double depotY = aDepotList[i][1];
    const int depotI = round(depotX/fDx);
    const int depotJ = round(depotY/fDy);
    (*fDepots)[depotI][depotJ] = true;
    if (aTotalRepairCostArray == nullptr) {
      (*fTotalRepairCost)[depotI][depotJ] = 0;
    } else {
      (*fTotalRepairCost)[depotI][depotJ] = (*aTotalRepairCostArray)[i];
    }
  }

  /** Record goal locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumGoals; ++i) {
    const double goalX = aGoalList[i][0];
    const double goalY = aGoalList[i][1];
    const int goalI = round(goalX/fDx);
    const int goalJ = round(goalY/fDy);
    (*fGoal)[goalI][goalJ] = true;
  }
  return;
}

/** Function based constructor for OnlyPartialBreakdowns (two-mode) solver. */
CRandomBreakdowns::CRandomBreakdowns(const int aNx, const int aNy,
                  const double aDx, const double aDy, const string aMode,
                  const int aMaxIterations, const double aRho,
                  double aDepotList[][2], const int aNumDepots,
                  double aGoalList[][2], const int aNumGoals,
                  function<double(double,double)> aBorderFunction,
                  function<double(double,double)> aWorkingSpeedFunction,
                  function<double(double,double)> aWorkingRunningCostFunction,
                  function<double(double,double)> aWorkingPartialBreakdownRateFunction,
                  function<double(double,double)> aDamagedSpeedFunction,
                  function<double(double,double)> aDamagedRunningCostFunction,
                  shared_ptr<vector<double>> aPartialRepairCostArray,
                  function<double(double,double)> aWorkingTotalBreakdownRateFunction,
                  function<double(double,double)> aDamagedTotalBreakdownRateFunction,
                  function<double(double,double)> aRepairSpeedFunction,
                  function<double(double,double)> aRepairRunningCostFunction,
                  shared_ptr<vector<double>> aTotalRepairCostArray,
                  function<double(double,double)> aFieldRepairCostFunction) {
  /** Initialize solver parameters and arrays. */
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;
  fMode = aMode;
  fConvergenceFlag = false;
  fTolerance = (fDx + fDy)/2;
  fMaxIterations = aMaxIterations;
  fRho = aRho;

  fDomain = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fDepots = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fGoal = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fWorkingSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fWorkingRunningCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fWorkingPartialBreakdownRate = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fDamagedSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fDamagedRunningCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fPartialRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  fGoalBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));

  if (fMode == "TwoBreakdownTypesVP") {
    fRepairSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fRepairRunningCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fWorkingTotalBreakdownRate = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fDamagedTotalBreakdownRate = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fTotalRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fFieldRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  }

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      if (aBorderFunction(x,y) >= 0 ) {
        // we are in the domain (includes the boundary)
        (*fDomain)[i][j] = true;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;

        (*fWorkingSpeed)[i][j] = aWorkingSpeedFunction(x,y);
        (*fWorkingRunningCost)[i][j] = aWorkingRunningCostFunction(x,y);
        (*fWorkingPartialBreakdownRate)[i][j] = aWorkingPartialBreakdownRateFunction(x,y);
        (*fDamagedSpeed)[i][j] = aDamagedSpeedFunction(x,y);
        (*fDamagedRunningCost)[i][j] = aDamagedRunningCostFunction(x,y);
        (*fPartialRepairCost)[i][j] = 0;

        if (fMode == "TwoBreakdownTypesVP") {
          (*fRepairSpeed)[i][j] = aRepairSpeedFunction(x,y);
          (*fRepairRunningCost)[i][j] = aRepairRunningCostFunction(x,y);
          (*fWorkingTotalBreakdownRate)[i][j] = aWorkingTotalBreakdownRateFunction(x,y);
          (*fDamagedTotalBreakdownRate)[i][j] = aDamagedTotalBreakdownRateFunction(x,y);
          (*fTotalRepairCost)[i][j] = 0;
          (*fFieldRepairCost)[i][j] = aFieldRepairCostFunction(x,y);
        }

      } else {
        // we are not in the domain
        (*fDomain)[i][j] = false;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;

        (*fWorkingSpeed)[i][j] = 0;
        (*fWorkingRunningCost)[i][j] = LARGE_NUMBER;
        (*fWorkingPartialBreakdownRate)[i][j] = 0;
        (*fDamagedSpeed)[i][j] = 0;
        (*fDamagedRunningCost)[i][j] = LARGE_NUMBER;
        (*fPartialRepairCost)[i][j] = 0;

        if (fMode == "TwoBreakdownTypesVP") {
          (*fRepairSpeed)[i][j] = 0;
          (*fRepairRunningCost)[i][j] = LARGE_NUMBER;
          (*fWorkingTotalBreakdownRate)[i][j] = 0;
          (*fDamagedTotalBreakdownRate)[i][j] = 0;
          (*fTotalRepairCost)[i][j] = 0;
          (*fFieldRepairCost)[i][j] = 0;
        }
      }
    }
  }

  /** Record depot locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumDepots; ++i) {
    const double depotX = aDepotList[i][0];
    const double depotY = aDepotList[i][1];
    const int depotI = round(depotX/fDx);
    const int depotJ = round(depotY/fDy);
    (*fDepots)[depotI][depotJ] = true;

    if (aPartialRepairCostArray == nullptr) {
      (*fPartialRepairCost)[depotI][depotJ] = 0;
    } else {
      (*fPartialRepairCost)[depotI][depotJ] = (*aPartialRepairCostArray)[i];
    }
    if (fMode == "TwoBreakdownTypesVP") {
      if (aTotalRepairCostArray == nullptr) {
        (*fTotalRepairCost)[depotI][depotJ] = 0;
      } else {
        (*fTotalRepairCost)[depotI][depotJ] = (*aTotalRepairCostArray)[i];
      }
    }
  }

  /** Record goal locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumGoals; ++i) {
    const double goalX = aGoalList[i][0];
    const double goalY = aGoalList[i][1];
    const int goalI = round(goalX/fDx);
    const int goalJ = round(goalY/fDy);
    (*fGoal)[goalI][goalJ] = true;
  }
  return;
}

/** Constructor for real-world data examples. */
CRandomBreakdowns::CRandomBreakdowns(const int aNx, const int aNy,
                  const double aDx, const double aDy, const string aMode,
                  const int aMaxIterations, const double aRho,
                  double aDepotList[][2], const int aNumDepots,
                  double aGoalList[][2], const int aNumGoals,
                  function<double(double,double)> aBorderFunction,
                  shared_ptr<array2D_t<double>> aWorkingSpeed,
                  shared_ptr<array2D_t<double>> aWorkingRunningCost,
                  shared_ptr<array2D_t<double>> aWorkingPartialBreakdownRate,
                  shared_ptr<array2D_t<double>> aDamagedSpeed,
                  shared_ptr<array2D_t<double>> aDamagedRunningCost,
                  shared_ptr<vector<double>> aPartialRepairCostArray,
                  shared_ptr<array2D_t<double>> aWorkingTotalBreakdownRate,
                  shared_ptr<array2D_t<double>> aDamagedTotalBreakdownRate,
                  shared_ptr<array2D_t<double>> aRepairSpeed,
                  shared_ptr<array2D_t<double>> aRepairRunningCost,
                  shared_ptr<vector<double>> aTotalRepairCostArray,
                  shared_ptr<array2D_t<double>> aFieldRepairCost) {
  /** Initialize solver parameters and arrays. */
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;
  fMode = aMode;
  fConvergenceFlag = false;
  fTolerance = (fDx + fDy)/2;
  fMaxIterations = aMaxIterations;
  fRho = aRho;

  fDomain = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fDepots = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fGoal = make_shared<array2D_t<bool>>(allocateArray2D<bool>(fNx, fNy));
  fWorkingSpeed = aWorkingSpeed;
  fWorkingRunningCost = aWorkingRunningCost;
  fWorkingPartialBreakdownRate = aWorkingPartialBreakdownRate;
  fDamagedSpeed = aDamagedSpeed;
  fDamagedRunningCost = aDamagedRunningCost;
  fRepairSpeed = aRepairSpeed;
  fRepairRunningCost = aRepairRunningCost;
  fWorkingTotalBreakdownRate = aWorkingTotalBreakdownRate;
  fDamagedTotalBreakdownRate = aDamagedTotalBreakdownRate;
  fFieldRepairCost = aFieldRepairCost;

  fPartialRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  if (aMode == "TwoBreakdownTypesVP") {
    fTotalRepairCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  }
  fGoalBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      if (aBorderFunction(x,y) >= 0 ) {
        // we are in the domain (includes the boundary)
        (*fDomain)[i][j] = true;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;
        (*fPartialRepairCost)[i][j] = 0;

        if (fMode == "TwoBreakdownTypesVP") {
          (*fTotalRepairCost)[i][j] = 0;
        }
      } else {
        // we are not in the domain
        (*fDomain)[i][j] = false;
        (*fDepots)[i][j] = false;
        (*fGoal)[i][j] = false;
        (*fGoalBoundaryValue)[i][j] = 0;

        (*fWorkingSpeed)[i][j] = 0;
        (*fWorkingRunningCost)[i][j] = LARGE_NUMBER;
        (*fWorkingPartialBreakdownRate)[i][j] = 0;
        (*fDamagedSpeed)[i][j] = 0;
        (*fDamagedRunningCost)[i][j] = LARGE_NUMBER;
        (*fPartialRepairCost)[i][j] = 0;

        if (fMode == "TwoBreakdownTypesVP") {
          (*fRepairSpeed)[i][j] = 0;
          (*fRepairRunningCost)[i][j] = LARGE_NUMBER;
          (*fWorkingTotalBreakdownRate)[i][j] = 0;
          (*fDamagedTotalBreakdownRate)[i][j] = 0;
          (*fTotalRepairCost)[i][j] = 0;
          (*fFieldRepairCost)[i][j] = 0;
        }
      }
    }
  }

  /** Record depot locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumDepots; ++i) {
    const double depotX = aDepotList[i][0];
    const double depotY = aDepotList[i][1];
    const int depotI = round(depotX/fDx);
    const int depotJ = round(depotY/fDy);
    (*fDepots)[depotI][depotJ] = true;

    if (aPartialRepairCostArray == nullptr) {
      (*fPartialRepairCost)[depotI][depotJ] = 0;
    } else {
      (*fPartialRepairCost)[depotI][depotJ] = (*aPartialRepairCostArray)[i];
    }
    if (fMode == "TwoBreakdownTypesVP") {
      if (aTotalRepairCostArray == nullptr) {
        (*fTotalRepairCost)[depotI][depotJ] = 0;
      } else {
        (*fTotalRepairCost)[depotI][depotJ] = (*aTotalRepairCostArray)[i];
      }
    }
  }

  /** Record goal locations, snapping to the nearest gridpoint as necessary */
  for (int i = 0; i < aNumGoals; ++i) {
    const double goalX = aGoalList[i][0];
    const double goalY = aGoalList[i][1];
    const int goalI = round(goalX/fDx);
    const int goalJ = round(goalY/fDy);
    (*fGoal)[goalI][goalJ] = true;
  }
  return;
}

/** ============================================================================
*   Main Compute Functions
==============================================================================*/
void CRandomBreakdowns::computeExpectedCost(const string aFilename) {
  /** Always initialize the solvers */
  initializeSolvers();

  /** The PDEs being solved depends on the mode */
  if (fMode == "OnlyTotalBreakdowns") {
    /** Solve repair vehicle problem */
    fRepairAgentFMM.march();

    /** Solve fully functional problem */
    for (int i = 0; i < fNx; ++i) {
      for (int j = 0; j < fNy; ++j) {
        (*fWorkingRunningCost)[i][j] = (*fWorkingRunningCost)[i][j]
          + getWorkingTotalBreakdownRate(i,j)*fRepairValue->getValue(i,j);
      }
    }

    fWorkingAgentFMM.march();
  } else if (fMode == "OnlyPartialBreakdownsV") {
    /** Iterate FMM solvers to compute mode 1 and mode 2 value functions */
    fNumIterations = 0;
    fMaxError = fTolerance + 1;

    cout << "Progress towards max iterations"<< endl;
    boost::progress_display iter_progress(fMaxIterations);
    while ((fMaxError > fTolerance) && (fNumIterations < fMaxIterations)) {
      computeValueIteration();

      fIterationProgress.push_back(fMaxError);
      fIterations.push_back(VALUE);

      ++fNumIterations;
      ++iter_progress;
    }
    cout << endl;

    if (fMaxError <= fTolerance) {
      fConvergenceMessage = "Converged with maximum error " + to_string(fMaxError) + " after " + to_string(fNumIterations) + " iterations";
      fConvergenceFlag = true;
    } else {
      fConvergenceMessage = "Max iterations reached with total error of " + to_string(fMaxError) + " compared to tolerance " + to_string(fTolerance);
      fConvergenceFlag = false;
    }
  } else if ((fMode == "OnlyPartialBreakdownsVP") || (fMode == "TwoBreakdownTypesVP")) {
    /** Iterate FMM solvers to compute travel and repair cost */
    fNumIterations = 0;

    if (fMode == "OnlyPartialBreakdownsVP") {
      computeValueIteration();
      ++fNumIterations;
      cout << "Completed " << fNumIterations << " iteration(s)." << endl;
    }

    double d = fMaxError;

    while ((fMaxError > fTolerance) && (fNumIterations < fMaxIterations)) {
      /** Check if value iterations has stagnated */
      while ((fMaxError > fRho*d) && (fNumIterations < fMaxIterations)){
        computeValueIteration();

        fIterationProgress.push_back(fMaxError);
        fIterations.push_back(VALUE);
        ++fNumIterations;
        cout << "Completed " << fNumIterations << " iteration(s)." << endl;
      }

      if (fMaxError > fTolerance) {
        cout << "Computing policy evaluation..." << endl;
        fPolicyEval.evaluate();
        cout << "Done." << endl;

        fMaxError = 0;
        for (int i = 0; i < fNx; ++i) {
          for (int j = 0; j < fNy; ++j) {
            const double currentError = abs(fWorkingValue->getValue(i,j) - getWorkingCost(i,j));
            if (currentError >= fMaxError) {
              fMaxError = currentError;
            }

            (*fWorkingCost)[i][j] = fWorkingValue->getValue(i,j);
            (*fDepotBoundaryValue)[i][j] = fWorkingValue->getValue(i,j) + getPartialRepairCost(i,j);
          }
        }
        d = fMaxError;
        fIterationProgress.push_back(fMaxError);
        fIterations.push_back(POLICY);
        ++fNumEvaluations;
      }
    }

    if (fMaxError <= fTolerance) {
      fConvergenceMessage = "Converged with maximum error " + to_string(fMaxError) + " after " + to_string(fNumIterations) + " iterations";
      fConvergenceFlag = true;
    } else {
      fConvergenceMessage = "Max iterations reached with total error of " + to_string(fMaxError) + " compared to tolerance " + to_string(fTolerance);
      fConvergenceFlag = false;
    }
  } else {
    cout << "Invalid mode selected" << endl;
    assert(false);
  }

  writeToFile(aFilename);
  cout << fConvergenceMessage << endl;

  return;
}

/** ============================================================================
*   Helper Functions
==============================================================================*/

void CRandomBreakdowns::initializeSolvers() {
  if (fMode == "OnlyTotalBreakdowns") {
    /** Allocate appropriate grids. */
    fRepairValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fRepairRunningCost, fRepairSpeed);
    fWorkingValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fWorkingRunningCost, fWorkingSpeed);

    /** Create FMM solvers */
    fRepairAgentFMM = CFMM(fRepairValue, fDomain, fDepots);
    fWorkingAgentFMM = CFMM(fWorkingValue, fDomain, fGoal);
  } else if ((fMode == "OnlyPartialBreakdownsV") || (fMode == "OnlyPartialBreakdownsVP")) {
    /** Allocate appropriate grids. */
    fDamagedValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fDamagedRunningCost, fDamagedSpeed);
    fWorkingValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fWorkingRunningCost,
                                       fWorkingSpeed, fWorkingPartialBreakdownRate);

    // for (int i = 0; i < fNx; ++i) {
    //   for (int j = 0; j < fNy; ++j) {
    //     cout << (*fDamagedRunningCost)[i][j] << ",";
    //   }
    //   cout << endl;
    // }

    fDamagedCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fWorkingCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fDepotBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    shared_ptr<array2D_t<double>> depotDistance = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));

    /** Compute depot boundary condition by solving 1D problem */
    initializeWorkingCost(depotDistance);
    for (int i = 0; i < fNx; ++i) {
      for (int j = 0; j < fNy; ++j) {
        if ((*fDepots)[i][j] == true) {
          (*fDepotBoundaryValue)[i][j] = getPartialRepairCost(i,j) + getWorkingCost(i,j);
        }
      }
    }

    /** Create FMM solvers for both modes */
    fDamagedAgentFMM = CFMM(fDamagedValue, fDomain, fDepots, fDepotBoundaryValue);
    fWorkingAgentFMM = CFMM(fWorkingValue, fDomain, fGoal, fGoalBoundaryValue,
                            fDamagedCost);

    if (fMode == "OnlyPartialBreakdownsVP") {
      fPolicyEval = CPolicyEvaluator(fDamagedValue, fWorkingValue, fGoal, fDepots,
                                     fGoalBoundaryValue, fPartialRepairCost, fDomain);
    }
  } else if (fMode == "TwoBreakdownTypesVP") {
    /** Allocate necessary grids */
    fRepairValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fRepairRunningCost, fRepairSpeed);
    fDamagedValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fDamagedRunningCost,
                                      fDamagedSpeed, fDamagedTotalBreakdownRate);
    fWorkingValue = make_shared<CGrid>(fNx, fNy, fDx, fDy, fWorkingRunningCost,
                                       fWorkingSpeed, fWorkingPartialBreakdownRate);

    /** Solve Eikonal equation to compute repair cost for total breakdowns*/
    CFMM tempRepairFMM = CFMM(fRepairValue, fDomain, fDepots, fTotalRepairCost);
    tempRepairFMM.march();

    /** Update running cost with the expected cost of total breakdowns */
    for (int i = 0; i < fNx; ++i) {
      for (int j = 0; j < fNy; ++j) {
        (*fWorkingRunningCost)[i][j] = getWorkingRunningCost(i,j)
          + getWorkingTotalBreakdownRate(i,j) * (fRepairValue->getValue(i,j)
          + getFieldRepairCost(i,j));
        (*fDamagedRunningCost)[i][j] = getDamagedRunningCost(i,j)
          + getDamagedTotalBreakdownRate(i,j) * (fRepairValue->getValue(i,j)
          + getFieldRepairCost(i,j));
      }
    }

    fDamagedCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fWorkingCost = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    fDepotBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
    shared_ptr<array2D_t<double>> depotDistance = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));

    /** Compute depot boundary condition by solving 1D problem */
    initializeWorkingCost(depotDistance);
    for (int i = 0; i < fNx; ++i) {
      for (int j = 0; j < fNy; ++j) {
        (*fDepotBoundaryValue)[i][j] = getPartialRepairCost(i,j) + getWorkingCost(i,j);
      }
    }

    /** Compute worst case cost using FMM solvers and modified problem */
    fDamagedAgentFMM = CFMM(fDamagedValue, fDomain, fDepots, fDepotBoundaryValue);
    fWorkingAgentFMM = CFMM(fWorkingValue, fDomain, fGoal, fGoalBoundaryValue,
                            fDamagedCost);
    computeValueIteration();

    /** Create FMM solvers for both modes */
    fDamagedAgentFMM = CFMM(fDamagedValue, fDomain, fDepots, fDepotBoundaryValue,
                           fWorkingCost);
    fPolicyEval = CPolicyEvaluator(fDamagedValue, fWorkingValue, fGoal, fDepots,
                                   fGoalBoundaryValue, fPartialRepairCost, fDomain);
  } else {
    cout << "Invalid mode selected" << endl;
    assert(false);
  }

  return;
}

void CRandomBreakdowns::initializeWorkingCost(shared_ptr<array2D_t<double>> aDepotDistance) {
  /** Compute worst-case parameters for 1D case */
  goalDepotDistance(aDepotDistance);
  double phiMax = 0;
  double kMaxW= 0;
  double kMaxB= 0;
  double fMinW = LARGE_NUMBER;
  double fMinB = LARGE_NUMBER;

  for (int i = 0; i < fNx; i++) {
    for (int j = 0; j < fNy; j++) {
      if ((*fDomain)[i][j] == true) {
        const double kW = getWorkingRunningCost(i,j);
        const double kB = getDamagedRunningCost(i,j);
        const double phi = getWorkingPartialBreakdownRate(i,j);
        const double fW = getWorkingSpeed(i,j);
        const double fB = getDamagedSpeed(i,j);

        if (phiMax < phi) {
          phiMax = phi;
        }
        if (kMaxW < kW) {
          kMaxW = kW;
        }
        if (kMaxB < kB) {
          kMaxB = kB;
        }
        if ((fMinW > fW) && (fW > 0)) {
          fMinW = fW;
        }
        if ((fMinB > fB) && (fB > 0)) {
          fMinB = fB;
        }
      }
    }
  }

  if (phiMax == 0) {
    phiMax = 0.001;
  }

  /** Use solution to the ODE to determine mode 1 value function at each depot */
  for (int i = 0; i < fNx; i++) {
    for (int j = 0; j < fNy; j++) {
      if((*fDepots)[i][j] == true) {
        const double L = (*aDepotDistance)[i][j];
        const double RD = getPartialRepairCost(i,j);

        (*fWorkingCost)[i][j] = (kMaxW/phiMax + (kMaxB*fMinW)/(phiMax*fMinB) + RD)
                                * (exp((phiMax/fMinW)*L) - 1)
                                - (kMaxB/fMinB)*L;
      }
    }
  }
}

void CRandomBreakdowns::goalDepotDistance(shared_ptr<array2D_t<double>> aDepotDistance) {
  /** Compute distance between depots and goal use FMM with constant speed */
  shared_ptr<array2D_t<double>> tempConstant2D = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx, fNy));
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      (*tempConstant2D)[i][j] = 1;
    }
  }

  shared_ptr<CGrid> distanceValueFunction = make_shared<CGrid>(fNx, fNy, fDx, fDy, tempConstant2D, tempConstant2D);
  CFMM distanceFMM = CFMM(distanceValueFunction, fDomain, fGoal);
  distanceFMM.march();

  /** Record the minimum distance to goal from each depot */
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      if ((*fDepots)[i][j] == true) {
        (*aDepotDistance)[i][j] = distanceValueFunction->getValue(i,j);
      }
    }
  }
}

void CRandomBreakdowns::computeValueIteration() {
  /** Solve mode 2 PDE and update value function*/
  fDamagedAgentFMM.march();

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      (*fDamagedCost)[i][j] = fDamagedValue->getValue(i,j);
    }
  }

  /** Solve mode 1 PDE and update value function */
  fWorkingAgentFMM.march();

  /** Update change between iterations */
  fMaxError = 0;
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double currentError = abs(fWorkingValue->getValue(i,j) - getWorkingCost(i,j));
      if (currentError >= fMaxError) {
        fMaxError = currentError;
      }

      (*fWorkingCost)[i][j] = fWorkingValue->getValue(i,j);
      (*fDepotBoundaryValue)[i][j] = fWorkingValue->getValue(i,j) + getPartialRepairCost(i,j);
    }
  }
}

/** ============================================================================
*   Write-to-file Functions
==============================================================================*/

void CRandomBreakdowns::writeToFile(const std::string aFilename) const {
  /** Grid sizes */
  std::vector<int> grid_sizes = {fNx, fNy};
  io::writeVectorToFile<int>(aFilename + "_GridSizes", grid_sizes);
  std::vector<double> step_sizes = {fDx, fDy};
  io::writeVectorToFile<double>(aFilename + "_StepSizes", step_sizes);

  /** Record iteration information */
  io::writeVectorToFile<iteration_t>(aFilename + "_IterationType", fIterations);
  io::writeVectorToFile<double>(aFilename + "_IterationProgress", fIterationProgress);

  /** Write arrays to file */
  io::writeToFile2D<bool>(aFilename + "_Domain", *fDomain);
  io::writeToFile2D<bool>(aFilename + "_Goal", *fGoal);
  io::writeToFile2D<bool>(aFilename + "_Depots", *fDepots);

  if (fMode == "OnlyTotalBreakdowns") {
    io::writeToFile2D<double>(aFilename + "_WorkingTotalBreakdownRate", *fWorkingTotalBreakdownRate);
    io::writeToFile2D<double>(aFilename + "_WorkingSpeed", *fWorkingSpeed);
    io::writeToFile2D<double>(aFilename + "_RepairSpeed", *fRepairSpeed);
    io::writeToFile2D<double>(aFilename + "_RepairRunningCost", *fRepairRunningCost);
    io::writeToFile2D<double>(aFilename + "_WorkingRunningCost", *fWorkingRunningCost);
    io::writeToFile2D<double>(aFilename + "_FieldRepairCost", *fFieldRepairCost);
    fRepairValue->writeGridToFile(aFilename + "_RepairCost");
    fWorkingValue->writeGridToFile(aFilename + "_WorkingCost");
  } else if ((fMode == "OnlyPartialBreakdownsV") || (fMode == "OnlyPartialBreakdownsVP")) {
    io::writeToFile2D<double>(aFilename + "_WorkingPartialBreakdownRate", *fWorkingPartialBreakdownRate);
    io::writeToFile2D<double>(aFilename + "_WorkingSpeed", *fWorkingSpeed);
    io::writeToFile2D<double>(aFilename + "_DamagedSpeed", *fDamagedSpeed);
    io::writeToFile2D<double>(aFilename + "_DamagedRunningCost", *fDamagedRunningCost);
    io::writeToFile2D<double>(aFilename + "_WorkingRunningCost", *fWorkingRunningCost);

    fDamagedValue->writeGridToFile(aFilename + "_DamagedCost");
    fWorkingValue->writeGridToFile(aFilename + "_WorkingCost");
  } else {
    io::writeToFile2D<double>(aFilename + "_WorkingTotalBreakdownRate", *fWorkingTotalBreakdownRate);
    io::writeToFile2D<double>(aFilename + "_WorkingPartialBreakdownRate", *fWorkingPartialBreakdownRate);
    io::writeToFile2D<double>(aFilename + "_DamagedTotalBreakdownRate", *fDamagedTotalBreakdownRate);
    io::writeToFile2D<double>(aFilename + "_WorkingSpeed", *fWorkingSpeed);
    io::writeToFile2D<double>(aFilename + "_DamagedSpeed", *fDamagedSpeed);
    io::writeToFile2D<double>(aFilename + "_RepairSpeed", *fRepairSpeed);
    io::writeToFile2D<double>(aFilename + "_RepairRunningCost", *fRepairRunningCost);
    io::writeToFile2D<double>(aFilename + "_WorkingRunningCost", *fWorkingRunningCost);
    io::writeToFile2D<double>(aFilename + "_DamagedRunningCost", *fDamagedRunningCost);
    io::writeToFile2D<double>(aFilename + "_PartialRepairCost", *fPartialRepairCost);
    io::writeToFile2D<double>(aFilename + "_TotalRepairCost", *fTotalRepairCost);
    io::writeToFile2D<double>(aFilename + "_FieldRepairCost", *fFieldRepairCost);

    fRepairValue->writeGridToFile(aFilename + "_RepairCost");
    fDamagedValue->writeGridToFile(aFilename + "_DamagedCost");
    fWorkingValue->writeGridToFile(aFilename + "_WorkingCost");
  }

  return;
}
