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
* File: CRandomBreakdowns.hpp
*
* Author: Marissa Gee
*
* Description: This is a class for solving the Random Breakdowns problem by
               iteratively solving a system of coupled PDEs.
* It has three modes: "OnlyTotalBreakdowns", "OnlyPartialBreakdowns", and "TwoBreakdownTypes".
* Each mode that requires an iterative solver can be solved by value iterations
*       value-policy iterations.
* This is determined by appending "V" or "VP" to the mode.
* This class can handle synthetic or real-world speed, running cost, and
*       breakdown rate.
*
* ==============================================================================
*/

#ifndef CRANDOMBREAKDOWNS_HPP
#define CRANDOMBREAKDOWNS_HPP


/** ------- Libraries ------------------------------------------------------- */
#include <string>

/** ------- Project-specific header files ----------------------------------- */
#include "CFMM.hpp"
#include "CPolicyEvaluator.hpp"
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"


/** ------- Typedefs -------------------------------------------------------- */


class CRandomBreakdowns
{
  public:
    /** ========================================================================
    *   Constructors
    ==========================================================================*/

    /** Default Constructor */
    CRandomBreakdowns() = default;

    /**
     * Function based constructor for OnlyTotalBreakdowns (non-iterative) solver.
     * Inputs are grid parameters, functions for environmental parameters
     *           defined thoughout the entire domain, and arrays for quantities
     *           defined only at the depots or goal.
     */
    CRandomBreakdowns(const int aNx, const int aNy,
                      const double aDx, const double aDy,
                      double aDepotList[][2], const int aNumDepots,
                      double aGoalList[][2], const int aNumGoals,
                      std::function<double(double,double)> aBorderFunction,
                      std::function<double(double,double)> aWorkingSpeedFunction,
                      std::function<double(double,double)> aWorkingRunningCostFunction,
                      std::function<double(double,double)> aWorkingTotalBreakdownRateFunction = &zero_constant2D,
                      std::function<double(double,double)> aRepairSpeedFunction = &constant2D,
                      std::function<double(double,double)> aRepairRunningCostFunction = &constant2D,
                      std::shared_ptr<std::vector<double>> aTotalRepairCostArray = nullptr,
                      std::function<double(double,double)> aFieldRepairCostFunction = &zero_constant2D);

    /**
     * Function based constructor for TwoBreakdownTypes (two-mode) solver.
     * Inputs are grid parameters, functions for environmental parameters
     *           defined thoughout the entire domain, and arrays for quantities
     *           defined only at the depots or goal.
     */
    CRandomBreakdowns(const int aNx, const int aNy,
                      const double aDx, const double aDy, const std::string aMode,
                      const int aMaxIterations, const double aRho,
                      double aDepotList[][2], const int aNumDepots,
                      double aGoalList[][2], const int aNumGoals,
                      std::function<double(double,double)> aBorderFunction,
                      std::function<double(double,double)> aWorkingSpeedFunction,
                      std::function<double(double,double)> aWorkingRunningCostFunction,
                      std::function<double(double,double)> aWorkingPartialBreakdownRateFunction,
                      std::function<double(double,double)> aDamagedSpeedFunction,
                      std::function<double(double,double)> aDamagedRunningCostFunction,
                      std::shared_ptr<std::vector<double>> aPartialRepairCostArray = nullptr,
                      std::function<double(double,double)> aWorkingTotalBreakdownRateFunction = &zero_constant2D,
                      std::function<double(double,double)> aDamagedTotalBreakdownRateFunction = &zero_constant2D,
                      std::function<double(double,double)> aRepairSpeedFunction = &constant2D,
                      std::function<double(double,double)> aRepairRunningCostFunction = &constant2D,
                      std::shared_ptr<std::vector<double>> aTotalRepairCostArray = nullptr,
                      std::function<double(double,double)> aFieldRepairCostFunction = &zero_constant2D);

    /**
     * Constructor for real-world data examples.
     * Inputs are grid parameters, arrays for quantities based on terrain data.
     */
    CRandomBreakdowns(const int aNx, const int aNy,
                      const double aDx, const double aDy, const std::string aMode,
                      const int aMaxIterations, const double aRho,
                      double aDepotList[][2], const int aNumDepots,
                      double aGoalList[][2], const int aNumGoals,
                      std::function<double(double,double)> aBorderFunction,
                      std::shared_ptr<memory::array2D_t<double>> aWorkingSpeed,
                      std::shared_ptr<memory::array2D_t<double>> aWorkingRunningCost,
                      std::shared_ptr<memory::array2D_t<double>> aWorkingPartialBreakdownRate,
                      std::shared_ptr<memory::array2D_t<double>> aDamagedSpeed,
                      std::shared_ptr<memory::array2D_t<double>> aDamagedRunningCost,
                      std::shared_ptr<std::vector<double>> aPartialRepairCostArray = nullptr,
                      std::shared_ptr<memory::array2D_t<double>> aWorkingTotalBreakdownRate = nullptr,
                      std::shared_ptr<memory::array2D_t<double>> aDamagedTotalBreakdownRate = nullptr,
                      std::shared_ptr<memory::array2D_t<double>> aRepairSpeed = nullptr,
                      std::shared_ptr<memory::array2D_t<double>> aRepairRunningCost = nullptr,
                      std::shared_ptr<std::vector<double>> aTotalRepairCostArray = nullptr,
                      std::shared_ptr<memory::array2D_t<double>> aFieldRepairCost = nullptr);

    /**
     * Main compute function.
     * Computes the expected cost-to-go at all points in the domain.
     * @param aFilename  string filename prefix used for storing results.
     */
    void computeExpectedCost(const std::string aFilename);

    /** ========================================================================
    *    Getters
    * ========================================================================*/
    double getTolerance() const;
    int getMaxIterations() const;
    double getError() const;
    double getRho() const;
    int getNumIterations() const;
    int getNumEvaluations() const;

    double getWorkingSpeed(const int aI, const int aJ) const;
    double getDamagedSpeed(const int aI, const int aJ) const;
    double getRepairSpeed(const int aI, const int aJ) const;
    double getWorkingRunningCost(const int aI, const int aJ) const;
    double getDamagedRunningCost(const int aI, const int aJ) const;
    double getRepairRunningCost(const int aI, const int aJ) const;
    double getWorkingTotalBreakdownRate(const int aI, const int aJ) const;
    double getWorkingPartialBreakdownRate(const int aI, const int aJ) const;
    double getDamagedTotalBreakdownRate(const int aI, const int aJ) const;
    double getWorkingCost(const int aI, const int aJ) const;
    double getDamagedCost(const int aI, const int aJ) const;
    double getTotalRepairCost(const int aI, const int aJ) const;
    double getPartialRepairCost(const int aI, const int aJ) const;
    double getFieldRepairCost(const int aI, const int aJ) const;

    /** ========================================================================
    *    Setters
    *=========================================================================*/
    void setMaxIterations(const int aMaxIterations);

    /**
     * Main write results function.
     * Writes problem parameters and solution to files.
     * @param aFilename  string the file prefix used when writing to file.
     */
    void writeToFile(const std::string aFilename) const;


  private:
    /** Grid Parameters */
    int fNx;
    int fNy;
    double fDx;
    double fDy;

    /** Solver parameters */
    std::string fMode;
    double fTolerance;
    int fMaxIterations;
    double fRho;

    /** Convergence information */
    double fMaxError;
    int fNumIterations;
    int fNumEvaluations;
    std::vector<iteration_t> fIterations;
    std::vector<double> fIterationProgress;
    bool fConvergenceFlag;
    std::string fConvergenceMessage;

    /** ========================================================================
    *    Arrays
    *=========================================================================*/
    /** Arrays storing speed for the different modes */
    std::shared_ptr<memory::array2D_t<double>> fWorkingSpeed;
    std::shared_ptr<memory::array2D_t<double>> fDamagedSpeed;
    std::shared_ptr<memory::array2D_t<double>> fRepairSpeed;

    /** Arrays storing running cost for the different modes */
    std::shared_ptr<memory::array2D_t<double>> fRepairRunningCost;
    std::shared_ptr<memory::array2D_t<double>> fDamagedRunningCost;
    std::shared_ptr<memory::array2D_t<double>> fWorkingRunningCost;

    /** Arrays storing breakdown rates for the different modes */
    std::shared_ptr<memory::array2D_t<double>> fWorkingTotalBreakdownRate;
    std::shared_ptr<memory::array2D_t<double>> fWorkingPartialBreakdownRate;
    std::shared_ptr<memory::array2D_t<double>> fDamagedTotalBreakdownRate;

    /** Arrays storing repair costs at depots and in the field */
    std::shared_ptr<memory::array2D_t<double>> fTotalRepairCost;
    std::shared_ptr<memory::array2D_t<double>> fPartialRepairCost;
    std::shared_ptr<memory::array2D_t<double>> fFieldRepairCost;

    /** Arrays storing boundary values at goal and depots */
    std::shared_ptr<memory::array2D_t<double>> fGoalBoundaryValue;
    std::shared_ptr<memory::array2D_t<double>> fDepotBoundaryValue;

    /** Arrays storing the value function for the different modes */
    std::shared_ptr<memory::array2D_t<double>> fWorkingCost;
    std::shared_ptr<memory::array2D_t<double>> fDamagedCost;

    /** Arrays storing location of domain, goal, and depots */
    std::shared_ptr<memory::array2D_t<bool>> fDomain;
    std::shared_ptr<memory::array2D_t<bool>> fDepots;
    std::shared_ptr<memory::array2D_t<bool>> fGoal;

    /** ========================================================================
    *    Solvers
    *=========================================================================*/
    /** FMM solver for the repair vehicle problem */
    CFMM fRepairAgentFMM;
    /** FMM solver for mode 1, fully functional */
    CFMM fWorkingAgentFMM;
    /** FMM solver for mode 2, partially broken down */
    CFMM fDamagedAgentFMM;
    /** FD solver for computing policy evaluations */
    CPolicyEvaluator fPolicyEval;

    /** ========================================================================
    *    Grids
    *=========================================================================*/
    /** CGrid storing the value function for the repair vehicle problem */
    std::shared_ptr<CGrid> fRepairValue;
    /** CGrid storing the value function for mode 1, fully functiona */
    std::shared_ptr<CGrid> fWorkingValue;
    /** CGrid storing the value function for mode 2, partially broken down */
    std::shared_ptr<CGrid> fDamagedValue;

    /** ========================================================================
    *   Helper Functions
    ==========================================================================*/
    /**
     * Random breakdowns solver initialization.
     * This function intializes the necessary FMM solvers based on the mode, as
     *      as initializing the value functions throughout the domain
     */
    void initializeSolvers();

    /**
     * Mode 1 value function initialization.
     * Computes the mode 1 value function at the depot locations using the
     *      solution to an ODE.
     */
    void initializeWorkingCost(std::shared_ptr<memory::array2D_t<double>> aDepotDistance);

    /**
     * Computes the maximum distance between a goal and a depot in the domain.
     */
    void goalDepotDistance(std::shared_ptr<memory::array2D_t<double>> aDepotDistance);
    /**
     * Computes one instance of value iteration, using FMM solvers and updating
     *      appropriate arrays.
     */
    void computeValueIteration();
    /**
     * Computes one instance of value iteration, using FMM solvers and updating
     *      appropriate arrays. Not currently implemented.
     */
    void computePolicyIteration();
};

/** ============================================================================
*   Inline Function Definitions
==============================================================================*/
inline
double CRandomBreakdowns::getTolerance() const {
  return fTolerance;
}

inline
int CRandomBreakdowns::getMaxIterations() const {
  return fMaxIterations;
}

inline
double CRandomBreakdowns::getError() const {
  return fMaxError;
}

inline
double CRandomBreakdowns::getRho() const {
  return fRho;
}

inline
int CRandomBreakdowns::getNumIterations() const {
  return fNumIterations;
}

inline
int CRandomBreakdowns::getNumEvaluations() const {
  return fNumEvaluations;
}

inline
void CRandomBreakdowns::setMaxIterations(const int aMaxIterations) {
  fMaxIterations = aMaxIterations;
}

inline
double CRandomBreakdowns::getWorkingSpeed(const int aI, const int aJ) const {
  return (*fWorkingSpeed)[aI][aJ];
}

inline
double CRandomBreakdowns::getDamagedSpeed(const int aI, const int aJ) const {
  return (*fDamagedSpeed)[aI][aJ];
}

inline
double CRandomBreakdowns::getRepairSpeed(const int aI, const int aJ) const {
  return (*fRepairSpeed)[aI][aJ];
}

inline
double CRandomBreakdowns::getWorkingRunningCost(const int aI, const int aJ) const {
  return (*fWorkingRunningCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getDamagedRunningCost(const int aI, const int aJ) const {
  return (*fDamagedRunningCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getRepairRunningCost(const int aI, const int aJ) const {
  return (*fRepairRunningCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getWorkingTotalBreakdownRate(const int aI, const int aJ) const {
  return (*fWorkingTotalBreakdownRate)[aI][aJ];
}

inline
double CRandomBreakdowns::getWorkingPartialBreakdownRate(const int aI, const int aJ) const {
  return (*fWorkingPartialBreakdownRate)[aI][aJ];
}

inline
double CRandomBreakdowns::getDamagedTotalBreakdownRate(const int aI, const int aJ) const {
  return (*fDamagedTotalBreakdownRate)[aI][aJ];
}

inline
double CRandomBreakdowns::getWorkingCost(const int aI, const int aJ) const {
  return (*fWorkingCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getDamagedCost(const int aI, const int aJ) const {
  return (*fDamagedCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getTotalRepairCost(const int aI, const int aJ) const {
  return (*fTotalRepairCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getPartialRepairCost(const int aI, const int aJ) const {
  return (*fPartialRepairCost)[aI][aJ];
}

inline
double CRandomBreakdowns::getFieldRepairCost(const int aI, const int aJ) const {
  return (*fFieldRepairCost)[aI][aJ];
}
#endif
