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
* File: CPolicyEvaluator.hpp
*
* Author: Marissa Gee
*
* Description: This is a class for solving a system of linear PDEs using finite
*           differences
* It works on a 2D regular grid on a square domain.
* It works for systems of 2 equations.
* It uses the CGrid and Eigen::Sparse::SparseLU classes.
*
* ==============================================================================
*/

#ifndef CPOLICYEVALUATOR_HPP
#define CPOLICYEVALUATOR_HPP

/** ------- Libraries ------------------------------------------------------- */
#include <string>
#include <Eigen/Sparse>

/** ------- Project-specific header files ----------------------------------- */
#include "CFMM.hpp"
#include "CGrid.hpp"
#include "MemoryAllocations.hpp"



class CPolicyEvaluator
{
  public:
    /** ========================================================================
      Constructors
    ==========================================================================*/
    /** Default Constructor */
    CPolicyEvaluator() = default;

    /**
     * The CPolicyEvaluator constructor
     * @param aDamagedValue  pointer to CGrid object containing the value
     *           function for mode 2
     * @param aTravelValue  pointer to CGrid object containing the value
     *           function for mode 1
     * @param aGoal  pointer to array defining the goal locations
     * @param aDepots  pointer to array defining depot locations
     * @param aGoalBoundaryValue  pointer to array defining boundary condition
     *           at the goal
     * @param aDepotBoundaryValue  pointer to array defining boundary condition
     *           at the depots
     * @param aDomain  pointer to array defining the domain on which we are
     *           solving the PDEs
     */
    CPolicyEvaluator(const std::shared_ptr<CGrid> aDamagedValue,
                     const std::shared_ptr<CGrid> aTravelValue,
                     const std::shared_ptr<memory::array2D_t<bool>> aGoal,
                     const std::shared_ptr<memory::array2D_t<bool>> aDepots,
                     const std::shared_ptr<memory::array2D_t<double>> aGoalBoundaryValue,
                     const std::shared_ptr<memory::array2D_t<double>> aDepotBoundaryValue,
                     const std::shared_ptr<memory::array2D_t<bool>> aDomain);

    /** Computes the PDE solution by solving finite difference equations
     *  throughout the domain. The matrix of equation is solved via sparse LU
     *  factorization.
     */
    void evaluate();

  private:
    /** Grid parameters */
    int fNx;
    int fNy;
    double fDx;
    double fDy;

    /** Matrix parameters */
    int fMatrixSize;

    /** A shared pointer to a CGrid for the value function in mode 2 */
    std::shared_ptr<CGrid> fDamagedValue;
    /** A shared pointer to a CGrid for the value function in mode 1 */
    std::shared_ptr<CGrid> fWorkingValue;
    /** Array representing Domain */
    std::shared_ptr<memory::array2D_t<bool>> fDomain;
    /** Array representing goal locations */
    std::shared_ptr<memory::array2D_t<bool>> fGoal;
    /** Array representing depot locations */
    std::shared_ptr<memory::array2D_t<bool>> fDepots;
    /** Array representing boundary condition at goal */
    std::shared_ptr<memory::array2D_t<double>> fGoalBoundaryValue;
    /** Array representing boundary condition at depot */
    std::shared_ptr<memory::array2D_t<double>> fDepotBoundaryValue;

    /** Finite differences initialization.
     *  This function constructs a matrix as a vector of value-position triples
     *      and initializes the right-hand side of the finite differences
     *      equations.
     * @param aEntries  pointer to a vector of value-position triples
     * @param aRHS  pointer to a vector for storing the right-hand side of the
     *      finite differences equations.
     */
    void initializeFDSystem(std::shared_ptr<std::vector<Eigen::Triplet<double>>> aEntries, std::shared_ptr<Eigen::VectorXd> aRHS);

    /** A helper function which computes the x component of the control at
    *       (aI, aJ).
    * @param aI  int x logical coordinate of grid point
    * @param aI  int y logical coordinate of grid point
    * @param aMode  string representing whether we are in mode 1 or 2
    */
    double computeXControl(const int aI, const int aJ, const std::string aMode) const;

    /** A helper function which computes the y component of the control at
    *       (aI, aJ).
    * @param aI  int x logical coordinate of grid point
    * @param aI  int y logical coordinate of grid point
    * @param aMode  string representing whether we are in mode 1 or 2
    */
    double computeYControl(const int aI, const int aJ, const std::string aMode) const;
};

#endif
