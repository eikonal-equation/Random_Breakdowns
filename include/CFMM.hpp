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
* File: CFMM.hpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
*
* Description: This is a class for solving stationary PDEs
*              using the Fast Marching Method.
* It has two modes: "Eikonal" and "RandomlyTerminated".
* It works on a 2D regular grid and a 5 point stencil.
* It allows for different spacing in the horizontal and vertical directions.
* The target set is determined by aDomain.
* The solution to the PDE is stored in the fPrimary grid.
* The FMM implementation uses the boost::heap and CGrid classes.
* It contains an augmented grid to compute the evolution of the boundary
*
* ==============================================================================
*/

#ifndef CFMM_HPP
#define CFMM_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <boost/heap/binomial_heap.hpp>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CGrid.hpp"

class CFMM
{
  public:
    /** ========================================================================
      Constructors
    ==========================================================================*/
    /** Default Constructor */
    CFMM() = default;

    /**
     * The Eikonal CFMM Constructor (RepairToYou)
     * @param aPrimary pointer to CGrid object containing the primary grid
     *           on which the Eikonal equation is solved
     * @param aDomain  pointer to array defining domain in which to solve
     *           Eikonal equation
     * @param aBoundarySet  pointer to array defining location of the goal set,
     *           assumes value of zero at boundary
     */
    CFMM(const std::shared_ptr<CGrid> aPrimary,
         const std::shared_ptr<memory::array2D_t<bool>> aDomain,
         const std::shared_ptr<memory::array2D_t<bool>> aBoundarySet);

     /**
      * The Eikonal CFMM Constructor (YouToRepair)
      * @param aPrimary pointer to CGrid object containing the primary grid
      *           on which the Eikonal equation is solved
      * @param aDomain  pointer to array defining domain in which to solve the
      *           Eikonal equation
      * @param aBoundarySet  pointer to array defining location of the goal set
      * @param aBoundaryValue pointer to array defining the boundary condition
      *           at the goal set
      */
     CFMM(const std::shared_ptr<CGrid> aPrimary,
          const std::shared_ptr<memory::array2D_t<bool>> aDomain,
          const std::shared_ptr<memory::array2D_t<bool>> aBoundarySet,
          const std::shared_ptr<memory::array2D_t<double>> aBoundaryValue);


    /**
    * The Randomly Terminated CFMM Constructor (YouToRepair, Combination)
    * @param aPrimary pointer to CGrid object containing the primary grid
    *           on which the Eikonal equation is solved
    * @param aDomain  pointer to array defining domain in which to solve the
    *           Eikonal equation
    * @param aBoundarySet  pointer to array defining location of the goal set
    * @param aBoundaryValue pointer to array defining the boundary condition
    *           at the goal set
    * @param aTerminalCost pointer to array defining the cost of breaking down
    */
    CFMM(const std::shared_ptr<CGrid> aPrimary,
         const std::shared_ptr<memory::array2D_t<bool>> aDomain,
         const std::shared_ptr<memory::array2D_t<bool>> aBoundarySet,
         const std::shared_ptr<memory::array2D_t<double>> aBoundaryValue,
         const std::shared_ptr<memory::array2D_t<double>> aTerminalCost);


    /** ========================================================================
      Compute functions
    ==========================================================================*/
    /**
     * This function computes the PDE solution using the Fast Marching Method.
     *
     * The solution is stored in fPrimary (accessible by getValue/setValue).
     */
    void march();

    /**
     * This function computes the policy improvement step of policy iteration
     *
     * The solution is stored in fPrimary (accessible by getValue/setValue).
     */
    void improve();

  private:
    /** A GridPoint Class to be used by the heap. */
    class CHeapGP {
      public:
        int fI;
        int fJ;
        double fValue;
        CHeapGP(int aI, int aJ, double aValue): fI(aI), fJ(aJ), fValue(aValue) {};
    };

    /** Struct comparison which is required by boost::heap */
    struct compare_CHeapGP
    {
      bool operator()(const CHeapGP& aPoint1, const CHeapGP& aPoint2) const
      {
        return aPoint1.fValue > aPoint2.fValue;
      }
    };

    /** Typedef Heap types to make it easier to read. */
    typedef boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> > CFMMHeap_t;
    typedef typename boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> >::handle_type handle_t;

    /**
     * A string defining the type of PDE being solved, Eikonal or
     * RandomlyTerminated
     */
    std::string fMode;

    /** A shared pointer to a CGrid for the solution of the Eikonal equation */
    std::shared_ptr<CGrid> fPrimary;

    /** Array representing Domain */
    std::shared_ptr<memory::array2D_t<bool>> fDomain;
    /** Array representing the locations the of the goal set */
    std::shared_ptr<memory::array2D_t<bool>> fBoundarySet;
    /** Array representing the value function at the boundary set*/
    std::shared_ptr<memory::array2D_t<double>> fBoundaryValue;
    /**
     * Array representing the cost of breakdown in the randomly terminated
     * problem
     */
    std::shared_ptr<memory::array2D_t<double>> fTerminalCost;

    /** The status (far/considered/accepted) of each grid point. */
    std::shared_ptr<memory::array2D_t<status_t>> fStatus;
    /** Backpointers for the heap */
    std::shared_ptr<memory::array2D_t<handle_t>> fHeapPointers;

    /**
     * FMM Initialization.
     * This function sets the status of all nodes
     *      and adds the border to the heap.
     * @param CFMMHeap_t& aCFMMHeap Boost::heap passed by value
     *      it is initialized by this function.
     */
    void initialize_cfmm(CFMMHeap_t& aCFMMHeap);

    /**
    * A helper function which updates neighbors and adds them to the heap.
    * @param aCFMMHeap boost heap object
    * @param aCurrent_i int x logical coordinate of grid point
    * @param aCurrent_j int y logical coordinate of grid point
    */
    void update_neighbors(CFMMHeap_t& aCFMMheap,
                          const int aCurrent_i, const int aCurrent_j);

    /**
     * A helper function to figure out if grid point (aI,aJ) is in the domain.
     * @param aI int x logical coordinate
     * @param aJ int y logical coordinate
     */
    bool in_domain(const int aI, const int aJ) const;

    /**
     * A helper function to figure out if grid point (aI,aJ) is a local minimum
     * of the terminal cost.
     *
     * @param aI int x logical coordinate
     * @param aJ int y logical coordinate
     */
    bool is_local_minimum(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two horizontal neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    int smaller_h_neighbor(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two vertical neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    int smaller_v_neighbor(const int aI, const int aJ) const;

    /** Computes the Eikonal update.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    * @param aVerticalNeighbor int direction of smaller vertical neighbor
    * @param aHorizontalNeighbor int direction of smaller horizontal neighbor
    */
    double compute_from_two_neighbors(const int aI, const int aJ,
                                      const int aVerticalNeighbor,
                                      const int aHorizontalNeighbor);

    /** Update a gridpoint. Compute the update, assign, and return if updated.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    bool update_gp(const int aI, const int aJ);
};

/** ============================================================================
  Inline function definitions
==============================================================================*/

/* Compute if gridpoint [aI,aJ] is inside the domain */
inline
bool CFMM::in_domain(const int aI, const int aJ) const {
  if (aI >= 0 && aI < fPrimary->getGridSizeX() && aJ >=0 && aJ < fPrimary->getGridSizeY()) {
    return (*fDomain)[aI][aJ];
  }
  return false;
}

#endif
