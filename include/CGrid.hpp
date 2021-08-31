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
* File: CGrid.hpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
*
* Description: This class is used to handle both the logical and physical
* representation of a 2D regular grid with fixed spacing in both dimensions.
* This class is an interface between the underlying
* data structures and the rest of the code.
* The grid values, cost function, speed function, and breakdown rates are stored
*           as arrays.
*
* ==============================================================================
*/

#ifndef CGRID_HPP
#define CGRID_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"
#include "GlobalConfiguration.hpp"

class CGrid
{
  private:
    /** Grid parameters */
    int fNx;
    int fNy;
    double fDx;
    double fDy;

    /* Grid of values. Stored as a 2D boost::multi_array */
    std::shared_ptr<memory::array2D_t<double>> fValues;

    /** Running cost function, stored as a 2D boost::multi_array.
    * This has to be given as an input.
    * Computes cost K at physical location (x,y) */
    std::shared_ptr<memory::array2D_t<double>> fCost;

    /** Speed function, stored as a 2D boost::multi_array.
    * This has to be given as an input.
    * Computes speed f at physical location (x,y) */
    std::shared_ptr<memory::array2D_t<double>> fSpeed;

    /** Breakdown rate function, stored as a 2D boost::multi_array.
    * This has to be given as an input.
    * Computes breakdown rate phi or lambda at physical location (x,y) */
    std::shared_ptr<memory::array2D_t<double>> fBreakdownRate;

  public:
    /** ========================================================================
    *    Constructors
    * ========================================================================*/

    /** Default constructor */
    CGrid(int aNx, int aNy, double aDx, double aDy);

    /**
     * The CGrid array-based constructor.
     */
    CGrid(int aNx, int aNy, double aDx, double aDy,
          std::shared_ptr<memory::array2D_t<double>> aCostArray,
          std::shared_ptr<memory::array2D_t<double>> aSpeedArray,
          std::shared_ptr<memory::array2D_t<double>> aBreadownRateArray);

    /**
     * The CGrid function-based constructor.
     */
    CGrid(int aNx, int aNy, double aDx, double aDy,
          std::function<double(double,double)> aCostFunction,
          std::function<double(double,double)> aSpeedFunction,
          std::function<double(double,double)> aBreakdownRateArray);

    /**
     * The CGrid array-based constructor.
     */
    CGrid(int aNx, int aNy, double aDx, double aDy,
          std::shared_ptr<memory::array2D_t<double>> aCostArray,
          std::shared_ptr<memory::array2D_t<double>> aSpeedArray);

    /**
     * The CGrid function-based constructor.
     */
    CGrid(int aNx, int aNy, double aDx, double aDy,
          std::function<double(double,double)> aCostFunction,
          std::function<double(double,double)> aSpeedFunction);

    /** ========================================================================
    *    Setters
    *=========================================================================*/
    void setValue(const int aI, const int aJ, const double aValue);
    void setCostArray(std::shared_ptr<memory::array2D_t<double>> aCostArray);
    void setSpeedArray(std::shared_ptr<memory::array2D_t<double>> aSpeedArray);
    void setBDRArray(std::shared_ptr<memory::array2D_t<double>> aBreakdownRateArray);
    void setCostFunctionPhysical(std::function<double(double,double)> aCostFunction);
    void setCostFunctionLogical(std::function<double(int,int)> aCostFunction);
    void setSpeedFunctionPhysical(std::function<double(double,double)> aSpeedFunction);
    void setSpeedFunctionLogical(std::function<double(int,int)> aSpeedFunction);
    void setBDRFunctionPhysical(std::function<double(double,double)> aBreakdownRateFunction);
    void setBDRFunctionLogical(std::function<double(int,int)> aBreakdownRateFunction);

    /** ========================================================================
    *    Getters
    * ========================================================================*/
    /** Logical-coordinate functions */
    double getValue(const int aI, const int aJ) const;
    double getCost(const int aI, const int aJ) const;
    double getSpeed(const int aI, const int aJ) const;
    double getBreakdownRate(const int aI, const int aJ) const;
    int getXDirection(const int aI, const int aJ) const;
    int getYDirection(const int aI, const int aJ) const;

    int getGridSizeX() const;
    int getGridSizeY() const;

    double getMaxX() const;
    double getMaxY() const;

    double getDx() const;
    double getDy() const;

    double computeXDerivative(const int aI, const int aJ) const;
    double computeYDerivative(const int aI, const int aJ) const;

    /** ========================================================================
    *    Other
    * ========================================================================*/
    /** Mapping back and forth between grid and physical */
    double xGridToPhysical(const int aI) const;
    double yGridToPhysical(const int aJ) const;

    /**
     * A I/O member which prints the value and cost grids to file.
     * @param aFilename a string which contains the prefix to the names
     *      of the files to which the grids will be printed.
     * The (cost/value) grids will be printed to
     *      files called "aFilename"+(Cost/Speed/Value)
     */
    void writeGridToFile(const std::string aFilename) const;

    /**
     * A I/O member which reads the value grid from a file.
     * @param aFilename a string which contains the prefix to the names
     *      of the files to which the grids will be printed.
     * The value grid will be updated with the contents of aFilename
     */
    void readGridFromFile(const std::string aFilename) const;
};

/** ============================================================================
*    Inline function definitions
* ============================================================================*/
/** Inline definitions must be in the header file.
 * These functions are used frequently. */

/** ------ Inline definition of setters ------------------------------*/
inline void CGrid::setValue(const int aI, const int aJ, const double aValue) {
  (*fValues)[aI][aJ] = aValue;
}

inline void CGrid::setCostArray(std::shared_ptr<memory::array2D_t<double>> aCostArray) {
  fCost = aCostArray;
}

inline void CGrid::setSpeedArray(std::shared_ptr<memory::array2D_t<double>> aSpeedArray) {
  fSpeed = aSpeedArray;
}

inline void CGrid::setBDRArray(std::shared_ptr<memory::array2D_t<double>> aBreakdownRateArray) {
  fBreakdownRate = aBreakdownRateArray;
}

/** ------ Inline definition of getters --------------------------------------*/
inline double CGrid::getValue(const int aI, const int aJ) const {
  return (*fValues)[aI][aJ];
}

inline double CGrid::getCost(const int aI, const int aJ) const {
  return (*fCost)[aI][aJ];
}

inline double CGrid::getSpeed(const int aI, const int aJ) const {
  return (*fSpeed)[aI][aJ];
}

inline double CGrid::getBreakdownRate(const int aI, const int aJ) const {
  return (*fBreakdownRate)[aI][aJ];
}

inline int CGrid::getGridSizeX() const {
  return fNx;
}

inline int CGrid::getGridSizeY() const {
  return fNy;
}

inline double CGrid::getMaxX() const {
  return fDx * (fNx-1);
}

inline double CGrid::getMaxY() const {
  return fDy * (fNy-1);
}

inline double CGrid::getDx() const {
  return fDx;
}

inline double CGrid::getDy() const {
  return fDy;
}

inline double CGrid::xGridToPhysical(const int aI) const {
  return (double)aI * fDx;
}

inline double CGrid::yGridToPhysical(const int aJ) const {
  return (double)aJ * fDy;
}


#endif
