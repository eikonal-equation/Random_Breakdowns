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
* File: CGrid.cpp
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
#include "CGrid.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <string>
#include <iostream>

/** ------ Project-specific header files -------------------------------------*/
#include "SimpleFunctions.hpp"
#include "WriteToFile.hpp"
#include "MemoryAllocations.hpp"
#include "GlobalConfiguration.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace memory;

/*==============================================================================
  Constructor
==============================================================================*/
/** Default constructor */
CGrid::CGrid(int aNx, int aNy, double aDx, double aDy) {
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;

  fValues = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
}


/**
 * The CGrid array-based constructor.
 */
CGrid::CGrid(int aNx, int aNy, double aDx, double aDy,
             std::shared_ptr<array2D_t<double>> aCostArray,
             std::shared_ptr<array2D_t<double>> aSpeedArray,
             std::shared_ptr<array2D_t<double>> aBreakdownRateArray) {
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;
  fCost = aCostArray;
  fSpeed = aSpeedArray;
  fBreakdownRate = aBreakdownRateArray;

  fValues = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
}

/**
 * The CGrid function-based constructor.
 */
CGrid::CGrid(int aNx, int aNy, double aDx, double aDy,
             std::function<double(double,double)> aCostFunction,
             std::function<double(double,double)> aSpeedFunction,
             std::function<double(double,double)> aBreakdownRateFunction) {
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;

  setCostFunctionPhysical(aCostFunction);
  setSpeedFunctionPhysical(aSpeedFunction);
  setBDRFunctionPhysical(aBreakdownRateFunction);

  fValues = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
}

/**
 * The CGrid array-based constructor for grids without a breakdown rate
 */
CGrid::CGrid(int aNx, int aNy, double aDx, double aDy,
             std::shared_ptr<array2D_t<double>> aCostArray,
             std::shared_ptr<array2D_t<double>> aSpeedArray) {
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;
  fCost = aCostArray;
  fSpeed = aSpeedArray;

  fValues = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
  fBreakdownRate = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      (*fBreakdownRate)[i][j] = 0;
    }
  }
}

/**
 * The CGrid function-based constructor for grids without breakdown rate
 */
CGrid::CGrid(int aNx, int aNy, double aDx, double aDy,
             std::function<double(double,double)> aCostFunction,
             std::function<double(double,double)> aSpeedFunction) {
  fNx = aNx;
  fNy = aNy;
  fDx = aDx;
  fDy = aDy;

  setCostFunctionPhysical(aCostFunction);
  setSpeedFunctionPhysical(aSpeedFunction);

  fValues = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
  fBreakdownRate = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));

  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      (*fBreakdownRate)[i][j] = 0;
    }
  }
}


/*==============================================================================
  Write grid to file
==============================================================================*/
void CGrid::writeGridToFile(std::string aFilename) const {
  /** Write array to file */
  io::writeToFile2D<double>(aFilename, *fValues);
}

void CGrid::readGridFromFile(std::string aFilename) const {
  /** Read array from file */
  io::readFromFile2D<double>(aFilename, *fValues);
}

/*==============================================================================
  Setters
==============================================================================*/
void CGrid::setCostFunctionPhysical(std::function<double(double,double)> aCostFunction) {
  array2D_t<double> cost_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      cost_array[i][j]  = aCostFunction(x, y);
    }
  }
  fCost = std::make_shared<array2D_t<double>>(cost_array);
}

void CGrid::setCostFunctionLogical(std::function<double(int,int)> aCostFunction) {
  array2D_t<double> cost_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      cost_array[i][j] = aCostFunction(i,j);
    }
  }
  fCost = std::make_shared<array2D_t<double>>(cost_array);
}

void CGrid::setSpeedFunctionPhysical(std::function<double(double,double)> aSpeedFunction) {
  array2D_t<double> speed_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      speed_array[i][j]  = aSpeedFunction(x, y);
    }
  }
  fSpeed = std::make_shared<array2D_t<double>>(speed_array);
}

void CGrid::setSpeedFunctionLogical(std::function<double(int,int)> aSpeedFunction) {
  array2D_t<double> speed_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      speed_array[i][j]  = aSpeedFunction(i,j);
    }
  }
  fSpeed = std::make_shared<array2D_t<double>>(speed_array);
}

void CGrid::setBDRFunctionPhysical(std::function<double(double,double)> aBreakdownRateFunction) {
  array2D_t<double> bdr_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      const double x = i * fDx;
      const double y = j * fDy;
      bdr_array[i][j]  = aBreakdownRateFunction(x, y);
    }
  }
  fBreakdownRate = std::make_shared<array2D_t<double>>(bdr_array);
}

void CGrid::setBDRFunctionLogical(std::function<double(int,int)> aBreakdownRateFunction) {
  array2D_t<double> bdr_array = allocateArray2D<double>(fNx,fNy);
  for (int i = 0; i < fNx; ++i) {
    for (int j = 0; j < fNy; ++j) {
      bdr_array[i][j]  = aBreakdownRateFunction(i,j);
    }
  }
  fBreakdownRate = std::make_shared<array2D_t<double>>(bdr_array);
}


/*==============================================================================
  Helper functions
==============================================================================*/
int CGrid::getXDirection(const int aI, const int aJ) const {
  int bestI;
  if (aI == 0) {
    bestI = 1;
  } else if (aI == fNx - 1) {
    bestI = -1;
  } else if ((*fValues)[aI-1][aJ] < (*fValues)[aI+1][aJ]) {
    bestI = -1;
  } else {
    bestI = 1;
  }

  if ((*fValues)[aI+bestI][aJ] > (*fValues)[aI][aJ]) {
    bestI = 0;
  }

  return bestI;
}

int CGrid::getYDirection(const int aI, const int aJ) const {
  int bestJ;
  if (aJ == 0) {
    bestJ = 1;
  } else if (aJ == fNy - 1) {
    bestJ = -1;
  } else if ((*fValues)[aI][aJ-1] < (*fValues)[aI][aJ+1]) {
    bestJ = -1;
  } else {
    bestJ = 1;
  }

  if ((*fValues)[aI][aJ+bestJ] > (*fValues)[aI][aJ]) {
    bestJ = 0;
  }

  return bestJ;
}

double CGrid::computeXDerivative(const int aI, const int aJ) const {
  const int xDir = getXDirection(aI, aJ);

  assert(aI + xDir < fNx);
  assert(aI + xDir >= 0);
  double xDerivative = xDir*((*fValues)[aI + xDir][aJ] - (*fValues)[aI][aJ])/fDx;

  return xDerivative;
}

double CGrid::computeYDerivative(const int aI, const int aJ) const {
  const int yDir = getYDirection(aI, aJ);
  double yDerivative;

  if ((aJ + yDir >= fNy) || (aJ + yDir < 0) || (yDir == 0)){
    yDerivative = 0;
  } else {
    assert(aJ + yDir < fNy);
    assert(aJ + yDir >= 0);
    yDerivative = yDir*((*fValues)[aI][aJ + yDir] - (*fValues)[aI][aJ])/fDy;
  }

  return yDerivative;
}
