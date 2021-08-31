/*
* ==============================================================================
*
*  Copyright (C) 2019  Elliot Cartee
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
* File: SimpleFunctions.hpp
*
* Authors: Elliot Cartee
*   (based on code by Marc Aurèle Gilles)
*
* Modified by Marissa Gee (3/8/10)
*
* Description: This file contains inline definitions of many small/simple
* functions used for speed or patrol density functions.
*
* ==============================================================================
*/

#ifndef SIMPLEFUNCTIONS_HPP
#define SIMPLEFUNCTIONS_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <vector>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"


/*==============================================================================
    Constant functions
==============================================================================*/
inline double constant1D(double x) {
  return 1.0;
}

inline double constant2D(double x, double y) {
  return 1.0;
}

inline double zero_constant2D(double x, double y) {
  return 0.0;
}

inline double var_constant2D(double x, double y, const double c) {
  return c;
}


/*==============================================================================
    Goal/Depot Functions
==============================================================================*/

inline double bar_at_zero(double x, double y) {
	return -0.5 + exp(-100*pow(x,2));
}

inline double bar_at_one(double x, double y) {
	return -0.5 + exp(-100*pow(x-1.0,2));
}

inline double point_gaussian(double x, double y, const double x_center, const double y_center) {
	return -0.5 + exp(-100*(pow(x-x_center,2) + pow(y-y_center,2)));
}

/*==============================================================================
    Speed functions
==============================================================================*/

inline double linear_speed(double x, double y) {
	return (1-x) + y + 1;
}

inline double speed_obstacle(double x, double y) {
	double f;
	if (0.25*0.25 - ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) > 0) {
		f = 0.0;
	} else {
		f = 1.0;
	}
	return f;
}

inline double horizontal_gaussian(double x, double y) {
	return exp(-10*pow(x-0.5,2));
}

inline double slow_horizontal_gaussian(double x, double y) {
	return 0.1*exp(-10*pow(x-0.5,2));
}

/*==============================================================================
    Benefit functions
==============================================================================*/
inline double circle_benefit(double x, double y) {
  return 2.0;
}

inline double circle_benefit_modified(double x, double y) {
  const double r = sqrt(pow(0.5-x,2) + pow(0.5-y,2));
  return 3.0 * (1 - 2*r);
}

inline double rectangle_benefit(double x, double y) {
  const double d1 = fmin(fmin(x, 3.0-x), fmin(y, 1.0-y));
  return 16*d1*(1.0 - d1);
}

/*==============================================================================
		Running cost functions
==============================================================================*/

/*==============================================================================
    Breakdown Chance Functions
==============================================================================*/
inline double gaussian(double x, double y, const double x_center, const double y_center) {
	return 7*exp(-10*(pow(x-x_center,2) + pow(y-y_center,2)));
}

inline double linear_breakdown_rate(double x, double y, const double phi) {
	return phi*((1-x) + y);
}

inline double circle_band_patrol(double x, double y) {
  const double d1 = sqrt(pow(x-0.5,2) + pow(y-0.5,2));
  if (d1 > 0.15 && d1 < 0.35) {
    return 1.0;
  }
  return 0.0;
}

inline double square_patrol(double x, double y) {
  const double d1 = fmin(fmin(x, 1.0-x), fmin(y, 1.0-y));
  return 1.0/(50*(d1-0.3)*(d1-0.3) + 0.5);
}

/*==============================================================================
    Domain functions
==============================================================================*/
inline double square_border(double x, double y) {
  if (x == 0 || x == 1 || y == 0 || y == 1) {
    return -1;
  } else {
    return 1;
  }
}

inline double rectangle_border(double x, double y) {
  if (x == 0 || x == 3 || y == 0 || y == 1) {
    return -1;
  } else {
    return 1;
  }
}

inline double circle_border(double x, double y) {
	return 0.5*0.5 - ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
}


#endif
