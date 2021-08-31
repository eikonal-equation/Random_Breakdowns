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
* File: GlobalConfiguration.hpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
*
* Description: This file contains the definitions of various numerical constants
* and types. It also contains all global parameters for all of the numerics.
*
* ==============================================================================
*/

#ifndef GLOBALCONFIG_HPP
#define GLOBALCONFIG_HPP

#include <functional>
#include <limits>

/*-----------------------------------------------------------------------------/
/-- Typedefs
/-----------------------------------------------------------------------------*/
enum status_t {FAR, CONSIDERED, ACCEPTED};
enum region_t {BORDER, DOMAIN, INITIAL};
enum iteration_t {VALUE, POLICY};


/*----------------------------------------------------------------------------//
//-- Global Constants
//----------------------------------------------------------------------------*/
/** Mathematical constants */
constexpr double INF = std::numeric_limits<double>::max();
constexpr double LARGE_NUMBER = 1.0e32;
constexpr double PI = 3.141592653589793;
constexpr double SQRT2 = 1.41421356237309504880168872420969807;

#endif
