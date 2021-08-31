/*
* ==============================================================================
*
*  Copyright (C) 2021  Marissa Gee
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
* File: CFMM.cpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
**
* Description: This is a class for solving stationary Eikonal equations using
* the Fast Marching Method.
* It works on a 2D regular grid and a 5 point stencil.
* It allows for different spacing in the horizontal and vertical directions.
* The target set is determined by aDomain.
* The solution to the Eikonal is stored in the fPrimary grid.
* The FMM implementation uses the boost::heap and CGrid classes.
* It contains an augmented grid to compute the evolution of the boundary
*
* (See also CFMM.hpp)
*
* ==============================================================================
*/
#include "CFMM.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <iostream>
#include <cassert>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CGrid.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace memory;

/** ============================================================================
  Eikonal Constructor (RepairToYou)
==============================================================================*/
CFMM::CFMM(const shared_ptr<CGrid> aPrimary,
           const shared_ptr<array2D_t<bool>> aDomain,
           const shared_ptr<array2D_t<bool>> aBoundarySet) {
  /** Initialize properties and allocate arrays. */
  fPrimary = aPrimary;
  fDomain = aDomain;
  fBoundarySet = aBoundarySet;
  fMode = "Eikonal";

  /** In this scenario the value function at the boundary set is always equal to
  *   zero
  */
  const int nx = aPrimary->getGridSizeX();
  const int ny = aPrimary->getGridSizeY();
  fBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      (*fBoundaryValue)[i][j] = 0;
    }
  }
}

/** ============================================================================
  Eikonal Constructor (YouToRepair)
==============================================================================*/
CFMM::CFMM(const shared_ptr<CGrid> aPrimary,
           const shared_ptr<array2D_t<bool>> aDomain,
           const shared_ptr<array2D_t<bool>> aBoundarySet,
           const shared_ptr<array2D_t<double>> aBoundaryValue) {
  /** Initialize properties and allocate arrays. */
  fPrimary = aPrimary;
  fDomain = aDomain;
  fBoundarySet = aBoundarySet;
  fBoundaryValue = aBoundaryValue;
  fMode = "Eikonal";
}

/** ============================================================================
  Randomly Terminated Constructor
==============================================================================*/
CFMM::CFMM(const shared_ptr<CGrid> aPrimary,
           const shared_ptr<array2D_t<bool>> aDomain,
           const shared_ptr<array2D_t<bool>> aBoundarySet,
           const shared_ptr<array2D_t<double>> aBoundaryValue,
           const shared_ptr<array2D_t<double>> aTerminalCost) {
  /** Initialize properties and allocate arrays. */
  fPrimary = aPrimary;
  fDomain = aDomain;
  fBoundarySet = aBoundarySet;
  fBoundaryValue = aBoundaryValue;
  fTerminalCost = aTerminalCost;
  fMode = "RandomlyTerminated";
}


/** ============================================================================
  Main compute functions
==============================================================================*/
void CFMM::march() {
  /** Initialize heap */
  CFMMHeap_t CFMMHeap;
  initialize_cfmm(CFMMHeap);

  /** Perform Djikstra's algorithm to update all grid points */
  while (!CFMMHeap.empty()) {
    CHeapGP current_GP = CFMMHeap.top();
    CFMMHeap.pop();
    (*fStatus)[current_GP.fI][current_GP.fJ] = ACCEPTED;
    update_neighbors(CFMMHeap, current_GP.fI, current_GP.fJ);
  }
}

void CFMM::initialize_cfmm(CFMMHeap_t& aCFMMHeap) {
  const int nx = fPrimary->getGridSizeX();
  const int ny = fPrimary->getGridSizeY();
  fStatus = make_shared<array2D_t<status_t>>(allocateArray2D<status_t>(nx,ny));
  fHeapPointers = make_shared<array2D_t<handle_t>>(allocateArray2D<handle_t>(nx,ny));

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (!in_domain(i,j)) {
        /**
        * All points outside the domain are set to ACCEPTED and have infinite cost.
        */
        fPrimary->setValue(i, j, INF);
        (*fStatus)[i][j] = ACCEPTED;
      } else {
        /** If speed is 0, then value has to be infinity there. */
        const double f = fPrimary->getSpeed(i,j);
        if (pow(f,2) == 0) {
          fPrimary->setValue(i,j,INF);
          (*fStatus)[i][j] = ACCEPTED;
        } else {
          if (fMode == "Eikonal") {
            /** Initialize gridpoints at the goal set */
            if ((*fBoundarySet)[i][j] == true) {
              const double tempVal = (*fBoundaryValue)[i][j];
              fPrimary->setValue(i, j, tempVal);
              (*fHeapPointers)[i][j] = aCFMMHeap.push(CHeapGP(i, j, tempVal));
              (*fStatus)[i][j] = CONSIDERED;
            } else {
              fPrimary->setValue(i, j, LARGE_NUMBER);
              (*fStatus)[i][j] = FAR;
            }
          } else if (fMode == "RandomlyTerminated") {
            /** Initialize gridpoints to be zero at the goal points and the
            *   repair cost at other points
            */
            if ((*fBoundarySet)[i][j] == true) {
              const double tempVal = (*fBoundaryValue)[i][j];
              fPrimary->setValue(i, j, tempVal);

              (*fHeapPointers)[i][j] = aCFMMHeap.push(CHeapGP(i, j, tempVal));
              (*fStatus)[i][j] = CONSIDERED;
            } else {
              /** Initialize grid with the terminal cost at local minima*/
              if ((fPrimary->getBreakdownRate(i,j) > 0) && (is_local_minimum(i,j))) {
                const double K = fPrimary->getCost(i,j);
                const double BDR = fPrimary->getBreakdownRate(i,j);
                const double tempVal = (*fTerminalCost)[i][j] + K/BDR;
                fPrimary->setValue(i, j, tempVal);

                CHeapGP gp =  CHeapGP(i, j, tempVal);
                (*fHeapPointers)[i][j] = aCFMMHeap.push(gp);
                (*fStatus)[i][j] = CONSIDERED;
              } else {
                fPrimary->setValue(i, j, LARGE_NUMBER);
                (*fStatus)[i][j] = FAR;
              }

            }
          } else {
            cout << "Invalid mode selected" << endl;
            assert(false);
          }
        }
      }
    }
  }
}


/** ============================================================================
  Helper functions
==============================================================================*/
void CFMM::update_neighbors(CFMMHeap_t& aCFMMHeap, const int aCurrent_i,
                           const int aCurrent_j) {
  /* Five-point Stencil */
  constexpr int stencil[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};

  /* Iterate over neighbors. */
  for (int k = 0; k < 4; ++ k) {
    const int nbr_i = aCurrent_i + stencil[k][0];
    const int nbr_j = aCurrent_j + stencil[k][1];

    if (in_domain(nbr_i,nbr_j)) {
      /* If accepted nothing to do, otherwise, try to update */
      if ((*fStatus)[nbr_i][nbr_j] != ACCEPTED) {
        const bool gp_was_updated = update_gp(nbr_i, nbr_j);

        /* If it was already considered and was updated then update heap */
        if ((*fStatus)[nbr_i][nbr_j] == CONSIDERED) {
          if (gp_was_updated) {
            CHeapGP gp = CHeapGP(nbr_i, nbr_j, fPrimary->getValue(nbr_i,nbr_j));
            aCFMMHeap.update((*fHeapPointers)[nbr_i][nbr_j], gp);
          }
        } else {
          /* Else add to heap */
          (*fStatus)[nbr_i][nbr_j] = CONSIDERED;
          CHeapGP gp = CHeapGP(nbr_i, nbr_j, fPrimary->getValue(nbr_i,nbr_j));
          (*fHeapPointers)[nbr_i][nbr_j] = aCFMMHeap.push(gp);
        }
      }
    }
  }
}

int CFMM::smaller_h_neighbor(const int aI, const int aJ) const {
  int bestI;
  if (aI == 0) {
    bestI = 1;
  } else if (aI == fPrimary->getGridSizeX() - 1) {
    bestI = -1;
  } else if (fPrimary->getValue(aI-1, aJ) < fPrimary->getValue(aI+1, aJ)) {
    bestI = -1;
  } else {
    bestI = 1;
  }

  return bestI;
}

int CFMM::smaller_v_neighbor(const int aI, const int aJ) const {
  int bestJ;
  if (aJ == 0) {
    bestJ = 1;
  } else if (aJ == fPrimary->getGridSizeY() - 1) {
    bestJ =  -1;
  } else if (fPrimary->getValue(aI, aJ-1) < fPrimary->getValue(aI, aJ+1)) {
    bestJ = -1;
  } else {
    bestJ = 1;
  }

  return bestJ;
}

bool CFMM::is_local_minimum(const int aI, const int aJ) const {
  /** Access values of Terminal cost function */
  const double q = (*fTerminalCost)[aI][aJ]
    + (1/fPrimary->getBreakdownRate(aI,aJ)*fPrimary->getCost(aI,aJ));

  /** Get values of neighbors of q*/
  double q_1 = INF;
  double q_2 = INF;
  double q_3 = INF;
  double q_4 = INF;
  if (in_domain(aI-1,aJ)) {
    if (fPrimary->getBreakdownRate(aI-1, aJ) > 0) {
      q_1 = (*fTerminalCost)[aI-1][aJ]
        + (1/fPrimary->getBreakdownRate(aI-1,aJ)*fPrimary->getCost(aI-1,aJ));
    }
  }
  if (in_domain(aI+1,aJ)) {
    if (fPrimary->getBreakdownRate(aI-1, aJ) > 0) {
      q_2 = (*fTerminalCost)[aI+1][aJ]
        + (1/fPrimary->getBreakdownRate(aI+1,aJ)*fPrimary->getCost(aI+1,aJ));
    }
  }
  if (in_domain(aI,aJ-1)) {
    if (fPrimary->getBreakdownRate(aI-1, aJ) > 0) {
      q_3 = (*fTerminalCost)[aI][aJ-1]
        + (1/fPrimary->getBreakdownRate(aI,aJ-1)*fPrimary->getCost(aI,aJ-1));
    }
  }
  if (in_domain(aI,aJ+1)) {
    if (fPrimary->getBreakdownRate(aI-1, aJ) > 0) {
      q_4 = (*fTerminalCost)[aI][aJ+1]
        + (1/fPrimary->getBreakdownRate(aI,aJ+1)*fPrimary->getCost(aI,aJ+1));
    }
  }

  return ((q <= q_1) && (q <= q_2) && (q <= q_3) && (q <= q_4));
}


/** ============================================================================
  Update functions
==============================================================================*/
double CFMM::compute_from_two_neighbors(const int aI, const int aJ,
                                      const int aVerticalNeighbor,
                                      const int aHorizontalNeighbor) {
  /** Allocate return value */
  double new_value;

  /** Get stepsizes */
  const double dx = fPrimary->getDx();
  const double dy = fPrimary->getDy();

  /** Get shorter names for local variables */
  const double f = fPrimary->getSpeed(aI,aJ);
  const double K = fPrimary->getCost(aI,aJ);

  const double u_h = fPrimary->getValue(aI + aHorizontalNeighbor, aJ);
  const double u_v = fPrimary->getValue(aI, aJ + aVerticalNeighbor);
  double R, phi;
  if (fMode == "RandomlyTerminated") {
    R = (*fTerminalCost)[aI][aJ];
    phi = fPrimary->getBreakdownRate(aI, aJ);
  }

  /** Variables containing possible updates to primary variable */
  /** Signature is {quad1, quad2, horizontal, vertical} */
  std::vector<double> updates = {LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER};
  std::vector<bool> upwinding = {false, false, false, false};

  /** Compute discriminant */
  double A, B, C;
  if (fMode == "Eikonal") {
    A = pow(f/dx,2) + pow(f/dy,2);
    B = -2*u_h*pow(f/dx,2) - 2*u_v*pow(f/dy,2);
    C = pow(f*u_h/dx,2) + pow(f*u_v/dy,2) - pow(K,2);
  } else if (fMode == "RandomlyTerminated") {
    A = pow(f/dx,2) + pow(f/dy,2) - pow(phi,2);
    B = -2*u_h*pow(f/dx,2) - 2*u_v*pow(f/dy,2) + 2*phi*(phi*R + K);
    C = pow(f*u_h/dx,2) + pow(f*u_v/dy,2) - pow(phi*R + K,2);
  } else {
    std::cout << "Invalid fMode string" << std::endl;
    assert(false);
  }
  const double discriminant = B*B - 4*A*C;

  if (discriminant >= 0) {
    /** Compute two-sided update */
    if (A == 0) {
      /** Degenerate case: Quadratic update is actually just a linear update */
      if (B == 0) {
        updates[0] = INF;
        updates[1] = INF;
      } else {
        updates[0] = -C / B;
        updates[1] = updates[0];
      }
    } else {
      /** Quadratic formula */
      updates[0] = (-B + sqrt(discriminant)) / (2*A);
      updates[1] = (-B - sqrt(discriminant)) / (2*A);
    }

    /** Check upwinding condition for both two-sided updates */
    if (updates[0] >= max(u_h, u_v)) {
      upwinding[0] = true;
    }
    if (updates[1] >= max(u_h, u_v)) {
      upwinding[1] = true;
    }
  }

  /** Compute one sided updates. */
  if ((fMode == "Eikonal")) {
    updates[2] = u_h + dx * K / f;
    updates[3] = u_v + dy * K / f;
  } else if (fMode == "RandomlyTerminated") {
    updates[2] = (f*u_h + K*dx + phi*R*dx)/(f + phi*dx);
    updates[3] = (f*u_v + K*dy + phi*R*dy)/(f + phi*dy);
  }

  /** Check upwinding condition for one-sided updates. */
  upwinding[2] = (updates[2] >= u_h);
  upwinding[3] = (updates[3] >= u_v);

  /** Find smallest upwind update */
  int min_upwind = -1;
  double min_update = INF;
  for (int i = 0; i < 4; i++) {
    if (upwinding[i] && (updates[i] < min_update)) {
      min_upwind = i;
      min_update = updates[i];
    }
  }
  if (min_upwind == -1) {
    /** No upwind updates */
    assert(false);
    new_value = LARGE_NUMBER;
  }

  /** Use smallest upwind update */
  new_value = updates[min_upwind];
  return new_value;
}

/*-----------------------------------------------------------------------------/
/-- Update gridpoint
/-----------------------------------------------------------------------------*/
bool CFMM::update_gp(const int aI, const int aJ) {
  /** Compute smaller vertical and horizontal neighbors. */
  const int hNeighbor = smaller_h_neighbor(aI,aJ);
  const int vNeighbor = smaller_v_neighbor(aI,aJ);

  /** Compute update based on smaller horizontal and vertical neighbors */
  const double new_value = compute_from_two_neighbors(aI, aJ, vNeighbor, hNeighbor);
  if (new_value < fPrimary->getValue(aI,aJ)) {
    fPrimary->setValue(aI, aJ, new_value);
    return true;
  } else {
    return false;
  }
}
