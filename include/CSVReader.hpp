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
* File: CSVReader.hpp
*
* Authors: Marissa Gee
*
* Description: This file contains inline definitions for reading numeric data
* from a CSV file into a Boost multi-array
*
* ==============================================================================
*/

#ifndef CSVREADER_HPP
#define CSVREADER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "MemoryAllocations.hpp"

/** Function for reading CSV file into boost::multi_array.
 * @param aFilename  string the name of the CSV file to be read.
 */
std::shared_ptr<memory::array2D_t<double>> read_csv(std::string aFilename) {
  int nx = 0;
  int ny = 0;

  std::ifstream file(aFilename);
  std::queue<double> values;

  if (file.good()) {
    std::string line;
    double dataPoint;
    std::string test;

    std::getline(file, line);
    ny++;
    std::stringstream lineStream(line);

    /** Extract first row and compute nx */
    // std::getline(lineStream, test, ',');
    while (lineStream >> dataPoint) {
      values.push(dataPoint);
      nx++;

      if(lineStream.peek() == ',') {
        lineStream.ignore();
      }
    }

    /** Extract the remaining rows */
    while (std::getline(file, line)) {
      ny++;
      std::stringstream lineStream(line);
      double dataPoint;
      while (lineStream >> dataPoint) {
        values.push(dataPoint);

        if(lineStream.peek() == ',') {
          lineStream.ignore();
        }
      }
    }
  } else {
    std::cout << "Could not open specified file." << std::endl;
    assert(false);
  }

  /** Format the results as a 2D array */
  std::shared_ptr<memory::array2D_t<double>> result =
    std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(nx, ny));

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      (*result)[j][i] = values.front();
      values.pop();
    }
  }
  return result;
}

#endif
