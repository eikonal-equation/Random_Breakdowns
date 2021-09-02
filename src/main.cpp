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
* File: main.cpp
*
* Author: Marissa Gee
*
* Description: This file contains the main functions to be executed from the
* command line.
*
* ==============================================================================
*/
/** ------ Libraries ---------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "SimpleFunctions.hpp"
#include "CRandomBreakdowns.hpp"
#include "WriteToFile.hpp"
#include "CSVReader.hpp"


/** ------ Namespaces --------------------------------------------------------*/
using namespace std::placeholders;

/*==============================================================================
  Environment 1: Constant
    - Square Domain
    - Speed, running cost, breakdown rate assumed constant
    - Assumes f_2 = 0.2*f_1, f_R = 0.1*f_1
    - Computes K_i = 0.5*(FuelCost + TimeValue*TimeCost)
    - Assumes lambda_1 = 0.1*phi, lambda_2 = 0.3*phi if total breakdowns are
          included
    - Assumes R_D is ten times higher for total breakdowns
==============================================================================*/
CRandomBreakdowns Environment_1(const int aNx, const int aNy,
                                const std::string aMode, const double aRho,
                                double aDepotList[][2], const int aNumDepots,
                                double aGoalList[][2], const int aNumGoals,
                                const bool aIncludeTotalBreakdowns = true,
                                const double aSpeed = 1.0,
                                const double aBreakdownRate = 5.0,
                                const double aFuelCost = 1.0,
                                const double aTimeCost = 1.0,
                                const double aTimeValue = 1.0,
                                std::shared_ptr<std::vector<double>> aRepairCostArray = nullptr,
                                std::shared_ptr<std::vector<double>> aRepairTimeArray = nullptr,
                                const double aFieldCost = 1.0) {
  /** Set grid parameters */
  const double dx = 1.0 / (aNx - 1);
  const double dy = 1.0 / (aNy - 1);
  const int max_iterations = 200;

  /** Create level set function for border */
  std::function<double(double,double)> border_function = &square_border;

  /** Construct environment parameter functions */
  std::function<double(double,double)> working_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, aSpeed);
  std::function<double(double,double)> damaged_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.2*aSpeed);
  std::function<double(double,double)> repair_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.1*aSpeed);

  std::function<double(double,double)> working_running_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.5*(aFuelCost + aTimeValue*aTimeCost));
  std::function<double(double,double)> damaged_running_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.5*(aFuelCost + aTimeValue*aTimeCost));
  std::function<double(double,double)> repair_running_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.5*(aFuelCost + aTimeValue*aTimeCost));

  std::function<double(double,double)> working_partial_breakdown_rate_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, aBreakdownRate);
  std::function<double(double,double)> working_total_breakdown_rate_f;
  std::function<double(double,double)> damaged_total_breakdown_rate_f;
  if (aIncludeTotalBreakdowns) {
    working_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.1*aBreakdownRate);
    damaged_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.3*aBreakdownRate);
  } else {
    working_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.0);
    damaged_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.0);
  }
  std::function<double(double,double)> field_repair_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, aFieldCost);
  std::shared_ptr<std::vector<double>> partial_repair_cost_a = std::make_shared<std::vector<double>>();
  std::shared_ptr<std::vector<double>> total_repair_cost_a = std::make_shared<std::vector<double>>();

  for (int i; i < aNumDepots; ++i) {
    double partial_repair_cost;
    double partial_repair_time;
    if (aRepairCostArray == nullptr) {
      partial_repair_cost = 0;
    } else {
      partial_repair_cost = (*aRepairCostArray)[i];
    }
    if (aRepairTimeArray == nullptr) {
      partial_repair_time = 0;
    } else {
      partial_repair_time = aTimeValue*(*aRepairTimeArray)[i];
    }
    partial_repair_cost_a->push_back(partial_repair_cost + partial_repair_time);
    total_repair_cost_a->push_back(10*(partial_repair_cost + partial_repair_time));
  }

  CRandomBreakdowns solver;
  if (aMode == "RepairToYou") {
    solver = CRandomBreakdowns(aNx, aNy, dx, dy, aDepotList, aNumDepots,
                               aGoalList, aNumGoals, border_function,
                               working_speed_f, working_running_cost_f,
                               working_total_breakdown_rate_f, repair_speed_f,
                               repair_running_cost_f, total_repair_cost_a,
                               field_repair_cost_f);
  } else {
    solver = CRandomBreakdowns(aNx, aNy, dx, dy, aMode, max_iterations, aRho,
                               aDepotList, aNumDepots, aGoalList, aNumGoals,
                               border_function, working_speed_f,
                               working_running_cost_f,
                               working_partial_breakdown_rate_f, damaged_speed_f,
                               damaged_running_cost_f, partial_repair_cost_a,
                               working_total_breakdown_rate_f,
                               damaged_total_breakdown_rate_f, repair_speed_f,
                               repair_running_cost_f, total_repair_cost_a,
                               field_repair_cost_f);
  }
  return solver;
}

/*==============================================================================
  Environment 2: Inhomogenous phi
    - Square Domain
    - Speed and running cost assumed constant
    - Assumes f_1 = 1, f_2 = 0.2, f_R = 0.1
    - Assumes K_1 = K_2 = K_R = 1
    - Assumes lambda_1 = 1, lambda_2 = 3 if total breakdowns are included
    - Assumes R_D is ten times higher for total breakdowns
==============================================================================*/
CRandomBreakdowns Environment_2(const int aNx, const int aNy,
                                const std::string aMode, const double aRho,
                                double aDepotList[][2], const int aNumDepots,
                                double aGoalList[][2], const int aNumGoals,
                                const bool aIncludeTotalBreakdowns = true,
                                const double aBreakdownRate = 5.0,
                                const double aTimeValue = 1.0,
                                std::shared_ptr<std::vector<double>> aRepairCostArray = nullptr,
                                std::shared_ptr<std::vector<double>> aRepairTimeArray = nullptr,
                                const double aFieldCost = 1.0) {
  /** Set grid parameters */
  const double dx = 1.0 / (aNx - 1);
  const double dy = 1.0 / (aNy - 1);
  const int max_iterations = 200;

  /** Create level set function for border */
  std::function<double(double,double)> border_function = &square_border;

  /** Construct environment parameter functions */
  std::function<double(double,double)> working_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 1);
  std::function<double(double,double)> damaged_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.2);
  std::function<double(double,double)> repair_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.1);

  std::function<double(double,double)> working_running_cost_f = &constant2D;
  std::function<double(double,double)> damaged_running_cost_f = &constant2D;
  std::function<double(double,double)> repair_running_cost_f = &constant2D;

  std::function<double(double,double)> working_partial_breakdown_rate_f =
    std::bind(gaussian, std::placeholders::_1, std::placeholders::_2, 0.5, 0.5);
  std::function<double(double,double)> working_total_breakdown_rate_f;
  std::function<double(double,double)> damaged_total_breakdown_rate_f;
  if (aIncludeTotalBreakdowns) {
    working_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 1.0);
    damaged_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 3.0);
  } else {
    working_total_breakdown_rate_f = &zero_constant2D;
    damaged_total_breakdown_rate_f = &zero_constant2D;
  }
  std::function<double(double,double)> field_repair_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, aFieldCost);
  std::shared_ptr<std::vector<double>> partial_repair_cost_a = std::make_shared<std::vector<double>>();
  std::shared_ptr<std::vector<double>> total_repair_cost_a = std::make_shared<std::vector<double>>();

  for (int i; i < aNumDepots; ++i) {
    double partial_repair_cost;
    double partial_repair_time;
    if (aRepairCostArray == nullptr) {
      partial_repair_cost = 0;
    } else {
      partial_repair_cost = (*aRepairCostArray)[i];
    }
    if (aRepairTimeArray == nullptr) {
      partial_repair_time = 0;
    } else {
      partial_repair_time = aTimeValue*(*aRepairTimeArray)[i];
    }
    partial_repair_cost_a->push_back(partial_repair_cost + partial_repair_time);
    total_repair_cost_a->push_back(10*partial_repair_cost + partial_repair_time);
  }

  CRandomBreakdowns solver;
  solver = CRandomBreakdowns(aNx, aNy, dx, dy, aMode, max_iterations, aRho,
                             aDepotList, aNumDepots, aGoalList, aNumGoals,
                             border_function, working_speed_f,
                             working_running_cost_f,
                             working_partial_breakdown_rate_f, damaged_speed_f,
                             damaged_running_cost_f, partial_repair_cost_a,
                             working_total_breakdown_rate_f,
                             damaged_total_breakdown_rate_f, repair_speed_f,
                             repair_running_cost_f, total_repair_cost_a,
                             field_repair_cost_f);
  return solver;
}

/*==============================================================================
  Environment 2b: Inhomogenous lambda
    - Square Domain
    - Speed and running cost assumed constant
    - Assumes f_2 = 0.2*f_1, f_R = 0.1*f_1
    - Assumes K_1 = K_2 = K_R = 1
    - Assumes phi = 5
    - Assumes R_D is ten times higher for total breakdowns
==============================================================================*/
CRandomBreakdowns Environment_2b(const int aNx, const int aNy,
                                const std::string aMode, const double aRho,
                                double aDepotList[][2], const int aNumDepots,
                                double aGoalList[][2], const int aNumGoals,
                                const bool aIncludeTotalBreakdowns = true,
                                const double aBreakdownRate = 5.0,
                                const double aTimeValue = 1.0,
                                std::shared_ptr<std::vector<double>> aRepairCostArray = nullptr,
                                std::shared_ptr<std::vector<double>> aRepairTimeArray = nullptr,
                                const double aFieldCost = 1.0) {
  /** Set grid parameters */
  const double dx = 1.0 / (aNx - 1);
  const double dy = 1.0 / (aNy - 1);
  const int max_iterations = 200;

  /** Create level set function for border */
  std::function<double(double,double)> border_function = &square_border;

  /** Construct environment parameter functions */
  /** Construct environment parameter functions */
  std::function<double(double,double)> working_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 1);
  std::function<double(double,double)> damaged_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.2);
  std::function<double(double,double)> repair_speed_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.1);

  std::function<double(double,double)> working_running_cost_f = &constant2D;
  std::function<double(double,double)> damaged_running_cost_f = &constant2D;
  std::function<double(double,double)> repair_running_cost_f = &constant2D;

  std::function<double(double,double)> working_partial_breakdown_rate_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 5.0);
  std::function<double(double,double)> working_total_breakdown_rate_f;
  std::function<double(double,double)> damaged_total_breakdown_rate_f;
  if (aIncludeTotalBreakdowns) {
    working_total_breakdown_rate_f =
      std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, 0.0);
    damaged_total_breakdown_rate_f =
      std::bind(gaussian, std::placeholders::_1, std::placeholders::_2, 0.5, 0.5);
  } else {
    working_total_breakdown_rate_f = &zero_constant2D;
    damaged_total_breakdown_rate_f = &zero_constant2D;
  }
  std::function<double(double,double)> field_repair_cost_f =
    std::bind(var_constant2D, std::placeholders::_1, std::placeholders::_2, aFieldCost);
  std::shared_ptr<std::vector<double>> partial_repair_cost_a = std::make_shared<std::vector<double>>();
  std::shared_ptr<std::vector<double>> total_repair_cost_a = std::make_shared<std::vector<double>>();

  for (int i; i < aNumDepots; ++i) {
    double partial_repair_cost;
    double partial_repair_time;
    if (aRepairCostArray == nullptr) {
      partial_repair_cost = 0;
    } else {
      partial_repair_cost = (*aRepairCostArray)[i];
    }
    if (aRepairTimeArray == nullptr) {
      partial_repair_time = 0;
    } else {
      partial_repair_time = aTimeValue*(*aRepairTimeArray)[i];
    }
    partial_repair_cost_a->push_back(partial_repair_cost + partial_repair_time);
    total_repair_cost_a->push_back(10*partial_repair_cost + partial_repair_time);
  }

  CRandomBreakdowns solver;
  solver = CRandomBreakdowns(aNx, aNy, dx, dy, aMode, max_iterations, aRho,
                             aDepotList, aNumDepots, aGoalList, aNumGoals,
                             border_function, working_speed_f,
                             working_running_cost_f,
                             working_partial_breakdown_rate_f, damaged_speed_f,
                             damaged_running_cost_f, partial_repair_cost_a,
                             working_total_breakdown_rate_f,
                             damaged_total_breakdown_rate_f, repair_speed_f,
                             repair_running_cost_f, total_repair_cost_a,
                             field_repair_cost_f);
  return solver;
}

/*==============================================================================
  Environment 3: Real-world terrain based
    - Square Domain
    - Takes filenames for quantities computed from realworld data
    - Assumes lambda_1 = lambda_2 = 0
    - Assumes K_1 = K_2 = 1
==============================================================================*/
CRandomBreakdowns Environment_3(const int aNx, const int aNy,
                                const std::string aMode, const double aRho,
                                double aDepotList[][2], const int aNumDepots,
                                double aGoalList[][2], const int aNumGoals,
                                const std::string aWorkingSpeedFilename,
                                const std::string adamagedSpeedFilename,
                                const std::string aBreakdownRateFilename,
                                const bool aIncludeTotalBreakdowns = false) {
  /** Set grid parameters */
  const double dx = 1.0 / (aNx - 1);
  const double dy = 1.0 / (aNy - 1);
  const int max_iterations = 400;

  /** Create level set function for border */
  std::function<double(double,double)> border_function = &square_border;

  /** Construct environment parameter functions */
  std::shared_ptr<memory::array2D_t<double>> working_speed_array = read_csv(aWorkingSpeedFilename);
  std::shared_ptr<memory::array2D_t<double>> damaged_speed_array = read_csv(adamagedSpeedFilename);
  std::shared_ptr<memory::array2D_t<double>> breakdown_rate_array = read_csv(aBreakdownRateFilename);

  std::shared_ptr<memory::array2D_t<double>> running_cost_array =
    std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(aNx, aNy));
  for (int i = 0; i < aNx; i++) {
    for (int j = 0; j < aNy; j++) {
      (*running_cost_array)[i][j] = 1.0;
    }
  }

  CRandomBreakdowns solver;
  solver = CRandomBreakdowns(aNx, aNy, dx, dy, aMode, max_iterations, aRho,
                             aDepotList, aNumDepots, aGoalList, aNumGoals,
                             border_function, working_speed_array,
                             running_cost_array, breakdown_rate_array,
                             damaged_speed_array, running_cost_array);
  return solver;
}

int main (int argc, char* argv[]) {
  /** Start timer */
  // auto t1 = std::chrono::high_resolution_clock::now();

  /** Select example based on command line argument */
  std::shared_ptr<CRandomBreakdowns> rand_breakdown_solver;
  std::string filename;
  if (argc < 2) {
    std::cout << "No argument given, terminating" << std::endl;
    return 0;
  }

  std::string arg1 = std::string(argv[1]);

  /** Set default reporting values, that will be changed if needed */
  const double rho = 0.9;

  if (arg1 == "Example1") {
    /** Create arrays for depots and goal points */
    const int depots_n = 1;
    const int goals_n = 1;
    double depots_list[depots_n][2] = {{0.5, 0.5}};
    double goals_list[goals_n][2] = {{0.5, 0.5}};

    /** Set repair cost and time for each depot */
    std::shared_ptr<std::vector<double>> repair_cost_a = std::make_shared<std::vector<double>>();
    std::shared_ptr<std::vector<double>> repair_time_a = std::make_shared<std::vector<double>>();
    repair_cost_a->push_back(0.0);
    repair_time_a->push_back(0.0);

    /** Set environment parameters */
    const double speed = 1.0;
    const double fuel_cost = 1.0;
    const double time_cost = 1.0;
    const double time_value = 1.0;
    const double breakdown_rate = 5.0;

    const int M = 11;

    std::vector<double> errors;
    std::vector<int> mValues;
    std::vector<int> numIterations;
    std::vector<int> numEvaluations;
    const std::string model = "TwoBreakdownTypes";
    const std::string mode = "VP";

    for (int m = 3; m <= M; m++) {
      const int n = pow(2,m) + 1;
      filename = "Example1_" + model + "_" + mode + "_" + std::to_string(m);

      rand_breakdown_solver = std::make_shared<CRandomBreakdowns>(
                                Environment_1(n, n, model+mode, rho, depots_list, depots_n,
                                              goals_list, goals_n, true,
                                              speed, breakdown_rate,
                                              fuel_cost, time_cost, time_value,
                                              repair_cost_a, repair_time_a)
                              );
      rand_breakdown_solver->computeExpectedCost(filename);

      errors.push_back(rand_breakdown_solver->getError());
      numIterations.push_back(rand_breakdown_solver->getNumIterations());
      numEvaluations.push_back(rand_breakdown_solver->getNumEvaluations());
      mValues.push_back(m);
    }
    filename = "Example1_" + model + "_" + mode;
    io::writeVectorToFile<double>(filename + "_Errors", errors);
    io::writeVectorToFile<int>(filename + "_NumIter", numIterations);
    io::writeVectorToFile<int>(filename + "_NumEval", numEvaluations);
    io::writeVectorToFile<int>(filename + "_mValues", mValues);
  } else if (arg1 == "Example2") {
    const int n = 501;

    /** Create arrays for depots and goal points */
    const int depots_n = 1;
    const int goals_n = 1;
    double depots_list[depots_n][2] = {{0.9, 0.5}};
    double goals_list[goals_n][2] = {{0.9, 0.5}};

    /** Run without total breakdowns */
    std::string model = "OnlyPartialBreakdowns";
    std::string mode = "VP";
    rand_breakdown_solver = std::make_shared<CRandomBreakdowns>(
                              Environment_2(n, n, model+mode, rho, depots_list, depots_n,
                                            goals_list, goals_n, true)
                            );
    filename = "Example2_" + model + "_" + mode;
    rand_breakdown_solver->computeExpectedCost(filename);

    /** Run with total breakdowns */
    model = "TwoBreakdownTypes";
    mode = "VP";
    rand_breakdown_solver = std::make_shared<CRandomBreakdowns>(
                              Environment_2(n, n, model+mode, rho, depots_list, depots_n,
                                            goals_list, goals_n, true)
                            );
    filename = "Example2_" + model + "_" + mode;
    rand_breakdown_solver->computeExpectedCost(filename);
  } else if (arg1 == "Example3") {
    const int n = 501;

    /** Create arrays for depots and goal points */
    const int depots_n = 3;
    const int goals_n = 1;
    double depots_list[depots_n][2] = {{0.3, 0.1}, {0.35, 0.55}, {0.8, 0.6}};
    double goals_list[goals_n][2] = {{0.9, 0.9}};

    /** Set repair cost and time for each depot */
    std::shared_ptr<std::vector<double>> repair_cost_a = std::make_shared<std::vector<double>>();
    std::shared_ptr<std::vector<double>> repair_time_a = std::make_shared<std::vector<double>>();
    repair_cost_a->push_back(1.0);
    repair_cost_a->push_back(1.0);
    repair_cost_a->push_back(1.0);
    repair_time_a->push_back(1.0);
    repair_time_a->push_back(1.0);
    repair_time_a->push_back(1.0);

    /** Set environment parameters */
    const double speed = 1.0;
    const double fuel_cost = 1.0;
    const double time_cost = 1.0;
    const double time_value = 1.0;

    /** Run mode without total breakdowns */
    const std::string model = "OnlyPartialBreakdowns";
    const std::string mode = "VP";
    const int b = 4;
    int bdrs[b] = {0, 1, 3, 5};
    for (int i = 0; i < b; ++i) {
      rand_breakdown_solver = std::make_shared<CRandomBreakdowns>(
                                Environment_1(n, n, model+mode, rho, depots_list, depots_n,
                                              goals_list, goals_n, false, speed, bdrs[i],
                                              fuel_cost, time_cost, time_value)
                              );
      filename = "Example3_" + model + "_" + mode + "_" + std::to_string(bdrs[i]);
      rand_breakdown_solver->computeExpectedCost(filename);
    }
  } else if (arg1 == "Example4") {
    const int n = 501;

    /** Create arrays for depots and goal points */
    const int depots_n = 1;
    const int goals_n = 1;
    double depots_list[depots_n][2] = {{0.1, 0.1}};
    double goals_list[goals_n][2] = {{0.1, 0.1}};

    /** Set environment parameters */
    const std::string working_speed_filename = "data/Ex4_f_1.csv";
    const std::string damaged_speed_filename = "data/Ex4_f_2.csv";
    const std::string breakdown_rate_filename = "data/Ex4_phi.csv";

    /** Run established mode */
    const std::string model = "OnlyPartialBreakdowns";
    const std::string mode = "VP";
    rand_breakdown_solver = std::make_shared<CRandomBreakdowns>(
                              Environment_3(n, n, model+mode, rho, depots_list,
                                            depots_n, goals_list, goals_n,
                                            working_speed_filename,
                                            damaged_speed_filename,
                                            breakdown_rate_filename)
                            );
    filename = "Example4_" + model + "_" + mode;
    rand_breakdown_solver->computeExpectedCost(filename);
  } else {
    std::cout << "Invalid testing input" << std::endl;
    assert(false);
  }

  // /** Stop timer */
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << "Complete, took "
  //           << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
  //           << " milliseconds" << std::endl;

  return 0;
}
