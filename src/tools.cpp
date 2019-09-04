#include "tools.h"
#include <iostream>
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // Parameter checking for invalid inputs
  if (!estimations.size() or estimations.size() != ground_truth.size()) {
      std::cout << "Invalid inputs" << std::endl;
      return rmse;
  }
  
  // Accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    // ... your code here
    for (unsigned int j = 0; j < 4; ++j) {
        rmse(j) += pow(estimations[i](j) - ground_truth[i](j), 2);
    }
  }

  // TODO: calculate the mean
  for (unsigned int j = 0; j < 4; ++j) {
    rmse(j) /= estimations.size();
    rmse(j) = sqrt(rmse(j));
  }

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */

  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  if (px == 0.0 and py == 0.0) {
      std::cout << "Division by zero!" << std::endl;
      Hj << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; // initilize to zero
      return Hj;
  }
  
  float px2 = px * px;
  float py2 = py * py;
  float sqrt2 = std::sqrt(px2 + py2);
  float p32 = (px2 + py2) * sqrt2;
  
  // compute the Jacobian matrix
  Hj(0,0) = px / sqrt2;
  Hj(0,1) = py / sqrt2;
  Hj(0,2) = 0;
  Hj(0,3) = 0;
  Hj(1,0) = - py / (px2 + py2);
  Hj(1,1) = px / (px2 + py2);
  Hj(1,2) = 0;
  Hj(1,3) = 0;
  Hj(2,0) = py * (vx*py - vy*px) / p32;
  Hj(2,1) = px * (vy*px - vx*py) / p32;
  Hj(2,2) = px / sqrt2;
  Hj(2,3) = py / sqrt2;
  
  return Hj;
}
