#include "kalman_filter.h"
#include "tools.h"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

static float getPI() {
  static float PI = std::atan(1.0) * 4.0;
  return PI;
}
// Function to normalize phi within -pi <= phi <= pi
static float NormalizePhi(float phi) {
  static float PI = getPI();
  static float PI2 = PI * 2.0;
  while (phi < -PI or phi > PI) {
    phi = (phi < 0)? phi + PI2: phi - PI2;
  }
  return phi;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  // Normalize phi of z using new vector nz
  VectorXd nz(3);
  nz << z[0], NormalizePhi(z[1]), z[2];
  
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  // Converting the cartesian state coords to the polar measurement ones
  float rho = std::sqrt(px * px + py * py);
  float phi = 0;
  if (px < 0.001 and px > -0.001) {
    phi = getPI() / 2.0;
    if (py < 0) phi *= -1;
  }
  else {
    phi = atan2(py, px);
    phi = NormalizePhi(phi);
  }
  float rhodot = (px*vx + py*vy) / rho;
  VectorXd hx(3);
  hx << rho, phi, rhodot;
  VectorXd y = nz - hx;
  // Normalize phi in y
  y[1] = NormalizePhi(y[1]);
  MatrixXd Hjt = H_.transpose();
  MatrixXd S = H_ * P_ * Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
