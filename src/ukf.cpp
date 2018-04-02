#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

float MICROSECONDS_PER_SECOND = 1000000.0;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  //Always initialized as false
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //Current state time in microseconds
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // this is probably too big
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //this one also seems too big 
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Mean Process noise vector
  mean_process_noise_ = VectorXd(2);
  mean_process_noise_(0) = 0;
  mean_process_noise_(1) = 0;

  // Standard deviation process noise vector
  std_process_noise_ = VectorXd(mean_process_noise_.size());
  std_process_noise_(0) = std_a_;
  std_process_noise_(1) = std_yawdd_;

  // State dimension is determined by size of x_
  n_x_ = x_.size();

  // Noise dimenions is determined by process_noise
  n_noise_ = mean_process_noise_.size();

  // Augmented state dimension is determined by the sum of process noise elements and x elements
  n_aug_ = n_x_ + n_noise_;
 
  // lambda is an experimental constant
  lambda_ = 3 - n_aug_;
 
  //Weights for sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_ + 1; i++) {
    double weight_ = 0.5/(lambda_ + n_aug_);
    weights_(i) = weight_;
  }

  //Sigma points matrix
  Xsig_pred_ = MatrixXd(n_aug_,n_aug_);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  MeasurementPackage::SensorType current_sensor = meas_package.sensor_type_;

  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    float x, y, velocity, phi, phi_delta;
   
    switch(current_sensor)
    { 
      case MeasurementPackage::SensorType::LASER : 
      {
        float meas_rho = meas_package.raw_measurements_[0];
        float meas_phi = meas_package.raw_measurements_[1]; //in radians probably
        float meas_rho_delta = meas_package.raw_measurements_[2];
      
        x = cos(phi)*meas_rho;
        y = sin(phi)*meas_rho;
        velocity = meas_rho_delta;
        phi = meas_phi;
        phi_delta = 0;
      
        break;
      } 
      case MeasurementPackage::SensorType::RADAR :
      { 
        float meas_x = meas_package.raw_measurements_[0];
        float meas_y = meas_package.raw_measurements_[1];
     
        x = meas_x;
        y = meas_y;
        velocity = 0;
        phi = 0;
        phi_delta = 0;
 
        break;
      }
    } 

    x_ << x, y, velocity, phi, phi_delta;

    //no need to predict or update
    is_initialized_ = true;
    return;
  }

  double elapsed_time = (meas_package.timestamp_ - time_us_)/MICROSECONDS_PER_SECOND; 
  time_us_ = meas_package.timestamp_;
 
  //Predicting state from previous state 
  Prediction(elapsed_time); 

  //Updating state based off latest measurement
  switch(current_sensor)
  { 
    case MeasurementPackage::SensorType::RADAR : {
      UpdateRadar(meas_package);
      break;
    }
    case MeasurementPackage::SensorType::LASER : {
      UpdateLidar(meas_package);
      break;
    }
  }

  return; 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  GenerateSigmaPoints(delta_t);
  
}

/**
 * Generates Sigma Points for current iteration of prediction and update
 */
void UKF::GenerateSigmaPoints(double delta_t) {
  //Generating square root of augmented version of covariance matrix
  MatrixXd A = MatrixXd(n_aug_, n_aug_);
  A.fill(0.0);
  A.topLeftCorner(n_x_,n_x_) = P_;
  for (int i = n_x_; i < n_aug_; i++) {
    A(i,i) = pow(std_process_noise_(i-n_x_),2); 
  }
  MatrixXd L = A.llt().matrixL();

  //Generating augmented state vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.head(n_x_) = x_;
  x_aug_.tail(n_noise_) = mean_process_noise_;

  //Filling out sigma points matrix
  Xsig_pred_.fill(0.0);
  Xsig_pred_.col(0) = x_aug_;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_pred_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_pred_.col(i+1+n_x_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  return; 

}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
