#include <iostream>
#include <math.h>
#include "ukf.h"
#include "tools.h"

using namespace std;

// for convenience

int main()
{
  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

        
  MeasurementPackage meas_package;
  meas_package.sensor_type_ = MeasurementPackage::LASER;
  meas_package.raw_measurements_ = VectorXd(2);
  float px = 1;
  float py = 1;
  meas_package.raw_measurements_ << px, py;
  float timestamp = 12345;
  meas_package.timestamp_ = timestamp;

  /*
  float x_gt = 2;
  float y_gt = 2;
  float vx_gt = 2;
  float vy_gt = 2;
  VectorXd gt_values(4);
  gt_values(0) = x_gt;
  gt_values(1) = y_gt; 
  gt_values(2) = vx_gt;
  gt_values(3) = vy_gt;
  ground_truth.push_back(gt_values);
  */
        
  //Call ProcessMeasurment(meas_package) for Kalman filter
  ukf.ProcessMeasurement(meas_package);    	  

  //Push the current estimated x,y positon from the Kalman filter's state vector

  VectorXd estimate(4);

  double p_x = ukf.x_(0);
  double p_y = ukf.x_(1);
  double v  = ukf.x_(2);
  double yaw = ukf.x_(3);

  double v1 = cos(yaw)*v;
  double v2 = sin(yaw)*v;

  estimate(0) = p_x;
  estimate(1) = p_y;
  estimate(2) = v1;
  estimate(3) = v2;
    	  
  estimations.push_back(estimate);

  //VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);

  //next two lines check if the values are instantiated correctly  
  if (px != p_x) {
    return 1;
  }
  if (py != p_y) {
    return 1;
  }

  return 0; 
}
