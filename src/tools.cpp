#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if(estimations.size() != ground_truth.size() 
	|| estimations.size() == 0){
	
	cout << "Invalid estimation or ground_truth data" << endl;
	return rmse;
  }

  //Calculating differences between corresponding vectors squared and adding results to rmse vector
  for(unsigned int i=0; i < estimations.size(); ++i){

      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      rmse += residual;
  }

  //calculating the mean from the total differences squared
  rmse = rmse/estimations.size();

  //calculating the squared root
  rmse = rmse.array().sqrt();

  //returning the result
  return rmse;
}
