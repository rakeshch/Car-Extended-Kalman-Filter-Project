#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
	  * predict the state
	*/
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update_L(const VectorXd &z) {
	/**
	  * update the state by using Kalman Filter equations
	*/
	//Calculate hx_ using current state and H_ matrix initialized in previous step
	VectorXd hx_ = H_ * x_;

	//Update state and covariance matrix
	UpdateEKF(z, hx_);
}

void KalmanFilter::Update_R(const VectorXd &z) {
	/**
	* update the state by using Extended Kalman Filter equations
	*/
	Tools tools;

	//Calculate Jacobian matrix using cartesian coordinates (px,py,vx,vy) derived from (rho, phi, rhodot) equations
	H_ = MatrixXd(3, 4);
	H_ = tools.CalculateJacobian(tools.PolarToCartesian(z));

	//Skip the measurement update step if px*px + py*py is zero or close to zero
	if(x_[0]*x_[0] + x_[1]*x_[1] >=0.0001)
	{
	   //Convert Cartesian coordinates back to Polar coordinates
	   VectorXd hx_ = VectorXd(4);
	   hx_ = tools.CartesianToPolar(x_);

	   //Update state and covariance matrix 
	   UpdateEKF(z, hx_);
	}	
}

void KalmanFilter::UpdateEKF(const VectorXd &z, VectorXd &hx_) {
	VectorXd y = z - hx_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//If measurement is from Radar: 
	//Normalize second value phi in the polar coordinate vector y=(rho, phi, rhodot) to be in between -pi and pi
	//using atan2 to return values between -pi and pi
	if (y.size() == 3)
		y(1) = atan2(sin(y(1)), cos(y(1)));

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
