#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);
  /*
  P_ << 0.5, 0.0, 0.0, 0,0,
     0.0, 0.5, 0.0, 0,0,
     0.0, 0.0, 10.0, 0,0,
     0.0, 0.0, 0.0, 10.0,0,
     0.0,0.0,0.0,0.0,10.0;
   */

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;//0.8;//0.4;//30;15
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;//0.9;//0.05;//30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  previous_timestamp_ = 1477010443000000;


  H_laser_ = MatrixXd(2, 5);
  H_laser_.fill(0.0);
  H_laser_ << 1.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0, 0.0;

  //measurement covariance
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0.0,
              0.0, std_laspy_*std_laspy_;

  // State dimension
  n_x_ = 5;
  // Augmented state dimension
  n_aug_ = n_x_ + 2;
  // number of augmented sigma points (state + procoess noise)
  n_sigma_ = 2 * n_aug_ + 1;
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  // the current NIS for radar
  NIS_radar_ = 0;
  // the current NIS for laser
  NIS_laser_ = 0;

  // create vector for weights
  weights_ = VectorXd(n_sigma_);
  // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < n_sigma_; i++)
    weights_(i) = (double) 0.5/(lambda_ + n_aug_);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //meas_package.raw_measurements_ << ro, theta, ro_dot;
      // px is cosine(theta)
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      /*
      * Thi s is the CTRV Model with 5 states
      */
      float px = cos(phi) * rho;
      float py = sin(phi) * rho;
      // vx
      float vx = 0;
      // vy
      float vy = 0;
      // v this should be constant
      float v = 0;//5;//0.01;
      // yaw is the orientation
      float yaw = 0;
      // yaw rate.  This should be constant
      float yawd = 0;//0.1;

      x_ << px, py, v, yaw, yawd;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //cout << meas_package.raw_measurements_ << endl;

      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      // v this should be constant
      float v = 0;//5;//0.01;
      // yaw is the orientation
      float yaw = 0;
      // yaw rate.  This should be constant
      float yawd = 0;//0.1;
      x_ << px, py, v, yaw, yawd;

    }
    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  //cout << "Prediction Radar with delta_t = " << dt << endl;
  Prediction(dt);

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
#if 1 
    UpdateLidar(meas_package);
#else    
    EKFUpdateLaser(meas_package);
#endif    
  }

}

double UKF::NormalizeAngle(double angle) {

  //angle normalization
  while (angle> M_PI) {
    angle-=2.*M_PI;
  }
  while (angle<-M_PI) {
    angle+=2.*M_PI;
  }

  return angle;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  //create augmented state vector per prediction iteration initialized by x_
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0; //accel has 0 mean noise.  This is always the case.
  x_aug(6) = 0.0; //yawdd has 0 mean noise.  This is always the case.

  //create augmented sigma point matrix per prediction iteration
  // n_aug_ x n_sigma_ matrix size 
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  Xsig_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_, n_x_) = std_yawdd_ * std_yawdd_;

  //calculate square root of P_aug
  MatrixXd A = P_aug.llt().matrixL();

  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+ n_aug_) * A.col(i);
  }

  //predict sigma points
  for (int i = 0; i< n_sigma_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    yaw = NormalizeAngle(yaw);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    yaw_p = NormalizeAngle(yaw_p);

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    yaw_p = NormalizeAngle(yaw_p);

    //cout << "old:" << yaw << " new:" << yaw_p << " yawd:" << yawd << endl;
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

 } 
}


MatrixXd UKF::S_laser_inverse() {
  return (H_laser_ * P_ * H_laser_.transpose() + R_laser_).inverse();
}

#if 0
void UKF::EKFUpdateLaser(MeasurementPackage meas_package) {
  /**
  DONE:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd y = meas_package.raw_measurements_ - H_laser_*x_;

  MatrixXd Ht = H_laser_.transpose();
  //MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S_laser_inverse(); // S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  VectorXd y_new = meas_package.raw_measurements_ - H_laser_*x_;
  MatrixXd Si_new = S_laser_inverse(); // (H_laser_ * P_ * Ht + R_laser_).inverse();
  NIS_laser_ = y_new.transpose()*Si_new*y_new;
}
#else

void UKF::EKFUpdateLaser(MeasurementPackage meas_package) {
  /**
  DONE:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd y = meas_package.raw_measurements_ - H_laser_*x_;

  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  VectorXd y_new = meas_package.raw_measurements_ - H_laser_*x_;
  MatrixXd Si_new = (H_laser_ * P_ * Ht + R_laser_).inverse();
  NIS_laser_ = y_new.transpose()* Si_new *y_new;
}
#endif

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

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  //create vector for predicted state meanå
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  //create covariance matrix for predicted state covariance matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  Zsig.fill(0.0);
  // predicted mean in measurement space
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  //add measurement noise covariance matrix
  MatrixXd Q = MatrixXd(n_z,n_z);
  Q.fill(0.0);
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    if (p_x == 0 && p_y == 0) {
      Zsig(0,i) = 0;
      Zsig(1,i) = 0;
    } else {
      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;
    }
  }

 // calculate the predicted mean
  for (int i=1; i < n_sigma_; i++)
  {
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate the predicated covariance
  for (int i=1; i < n_sigma_; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    x_diff(3) = NormalizeAngle(x_diff(3));

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //z_diff(1) = NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  Q << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  S = S + Q;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff =  z - z_pred;

  //z_diff(1) = NormalizeAngle(z_diff(1));


  //update state mean and covariance matrix
  x_ = x_pred + K * z_diff;
  P_ = P_pred - K*S*K.transpose();
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  //cout << "NIS_laser_ = " << NIS_laser_ << endl;

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

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create vector for predicted state meanå
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  //create covariance matrix for predicted state covariance matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  Zsig.fill(0.0);
  // predicted mean in measurement space
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v; 

    // measurement model
    if (p_x == 0 && p_y == 0) {
      Zsig(0,i) = 0;
      Zsig(1,i) = 0;
      Zsig(2,i) = 0;
    } else {
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
  }

  // calculate the predicted mean
  for (int i=1; i < n_sigma_; i++)
  {
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate the predicated covariance
  for (int i=1; i < n_sigma_; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    x_diff(3) = NormalizeAngle(x_diff(3));

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff =  z - z_pred;

  z_diff(1) = NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_pred + K * z_diff;
  P_ = P_pred - K*S*K.transpose();
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  //cout << "NIS_radar_ = " << NIS_radar_ << endl;

}
