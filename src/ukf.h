#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

	///* initially set to false, set to true in first call of ProcessMeasurement
	bool is_initialized_;

	///* if this is false, laser measurements will be ignored (except for init)
	bool use_laser_;

	///* if this is false, radar measurements will be ignored (except for init)
	bool use_radar_;

	///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	VectorXd x_;

	///* state covariance matrix
	MatrixXd P_;

	///* predicted sigma points matrix
	MatrixXd Xsig_pred_;

	// previous timestamp
	long previous_timestamp_;

	///* Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a_;

	///* Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd_;

	///* Laser measurement noise standard deviation position1 in m
	double std_laspx_;

	///* Laser measurement noise standard deviation position2 in m
	double std_laspy_;

	///* Radar measurement noise standard deviation radius in m
	double std_radr_;

	///* Radar measurement noise standard deviation angle in rad
	double std_radphi_;

	///* Radar measurement noise standard deviation radius change in m/s
	double std_radrd_ ;

	///* Weights of sigma points
	VectorXd weights_;

	///* State dimension
	int n_x_;

	///* Augmented state dimension
	int n_aug_;

	///* Sigma point spreading parameter
	double lambda_;

	///* the current NIS for radar
	double NIS_radar_;

	///* the current NIS for laser
	double NIS_laser_;

	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z_;

	/**
	 * Constructor
	 */
	UKF(bool use_laser, bool use_radar);

	/**
	 * Destructor
	 */
	virtual ~UKF();

	/**
	 * ProcessMeasurement
	 * @param meas_package The latest measurement data of either radar or laser
	 */
	void ProcessMeasurement(MeasurementPackage measurement_pack);

	/**
	 * FindFirstMeasurement
	 * @param meas_package The latest measurement data of either radar or laser
	 */
	bool FindFirstMeasurement(const MeasurementPackage &measurement_pack);

	/**
	 * Prediction Predicts sigma points, the state, and the state covariance
	 * matrix
	 * @param delta_t Time between k and k+1 in s
	 */
	void Prediction(MatrixXd &Xsig_pred, double delta_t);

	/**
	 * Updates the state and the state covariance matrix using a laser measurement
	 * @param meas_package The measurement at k+1
	 */
	void UpdateLidar(const MeasurementPackage &meas_package);

	/**
	 * Updates the state and the state covariance matrix using a radar measurement
	 * @param meas_package The measurement at k+1
	 */
	void UpdateRadar(const MeasurementPackage &meas_package, const MatrixXd &Xsig_pred);
private:
	void AugmentedSigmaPoints(MatrixXd &Xsig_aug);
	void SigmaPointPrediction(MatrixXd &Xsig_pred, double delta_t, const MatrixXd &Xsig_aug);
	void PredictMeanAndCovariance(const MatrixXd &Xsig_pred);
	void PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig, const MatrixXd &Xsig_pred);
	void UpdateRadarState(const MeasurementPackage &meas_package, const MatrixXd &Xsig_pred, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S);
	void NormalizeAngle(double &angle);
	double ComputeNIS(const VectorXd &z_pred, const MatrixXd &S, const VectorXd &z);
};

#endif /* UKF_H */
