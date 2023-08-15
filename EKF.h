/*
 * EKF.h
 *
 *  Created on: Apr 6, 2023
 *      Author: Emre Emir Fidan
 */

#ifndef INC_EKF_H_
#define INC_EKF_H_

#include "stm32l4xx_hal.h"

typedef struct {
	float P[4][4];	// Prediction error covariance matrix
	float Q[4][4];  // Process noise covariance matrix
	float R[6][6];	// Measurement noise covariance matrix

	float K[4][6];

	float x[4];
	float xp[4];

	float Pp[4][4];

	// Magnetic Vector References
	float ref_mx;
	float ref_my;
	float ref_mz;
}ekf_t;

void EKF_init(ekf_t* ekf, float ref_mx, float ref_my, float ref_mz, float N_Q, float N_P, float N_R);

void EKF_update(ekf_t* ekf, float euler[3], float ax, float ay, float az, float p, float q, float r, float mx, float my, float mz, float dt);

uint8_t inverse_matrix(float a[6][6], float a_inv[6][6]);

void EKFquaternionToEuler(float q[4], float euler[3]);


#endif /* INC_EKF_H_ */
