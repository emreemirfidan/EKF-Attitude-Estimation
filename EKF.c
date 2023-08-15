/*
 * EKF.c
 *
 *  Created on: Apr 6, 2023
 *      Author: Emre Emir Fidan
 */

#include "EKF.h"
#include "math.h"

void EKF_init(ekf_t* ekf, float ref_mx, float ref_my, float ref_mz, float N_Q, float N_P, float N_R) {
	// Initialization:
	// Prediction error covariance matrix
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			if (i == j) {
				ekf->P[i][j] = N_P;
			}
			else {
				ekf->P[i][j] = 0.0f;
			}
		}
	}

	// Process noise covariance matrix
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			if (i == j) {
				ekf->Q[i][j] = N_Q;
			}
			else {
				ekf->Q[i][j] = 0.0f;
			}
		}
	}

	// Measurement noise covariance matrix
	for (uint8_t i = 0; i < 6; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			if (i == j) {
				ekf->R[i][j] = N_R;
			}
			else {
				ekf->R[i][j] = 0.0f;
			}
		}
	}

	ekf->x[0] = 1;
	ekf->x[1] = 0;
	ekf->x[2] = 0;
	ekf->x[3] = 0;


	// Normalize Reference Magnetic Vector
	float M = sqrtf(ref_mx * ref_mx + ref_my * ref_my + ref_mz * ref_mz);
	ekf->ref_mx = ref_mx / M;
	ekf->ref_my = ref_my / M;
	ekf->ref_mz = ref_mz / M;
}


void EKF_update(ekf_t* ekf, float euler[3], float ax, float ay, float az, float p, float q, float r, float mx, float my, float mz, float dt) {
	// Variable Definitions
	float F[4][4];	// Jacobian matrix of F
	float H[6][4];	// Jacobian matrix of H
	float FP[4][4];
	float FPFt[4][4];
	float HPp[6][4];
	float HPpHt[6][6];
	float PpHt[4][6];
	float S_inv[6][6];
	float Hxp[6];
	float z[6];
	float zmHxp[6];
	float KzmHxp[4];
	float KH[4][4];
	float KHPp[4][4];
	float G;
	float M;
	uint8_t mat_error;


	// Normalization
	G = sqrtf(ax * ax + ay * ay + az * az);
	M = sqrtf(powf(mx, 2) + powf(my, 2) + powf(mz, 2));
	ax = ax / G;
	ay = ay / G;
	az = az / G;
	mx = mx / M;
	my = my / M;
	mz = mz / M;


	// Calculate the Jacobian matrix of F
	F[0][0] = 1;
	F[0][1] = -p * dt / 2;
	F[0][2] = -q * dt / 2;
	F[0][3] = -r * dt / 2;

	F[1][0] = p * dt / 2;
	F[1][1] = 1;
	F[1][2] = r * dt / 2;
	F[1][3] = -q * dt / 2;

	F[2][0] = q * dt / 2;
	F[2][1] = -r * dt / 2;
	F[2][2] = 1;
	F[2][3] = p * dt / 2;

	F[3][0] = r * dt / 2;
	F[3][1] = q * dt / 2;
	F[3][2] = -p * dt / 2;
	F[3][3] = 1;


	// Prediction of x
	ekf->xp[0] = F[0][0] * ekf->x[0] + F[0][1] * ekf->x[1] + F[0][2] * ekf->x[2] + F[0][3] * ekf->x[3];
	ekf->xp[1] = F[1][0] * ekf->x[0] + F[1][1] * ekf->x[1] + F[1][2] * ekf->x[2] + F[1][3] * ekf->x[3];
	ekf->xp[2] = F[2][0] * ekf->x[0] + F[2][1] * ekf->x[1] + F[2][2] * ekf->x[2] + F[2][3] * ekf->x[3];
	ekf->xp[3] = F[3][0] * ekf->x[0] + F[3][1] * ekf->x[1] + F[3][2] * ekf->x[2] + F[3][3] * ekf->x[3];


	// Prediction of P
	/* Pp = F*P*F' + Q; */
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			FP[i][j] = F[i][0] * ekf->P[0][j] + F[i][1] * ekf->P[1][j] + F[i][2] * ekf->P[2][j] + F[i][3] * ekf->P[3][j];
		}
	}

	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			FPFt[i][j] = FP[i][0] * F[j][0] + FP[i][1] * F[j][1] + FP[i][2] * F[j][2] + FP[i][3] * F[j][3];
		}
	}

	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			ekf->Pp[i][j] = FPFt[i][j] + ekf->Q[i][j];
		}
	}


	// Calculate the Jacobian matrix of h
	H[0][0] = -ekf->x[2];
	H[0][1] = ekf->x[3];
	H[0][2] = -ekf->x[0];
	H[0][3] = ekf->x[1];

	H[1][0] = ekf->x[1];
	H[1][1] = ekf->x[0];
	H[1][2] = ekf->x[3];
	H[1][3] = ekf->x[2];

	H[2][0] = ekf->x[0];
	H[2][1] = -ekf->x[1];
	H[2][2] = -ekf->x[2];
	H[2][3] = ekf->x[3];

	H[3][0] = ekf->x[0] * ekf->ref_mx + ekf->x[3] * ekf->ref_my - ekf->x[2] * ekf->ref_mz;
	H[3][1] = ekf->x[1] * ekf->ref_mx + ekf->x[2] * ekf->ref_my + ekf->x[3] * ekf->ref_mz;
	H[3][2] = -ekf->x[2] * ekf->ref_mx + ekf->x[1] * ekf->ref_my - ekf->x[0] * ekf->ref_mz;
	H[3][3] = -ekf->x[3] * ekf->ref_mx + ekf->x[0] * ekf->ref_my + ekf->x[1] * ekf->ref_mz;

	H[4][0] = -ekf->x[3] * ekf->ref_mx + ekf->x[0] * ekf->ref_my + ekf->x[1] * ekf->ref_mz;
	H[4][1] = ekf->x[2] * ekf->ref_mx - ekf->x[1] * ekf->ref_my + ekf->x[0] * ekf->ref_mz;
	H[4][2] = ekf->x[1] * ekf->ref_mx + ekf->x[2] * ekf->ref_my + ekf->x[3] * ekf->ref_mz;
	H[4][3] = -ekf->x[0] * ekf->ref_mx - ekf->x[3] * ekf->ref_my + ekf->x[2] * ekf->ref_mz;

	H[5][0] = ekf->x[2] * ekf->ref_mx - ekf->x[1] * ekf->ref_my + ekf->x[0] * ekf->ref_mz;
	H[5][1] = ekf->x[3] * ekf->ref_mx - ekf->x[0] * ekf->ref_my - ekf->x[1] * ekf->ref_mz;
	H[5][2] = ekf->x[0] * ekf->ref_mx + ekf->x[3] * ekf->ref_my - ekf->x[2] * ekf->ref_mz;
	H[5][3] = ekf->x[1] * ekf->ref_mx + ekf->x[2] * ekf->ref_my + ekf->x[3] * ekf->ref_mz;


	// S = (H*Pp*H' + R);
	// H*Pp
	for (uint8_t i = 0; i < 6; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			HPp[i][j] = H[i][0] * ekf->Pp[0][j] + H[i][1] * ekf->Pp[1][j] + H[i][2] * ekf->Pp[2][j] + H[i][3] * ekf->Pp[3][j];
		}
	}

	// H*Pp*H'
	for (uint8_t i = 0; i < 6; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			HPpHt[i][j] = HPp[i][0] * H[j][0] + HPp[i][1] * H[j][1] + HPp[i][2] * H[j][2] + HPp[i][3] * H[j][3];
		}
	}

	// H*Pp*H' + R
	for (uint8_t i = 0; i < 6; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			HPpHt[i][j] = HPpHt[i][j] + ekf->R[i][j]; // S
		}
	}


	// K = Pp*H'*(S.inverse);
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			PpHt[i][j] = ekf->Pp[i][0] * H[j][0] + ekf->Pp[i][1] * H[j][1] + ekf->Pp[i][2] * H[j][2] + ekf->Pp[i][3] * H[j][3];
		}
	}

	mat_error = inverse_matrix(HPpHt, S_inv);
	if (!mat_error) {
		return;
	}

	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			ekf->K[i][j] = PpHt[i][0] * S_inv[0][j] + PpHt[i][1] * S_inv[1][j] + PpHt[i][2] * S_inv[2][j] + PpHt[i][3] * S_inv[3][j] + PpHt[i][4] * S_inv[4][j] + PpHt[i][5] * S_inv[5][j];
		}
	}


	// x = xp + K*(z - H*xp);
	// H*xp
	for (uint8_t i = 0; i < 6; i++) {
		Hxp[i] = ekf->xp[0] * H[i][0] + ekf->xp[1] * H[i][1] + ekf->xp[2] * H[i][2] + ekf->xp[3] * H[i][3];
	}

	z[0] = ax;
	z[1] = ay;
	z[2] = az;
	z[3] = mx;
	z[4] = my;
	z[5] = mz;

	// (z - H*xp)
	for (uint8_t i = 0; i < 6; i++) {
		zmHxp[i] = (z[i] - Hxp[i]);
	}

	// K*(z - H*xp)
	for (uint8_t i = 0; i < 4; i++) {
		KzmHxp[i] = 0;
		for (uint8_t j = 0; j < 6; j++) {
			KzmHxp[i] += ekf->K[i][j] * zmHxp[j];
		}
	}

	// xp + K*(z - H*xp)
	for (uint8_t i = 0; i < 4; i++) {
		ekf->x[i] = ekf->xp[i] + KzmHxp[i];
	}


	// P = Pp - K*H*Pp;
	// K*H
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			KH[i][j] = ekf->K[i][0] * H[0][j] + ekf->K[i][1] * H[1][j] + ekf->K[i][2] * H[2][j] + ekf->K[i][3] * H[3][j] + ekf->K[i][4] * H[4][j] + ekf->K[i][5] * H[5][j];
		}
	}

	// K*H*Pp
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			KHPp[i][j] = KH[i][0] * ekf->Pp[0][j] + KH[i][1] * ekf->Pp[1][j] + KH[i][2] * ekf->Pp[2][j] + KH[i][3] * ekf->Pp[3][j];
		}
	}

	// P = Pp - K*H*Pp
	for (uint8_t i = 0; i < 4; i++) {
		for (uint8_t j = 0; j < 4; j++) {
			ekf->P[i][j] = ekf->Pp[i][j] - KHPp[i][j];
		}
	}


	EKFquaternionToEuler(ekf->x, euler);

}

uint8_t inverse_matrix(float a[6][6], float a_inv[6][6]) {
	// Initialize matrix
	for (uint8_t i = 0; i < 6; i++) {
		for (uint8_t j = 0; j < 6; j++) {
			if (i != j)
				a_inv[i][j] = 0;
			else
				a_inv[i][j] = 1;
		}
	}

	// Applying Gauss Jordan Elimination
	float ratio = 1.0;
	for (uint8_t i = 0; i < 6; i++)
	{
		if (a[i][i] == 0.0)
		{
			return 0;
		}
		for (uint8_t j = 0; j < 6; j++)
		{
			if (i != j)
			{
				ratio = a[j][i] / a[i][i];
				for (uint8_t k = 0; k < 6; k++)
				{
					a[j][k] = a[j][k] - ratio * a[i][k];
					a_inv[j][k] = a_inv[j][k] - ratio * a_inv[i][k];
				}
			}
		}
	}
	// Row Operation to Make Principal Diagonal to 1
	for (uint8_t i = 0; i < 6; i++)
	{
		for (uint8_t j = 0; j < 6; j++)
		{
			a_inv[i][j] = a_inv[i][j] / a[i][i];
		}
	}
	return 1;
}

void EKFquaternionToEuler(float q[4], float euler[3]) {
	float phi = atan2f(2 * (q[2] * q[3] + q[0] * q[1]), 1 - 2 * (q[1] * q[1] + q[2] * q[2]));

	float sinp = 2 * (q[1] * q[3] - q[0] * q[2]);
	float theta = 0.0;

	if (fabsf(sinp) >= 1) {
		theta = (sinp >= 0) ? M_PI / 2 : -M_PI / 2;
	}
	else {
		theta = -asinf(sinp);
	}

	float psi = atan2f(2 * (q[1] * q[2] + q[0] * q[3]), 1 - 2 * (q[2] * q[2] + q[3] * q[3]));

	euler[0] = phi;
	euler[1] = theta;
	euler[2] = psi;
}
