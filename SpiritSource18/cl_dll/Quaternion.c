//quaternion.c
// Written by Jim "Entropy" Hunter

#include "quaternion.h"

// up / down
#define	PITCH	0
// left / right
#define	YAW		1
// fall over
#define	ROLL	2 

/////////////////////////
// Quaternion functions
/////////////////////////
// QuaternionNormalize //
/////////////////////////
void QuaternionNormalize(quaternion_t * q)
{
	double mag;

	mag = q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z;
	mag = sqrt(mag);

	q->w /= mag;
	q->x /= mag;
	q->y /= mag;
	q->z /= mag;
}

///////////////////////////////////
//  EulerToQuaternion			 //
// Assumes pitch, yaw, then roll //
///////////////////////////////////
void EulerToQuaternion(const float * angles, quaternion_t * q)
{
	double cp, cy, cr, sp, sy, sr;
	
	// need angle in radians divided by 2
	cp = cos(angles[PITCH] * M_PI / 360);
	cy = cos(angles[YAW] * M_PI / 360);
	cr = cos(angles[ROLL] * M_PI / 360);
		
	sp = sin(angles[PITCH] * M_PI / 360);
	sy = sin(angles[YAW] * M_PI / 360);
	sr = sin(angles[ROLL] * M_PI / 360);

	q->w = cp*cy*cr + sp*sy*sr;
	q->x = cp*cy*sr - sp*sy*cr;
	q->y = cr*sp*cy + sr*cp*sy;
	q->z = cr*cp*sy - sr*sp*cy;

	QuaternionNormalize(q);
}

//////////////////////////////////
//  QuaternionToEuler			//
// Assumes Pitch, yaw then roll //
//////////////////////////////////
void QuaternionToEuler(const quaternion_t q, float * angles)
{
	double yaw, pitch, roll;
	double sqw = q.w*q.w;
	double sqx = q.x*q.x;
	double sqy = q.y*q.y;
	double sqz = q.z*q.z;

	//yaw = atan(2.0 * (q.x * q.y + q.z * q.w) / (sqx - sqy - sqz + sqw));
	yaw = atan2(2.0 * (q.x * q.y + q.z * q.w), (sqx - sqy - sqz + sqw));
	//roll = atan(2.0 * (q.y * q.z + q.x * q.w) / (-sqx - sqy + sqz + sqw));
	roll = atan2(2.0 * (q.y * q.z + q.x * q.w), (-sqx - sqy + sqz + sqw));
	pitch = asin(-2.0 * (q.x * q.z - q.y * q.w));

	//Convert angles to degrees
	angles[YAW] = (float)(yaw * 180 / M_PI);
	angles[ROLL] = (float)(roll * 180 / M_PI);
	angles[PITCH] = (float)(pitch * 180 / M_PI);
}

////////////////////////////
// QuaternionMultiply //
////////////////////////////
void QuaternionMultiply(const quaternion_t q1, const quaternion_t q2, quaternion_t * result)
{	
	result->w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
	result->x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
	result->y = q1.w*q2.y + q1.y*q2.w + q1.z*q2.x - q1.x*q2.z;
	result->z = q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x;
}

//////////////////////////
// QuaternionConjugate  //
//////////////////////////
void QuaternionConjugate(const quaternion_t q1, quaternion_t * q2)
{
	q2->w = q1.w;
	q2->x = -q1.x;
	q2->y = -q1.y;
	q2->z = -q1.z;
}

////////////////////////////////
// QuaternionFromAxisAndAngle //
// angle is in degrees		  //
// axis is a vector and		  //
// should be normalized		  //
////////////////////////////////
void QuaternionFromAxisAndAngle(const float * axis, const float angle, quaternion_t * result)
{
	double sa, ca;
	
	sa = sin( angle * M_PI / 360.0 );
    ca = cos( angle * M_PI / 360.0 );

    result->x = axis[0] * sa;
    result->y = axis[1] * sa;
    result->z = axis[2] * sa;
    result->w = ca;

	QuaternionNormalize(result);
}

void PointRotate(const float * point, const quaternion_t * q, float * result)
{
	quaternion_t temp;
	quaternion_t qpoint;
	quaternion_t iq;

	qpoint.w = 0;
	qpoint.x = point[0];
	qpoint.y = point[1];
	qpoint.z = point[2];

	QuaternionMultiply(*q, qpoint, &temp);
	QuaternionConjugate(*q, &iq);
	QuaternionMultiply(temp, iq, &qpoint);

	result[0] = (float)qpoint.x;
	result[1] = (float)qpoint.y;
	result[2] = (float)qpoint.z;
}

void QuaternionToAxisAndAngle(const quaternion_t q, float * angle, float * axis)
{
	double sina2 = 1-q.w*q.w;
	
	*angle = (float)(acos(q.w) * 2 * 180 / M_PI);
	if (fabs(sina2) < 0.0001)
		sina2 = 1;
	axis[0] = (float)(q.x / sina2);
	axis[1] = (float)(q.y / sina2);
	axis[2] = (float)(q.z / sina2);
}

void AllVectorsToAngles(float * forward, float * right, float * up, float * angles)
{
	float	tmp, yaw, pitch;
	float planeup[3], roll;

	// Clamp |x| and |y| to > 0.0001
	tmp = (float)fabs(forward[0]);
	if (tmp < 0.0001)
		forward[0] = (float)(forward[0]/tmp * 0.0001);
	tmp = (float)fabs(forward[1]);
	if (tmp < 0.0001)
		forward[1] = (float)(forward[1]/tmp * 0.0001);
	
	if (forward[1] == 0 && forward[0] == 0)
	{
		yaw = 0;
		if (forward[2] > 0)
			pitch = 90;
		else
			pitch = 270;
	}
	else
	{
		yaw = (float)(atan2(forward[1], forward[0]) * 180 / M_PI);
		if (yaw < 0)
			yaw += 360;

		tmp = (float)sqrt (forward[0]*forward[0] + forward[1]*forward[1]);
		pitch = (float)-(atan2(forward[2], tmp) * 180 / M_PI);
		if (pitch < 0)
			pitch += 360;
	}
	
	angles[0] = pitch;
	angles[1] = yaw;
	
	angles[2] = 0;
	AngleVectors(angles, forward, 0, planeup);
	roll = AngleBetweenVectors(up, planeup);
	angles[ROLL] = roll;
}

void QuaternionToMatrix(const quaternion_t q, float matrix[3][4])
{
	matrix[0][0] = (float)(1 - 2*q.y*q.y - 2*q.z*q.z);
	matrix[1][0] = (float)(2*q.x*q.y - 2*q.z*q.w);
	matrix[2][0] = (float)(2*q.x*q.z + 2*q.y*q.w);

	matrix[0][1] = (float)(2*q.x*q.y + 2*q.z*q.w);
	matrix[1][1] = (float)(1 - 2*q.x*q.x - 2*q.z*q.z);
	matrix[2][1] = (float)(2*q.y*q.z - 2*q.x*q.w);

	matrix[0][2] = (float)(2*q.x*q.z - 2*q.y*q.w);
	matrix[1][2] = (float)(2*q.y*q.z + 2*q.x*q.w);
	matrix[2][2] = (float)(1 - 2*q.x*q.x - 2*q.y*q.y);

	matrix[0][3] = matrix[1][3] = matrix[2][3] = 0;
}

void QuaternionCopy(const quaternion_t qin, quaternion_t * qout)
{
	qout->w = qin.w;
	qout->x = qin.x;
	qout->y = qin.y;
	qout->z = qin.z;
}

void QuaternionRotate(const float * anglesIn, const float * deltaAngles, float * anglesOut)
{
	static quaternion_t qIn, qYaw, qPitch, qRoll, qDelta, qResult, qTemp;
	vec3_t forward, right, up;

	// Convert the input angles to quaternion notation and find the axes
	EulerToQuaternion(anglesIn, &qIn);
	AngleVectors(anglesIn, forward, right, up);

	// Calculate the additional rotations from the delta angles and the axes
	QuaternionFromAxisAndAngle(up, deltaAngles[YAW], &qYaw);
	QuaternionFromAxisAndAngle(right, -deltaAngles[PITCH], &qPitch);
	QuaternionFromAxisAndAngle(forward, deltaAngles[ROLL], &qRoll);

	// Get the total additional rotation
	QuaternionMultiply(qPitch, qYaw, &qTemp);
	QuaternionMultiply(qTemp, qRoll, &qDelta);

	// Calculate the new angles
	QuaternionMultiply(qDelta, qIn, &qResult);
	QuaternionToEuler(qResult, anglesOut);
}	
