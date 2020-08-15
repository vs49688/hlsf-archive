// quaternion.h
// Written by Jim "Entropy" Hunter

#ifndef QUATERNION_H
#define QUATERNION_H

#include <math.h>
#include "mathlib.h"

float AngleBetweenVectors( const vec3_t v1, const vec3_t v2 );

typedef struct {
	double w;
	double x;
	double y;
	double z;
} quaternion_t;

void QuaternionNormalize(quaternion_t * q);
void EulerToQuaternion(const float * angles, quaternion_t * q);
void QuaternionToEuler(const quaternion_t q, float * angles);
void QuaternionMultiply(const quaternion_t q1, const quaternion_t q2, quaternion_t * result);
void QuaternionConjugate(const quaternion_t q1, quaternion_t * q2);
void QuaternionFromAxisAndAngle(const float * axis, const float angle, quaternion_t * result);
void PointRotate(const float * point, const quaternion_t * q, float * result);
void QuaternionToMatrix(const quaternion_t q, float matrix[3][4]);

#endif