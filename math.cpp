#include "stdafx.h"
#include "math.h"

float Deg2Rad(float deg)
{
	return deg*PI / 180.0f;
}

float Rad2Deg(float rad)
{
	return rad*180.0f / PI;
}

double CalcWENO(double v1, double v2, double v3, double v4, double v5)
{
	double phi1 = v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6;
	double phi2 = -v2 / 6 + 5 * v3 / 6 + v4 / 3;
	double phi3 = v3 / 3 + 5 * v4 / 6 - v5 / 6;

	double s1 = 13 / 12.0*square(v1 - 2 * v2 + v3) + 0.25*square(v1 - 4 * v2 + 3 * v3);
	double s2 = 13 / 12.0*square(v2 - 2 * v3 + v4) + 0.25*square(v2 - v4);
	double s3 = 13 / 12.0*square(v3 - 2 * v4 + v5) + 0.25*square(3 * v3 - 4 * v4 + v5);

	using std::max;
	double maxv = max(square(v1), max(square(v2), max(square(v3), max(square(v4), square(v5)))));

	double epsilon = 1e-6 * maxv + 1e-99;
	double alpha1 = 0.1 / square(s1 + epsilon);
	double alpha2 = 0.6 / square(s2 + epsilon);
	double alpha3 = 0.3 / square(s3 + epsilon);

	double alphaDenom = alpha1 + alpha2 + alpha3;
	double w1 = alpha1 / alphaDenom;
	double w2 = alpha2 / alphaDenom;
	double w3 = alpha3 / alphaDenom;

	return w1 * phi1 + w2 * phi2 + w3 * phi3;
}

double square(double x)
{
	return x*x;
}