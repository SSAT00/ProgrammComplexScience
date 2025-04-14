#include"Sourse.h"

extern "C" __declspec(dllexport) void Kepler2Pos(double* params, double* pos) {
	double M = params[0];
	double a = params[1];
	double e = params[2];
	double i = params[3];
	double theta = params[4];
	double w = params[5];
	double mu = 3.986004418E5;

	double E = 0.0, E1 = M;
	double diff = 1.0;
	while (abs(diff) > 1e-9) {
		E = E1 - (E1 - e * sin(E1) - M) / (1 - e * cos(E1));
		diff = E - E1;
		E1 = E;
	}

	double tz = a * (cos(E) - e); //
	double et = a * sqrt(1 - pow(e, 2)) * sin(E); //

	double a1 = cos(w) * cos(theta) - sin(w) * sin(theta) * cos(i); //
	double a2 = cos(w) * sin(theta) + sin(w) * cos(theta) * cos(i); //
	double a3 = sin(w) * sin(i); //

	double b1 = -sin(w) * cos(theta) - cos(w) * sin(theta) * cos(i); //
	double b2 = -sin(w) * sin(theta) + cos(w) * cos(theta) * cos(i); //
	double b3 = cos(w) * sin(i); //

	pos[0] = tz * a1 + et * b1;
	pos[1] = tz * a2 + et * b2;
	pos[2] = tz * a3 + et * b3;

	double r = a * (1 - e * cos(E)); //
	double p = a * (1 - pow(e, 2)); //
	double v = 2 * atan2(sqrt((1 + e) / (1 - e)) * tan(E / 2), 1.0); //
	double Vr = sqrt(mu / p) * e * sin(v); //
	double Vn = sqrt(mu / p) * (1 + e * cos(v)); //
	double u = v + w;

	double g1 = -sin(u) * cos(theta) - cos(u) * sin(theta) * cos(i); //
	double g2 = -sin(u) * sin(theta) + cos(u) * cos(theta) * cos(i); //
	double g3 = cos(u) * sin(i); //

	pos[3] = pos[0] * Vr / r + g1 * Vn;
	pos[4] = pos[1] * Vr / r + g2 * Vn;
	pos[5] = pos[2] * Vr / r + g3 * Vn;
}

extern "C" __declspec(dllexport) void Pos2Kepler(double* pos, double* params) {
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
	double vx = pos[3];
	double vy = pos[4];
	double vz = pos[5];
	double mu = 3.986004418E5;

	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	double v2 = pow(vx, 2) + pow(vy, 2) + pow(vz, 2);

	double c1 = y * vz - z * vy; //
	double c2 = z * vx - x * vz; //
	double c3 = x * vy - y * vx; //

	double l1 = -mu * x / r + vy * c3 - vz * c2; //
	double l2 = -mu * y / r + vz * c1 - vx * c3; //
	double l3 = -mu * z / r + vx * c2 - vy * c1; //

	double c = sqrt(pow(c1, 2) + pow(c2, 2) + pow(c3, 2)); //
	double l = sqrt(pow(l1, 2) + pow(l2, 2) + pow(l3, 2)); //
	double h = v2 / 2 - mu / r; // 

	params[1] = -mu / (2 * h); // a //
	params[2] = l / mu; // e //
	params[4] = atan2(-c1 / c2, 1.0); // theta //
	params[3] = asin(c1 / (c * sin(params[4]))); // i //
	params[5] = atan2(l3 / (sin(params[3]) * (l1 * cos(params[4]) + l2 * sin(params[4]))), 1.0); // w //
	double u = atan2(z / (sin(params[3]) * (x * cos(params[4]) + y * sin(params[4]))), 1.0); //
	double v = u - params[5]; //
	double E = 2 * atan2(sqrt((1 - params[2]) / (1 + params[2])) * tan(v / 2), 1.0); // 
	params[0] = E - params[2] * sin(E); // M
}
