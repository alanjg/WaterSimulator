#pragma once

class vec2d
{
public:
	vec2d(double x = 0, double y = 0) { v[0] = x; v[1] = y; }
	vec2d(const vec2d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; }
	vec2d(double vec[2]) { v[0] = vec[0]; v[1] = vec[1]; }
	void Set(double x = 0, double y = 0) { v[0] = x; v[1] = y; }
	double& operator[] (int index) { return v[index]; }
	const double& operator[] (int index) const { return v[index]; }
	operator double*() { return &v[0]; }
	double Dot(const vec2d& rhs) const { return v[0] * rhs[0] + v[1] * rhs[1]; }
	vec2d &Normalize() { double z = sqrt(v[0] * v[0] + v[1] * v[1]); if (z != 0) { v[0] /= z; v[1] /= z; } return *this; }
	double Magnitude() const { return sqrt(v[0] * v[0] + v[1] * v[1]); }
	vec2d& operator*=(double rhs) { v[0] *= rhs; v[1] *= rhs; return *this; }
	vec2d& operator/=(double rhs) { v[0] /= rhs; v[1] /= rhs; return *this; }
	vec2d operator*(double rhs) { return vec2d(*this) *= rhs; }
	vec2d operator/(double rhs) { return vec2d(*this) /= rhs; }
	vec2d& operator+=(vec2d rhs) { v[0] += rhs[0]; v[1] += rhs[1]; return *this; }
	vec2d& operator-=(vec2d rhs) { v[0] -= rhs[0]; v[1] -= rhs[1]; return *this; }
	vec2d operator+(vec2d rhs) { return vec2d(*this) += rhs; }
	vec2d operator-(vec2d rhs) { return vec2d(*this) -= rhs; }
	bool equals(vec2d rhs, double dif) { return (*this - rhs).Magnitude()<dif; }
	friend std::ostream& operator<<(std::ostream&, const vec2d&);
	friend std::istream& operator>>(std::istream&, vec2d&);
private:
	double v[2];
};

class vec3d
{
public:
	vec3d(double x = 0, double y = 0, double z = 0) { v[0] = x; v[1] = y; v[2] = z; }
	vec3d(const vec3d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; v[2] = rhs[2]; }
	vec3d(double vec[3]) { v[0] = vec[0]; v[1] = vec[1]; v[2] = vec[2]; }
	vec3d(const vec2d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; v[2] = 0; }
	operator vec2d() { return vec2d(v[0], v[1]); }
	void Set(double x = 0, double y = 0, double z = 0) { v[0] = x; v[1] = y; v[2] = z; }
	double& operator[] (int index) { return v[index]; }
	const double& operator[] (int index) const { return v[index]; }
	operator double*() { return &v[0]; }
	operator const double*() { return &v[0]; }
	vec3d Cross(const vec3d& rhs) const { return vec3d(v[1] * rhs[2] - v[2] * rhs[1], v[2] * rhs[0] - v[0] * rhs[2], v[0] * rhs[1] - v[1] * rhs[0]); }
	double Dot(const vec3d& rhs) const { return v[0] * rhs[0] + v[1] * rhs[1] + v[2] * rhs[2]; }
	vec3d& Normalize() { double z = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); if (z != 0) { v[0] /= z; v[1] /= z; v[2] /= z; } return *this; }
	double Magnitude() const { return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
	vec3d& operator*=(double rhs) { v[0] *= rhs; v[1] *= rhs; v[2] *= rhs; return *this; }
	vec3d& operator/=(double rhs) { v[0] /= rhs; v[1] /= rhs; v[2] /= rhs; return *this; }
	vec3d operator*(double rhs) const { return vec3d(*this) *= rhs; }
	vec3d operator/(double rhs) const { return vec3d(*this) /= rhs; }
	vec3d& operator+=(vec3d rhs) { v[0] += rhs[0]; v[1] += rhs[1]; v[2] += rhs[2]; return *this; }
	vec3d& operator-=(vec3d rhs) { v[0] -= rhs[0]; v[1] -= rhs[1]; v[2] -= rhs[2]; return *this; }
	vec3d operator+(vec3d rhs) const { return vec3d(*this) += rhs; }
	vec3d operator-(vec3d rhs) const { return vec3d(*this) -= rhs; }
	vec3d operator-() const { return vec3d(-v[0], -v[1], -v[2]); }
	bool equals(vec3d rhs, double dif) { return (*this - rhs).Magnitude()<dif; }
	friend std::ostream& operator<<(std::ostream&, const vec3d&);
	friend std::istream& operator>>(std::istream&, vec3d&);
private:
	double v[3];
};

class vec4d
{
public:
	vec4d(double x = 0, double y = 0, double z = 0, double w = 0) { v[0] = x; v[1] = y; v[2] = z; v[3] = w; }
	vec4d(const vec4d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; v[2] = rhs[2]; v[3] = rhs[3]; }
	vec4d(double vec[4]) { v[0] = vec[0]; v[1] = vec[1]; v[2] = vec[2]; v[3] = vec[3]; }

	vec4d(const vec2d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; v[2] = 0; v[3] = 0; }
	vec4d(const vec3d& rhs) { v[0] = rhs[0]; v[1] = rhs[1]; v[2] = rhs[2]; v[3] = 0; }
	operator vec3d() { return vec3d(v[0], v[1], v[2]); }
	operator vec2d() { return vec2d(v[0], v[1]); }
	void Set(double x = 0, double y = 0, double z = 0, double w = 0) { v[0] = x; v[1] = y; v[2] = z; v[3] = w; }
	double& operator[] (int index) { return v[index]; }
	const double& operator[] (int index) const { return v[index]; }
	operator double*() { return &v[0]; }
	double Dot(const vec4d& rhs) const { return v[0] * rhs[0] + v[1] * rhs[1] + v[2] * rhs[2] + v[3] * rhs[3]; }
	vec4d& Normalize() { double z = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]); if (z != 0) { v[0] /= z; v[1] /= z; v[2] /= z; v[3] /= z; } return *this; }
	double Magnitude() const { return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]); }
	vec4d& operator*=(double rhs) { v[0] *= rhs; v[1] *= rhs; v[2] *= rhs; v[3] *= rhs; return *this; }
	vec4d& operator/=(double rhs) { v[0] /= rhs; v[1] /= rhs; v[2] /= rhs; v[3] /= rhs; return *this; }
	vec4d operator*(double rhs) const { return vec4d(*this) *= rhs; }
	vec4d operator/(double rhs) const { return vec4d(*this) /= rhs; }
	vec4d& operator+=(vec4d rhs) { v[0] += rhs[0]; v[1] += rhs[1]; v[2] += rhs[2]; v[3] += rhs[3]; return *this; }
	vec4d& operator-=(vec4d rhs) { v[0] -= rhs[0]; v[1] -= rhs[1]; v[2] -= rhs[2]; v[3] -= rhs[3]; return *this; }
	vec4d operator+(vec4d rhs) const { return vec4d(*this) += rhs; }
	vec4d operator-(vec4d rhs) const { return vec4d(*this) -= rhs; }
	bool equals(vec4d rhs, double dif) { return (*this - rhs).Magnitude()<dif; }
	friend std::ostream& operator<<(std::ostream&, const vec4d&);
	friend std::istream& operator>>(std::istream&, vec4d&);
private:
	double v[4];
};

struct line2d
{
	vec2d v[2];
};


struct triangle3d
{
	vec3d v[3];
};