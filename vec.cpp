#include "stdafx.h"
#include "vec.h"


#include "vec.h"

std::ostream& operator<<(std::ostream& out, const vec2d& v)
{
	return out << v[0] << " " << v[1];
}

std::istream& operator>>(std::istream& in, vec2d& v)
{
	return in >> v[0] >> v[1];
}

std::ostream& operator<<(std::ostream& out, const vec3d& v)
{
	return out << v[0] << " " << v[1] << " " << v[2];
}
std::istream& operator>>(std::istream& in, vec3d& v)
{
	return in >> v[0] >> v[1] >> v[2];
}

std::ostream& operator<<(std::ostream& out, const vec4d& v)
{
	return out << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
}

std::istream& operator>>(std::istream& in, vec4d& v)
{
	return in >> v[0] >> v[1] >> v[2] >> v[3];
}
