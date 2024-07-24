#include "vector_math.h"
#include <cmath>

Vector3d Vector3d::operator*(double k) const {
	return Vector3d(k * this->x(), k * this->y(), k * this->z());
}

Vector3d Vector3d::operator*(const Vector3d& v) const {
	return Vector3d(
		this->y() * v.z() - this->z() * v.y(),
		this->z() * v.x() - this->x() * v.z(),
		this->x() * v.y() - this->y() * v.x()
	);
}

Vector3d Vector3d::operator-(const Vector3d& v) const {
	return Vector3d(
		this->x() - v.x(),
		this->y() - v.y(),
		this->z() - v.z()
	);
}
Vector3d Vector3d::operator+(const Vector3d& v) const {
	return Vector3d(
		this->x() + v.x(),
		this->y() + v.y(),
		this->z() + v.z()
	);
}

double Vector3d::operator^(const Vector3d& v) const {
	return this->x() * v.x() + this->y() * v.y() + this->z() * v.z();
}

double Vector3d::lenSqr() const {
	return (*this) ^ (*this);
}

double Vector3d::len() const {
	return sqrtl(this->lenSqr());
}

Vector2d Vector2d::operator*(double k) const {
	return Vector2d(p_x * k, p_y * k);
}

double Vector2d::operator*(const Vector2d& v) const {
	return p_x * v.p_y - p_y * v.p_x;
}

Vector2d Vector2d::operator-(const Vector2d& v) const {
	return Vector2d(p_x - v.p_x, p_y - v.p_y);
}

Vector2d Vector2d::operator+(const Vector2d& v) const {
	return Vector2d(p_x + v.p_x, p_y + v.p_y);
}

double Vector2d::operator^(const Vector2d& v) const {
	return p_x * v.p_x + p_y * v.p_y;
}

bool Vector2d::operator< (const Vector2d& p) const {
	return p_x < p.p_x - 1e-9 || abs(p_x - p.p_x) < 1e-9 && p_y < p.p_y - 1e-9;
}

double Vector2d::lenSqr() const {
	return p_x * p_x + p_y * p_y;
}

double Vector2d::len() const {
	return std::sqrt(lenSqr());
}

Vector2d Vector2d::rotate(double ang) const
{
	return Vector2d(Vector2d(cos(ang), -sin(ang)) ^ (*this), Vector2d(sin(ang), cos(ang)) ^ (*this));
}
