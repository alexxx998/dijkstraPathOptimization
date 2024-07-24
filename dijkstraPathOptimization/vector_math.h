#pragma once

class Vector3d {
public:
	Vector3d() {};
	Vector3d(double x, double y, double z) : p_x(x), p_y(y), p_z(z) {};

    double& x() { return p_x; }
    double x() const { return p_x; }
    double& y() { return p_y; }
    double y() const { return p_y; }
    double& z() { return p_z; }
    double z() const { return p_z; }

	Vector3d operator*(double k) const;
	Vector3d operator*(const Vector3d& v) const;
	Vector3d operator-(const Vector3d& v) const;
	Vector3d operator+(const Vector3d& v) const;
	double operator^(const Vector3d& v) const;
	double lenSqr() const;
	double len() const;
	
private:
	double p_x, p_y, p_z;
};

class Vector2d {
public:
    Vector2d() : p_x(0.0), p_y(0.0) {};
    Vector2d(double x, double y) : p_x(x), p_y(y) {};

    double& x() { return p_x; }
    double x() const { return p_x; }
    double& y() { return p_y; }
    double y() const { return p_y; }

    Vector2d operator*(double k) const;
    double operator*(const Vector2d& v) const;
    Vector2d operator-(const Vector2d& v) const;
    Vector2d operator+(const Vector2d& v) const;
    double operator^(const Vector2d& v) const;
    bool operator< (const Vector2d& p) const;
    double lenSqr() const;
    double len() const;
    Vector2d rotate(double ang) const;
private:
    double p_x, p_y;
};