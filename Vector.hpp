#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>
#include <algorithm>

class Vector3 {
public:
    double x;
    double y;
    double z;

    constexpr Vector3() : x(0.0), y(0.0), z(0.0) {}
    constexpr Vector3(double x_, double y_, double z_)
        : x(x_), y(y_), z(z_) {}

    // =========================
    // Basic Arithmetic Operators
    // =========================

    Vector3 operator+(const Vector3& v) const {
        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    Vector3 operator-(const Vector3& v) const {
        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    Vector3 operator*(double s) const {
        return Vector3(x * s, y * s, z * s);
    }

    Vector3 operator/(double s) const {
        return Vector3(x / s, y / s, z / s);
    }

    Vector3& operator+=(const Vector3& v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    Vector3& operator-=(const Vector3& v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Vector3& operator*=(double s) {
        x *= s; y *= s; z *= s;
        return *this;
    }

    Vector3& operator/=(double s) {
        x /= s; y /= s; z /= s;
        return *this;
    }

    // =========================
    // Vector Math
    // =========================

    double dot(const Vector3& v) const {
        return x*v.x + y*v.y + z*v.z;
    }

    Vector3 cross(const Vector3& v) const {
        return Vector3(
            y*v.z - z*v.y,
            z*v.x - x*v.z,
            x*v.y - y*v.x
        );
    }

    double lengthSquared() const {
        return x*x + y*y + z*z;
    }

    double length() const {
        return std::sqrt(lengthSquared());
    }

    // =========================
    // Normalization
    // =========================

    Vector3 normalized() const {
        double len = length();
        if (len < 1e-12)
            return Vector3();
        return *this / len;
    }

    Vector3& normalize() {
        double len = length();
        if (len < 1e-12)
            return *this;
        return (*this /= len);
    }

    // =========================
    // Projection & Rejection
    // =========================

    Vector3 projectOnto(const Vector3& b) const {
        double denom = b.lengthSquared();
        if (denom < 1e-12)
            return Vector3();
        return (dot(b) / denom) * b;
    }

    Vector3 rejectFrom(const Vector3& b) const {
        return *this - projectOnto(b);
    }

    // =========================
    // Reflection
    // =========================

    Vector3 reflect(const Vector3& normal) const {
        return *this - normal * (2.0 * this->dot(normal));
    }

    // =========================
    // Refraction (Snell's Law)
    // =========================

    Vector3 refract(const Vector3& normal, double eta) const {
        double cosi = std::clamp(this->dot(normal), -1.0, 1.0);
        double k = 1.0 - eta*eta * (1.0 - cosi*cosi);

        if (k < 0.0)
            return Vector3();

        return (*this * eta) - normal * (eta * cosi + std::sqrt(k));
    }

    // =========================
    // Angle Between
    // =========================

    double angleBetween(const Vector3& v) const {
        double denom = length() * v.length();
        if (denom < 1e-12)
            return 0.0;

        double c = dot(v) / denom;
        c = std::clamp(c, -1.0, 1.0);
        return std::acos(c);
    }

    // =========================
    // Scalar Triple Product
    // =========================

    double tripleProduct(const Vector3& b, const Vector3& c) const {
        return this->dot(b.cross(c));
    }

    // =========================
    // Orthonormal Basis
    // =========================

    static void orthonormalBasis(
        const Vector3& v,
        Vector3& e1,
        Vector3& e2,
        Vector3& e3)
    {
        e1 = v.normalized();

        Vector3 helper =
            (std::fabs(e1.x) < 0.9)
            ? Vector3(1.0, 0.0, 0.0)
            : Vector3(0.0, 1.0, 0.0);

        e2 = e1.cross(helper).normalized();
        e3 = e1.cross(e2);
    }
};

inline Vector3 operator*(double s, const Vector3& v) {
    return v * s;
}

#endif