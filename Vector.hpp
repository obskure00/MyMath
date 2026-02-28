#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <algorithm>
#include "Matrix.hpp"

//=======================================
//---------------2D Vector---------------
//=======================================

class Vector2 {
public:
    double x;
    double y;

    constexpr Vector2() : x(0.0), y(0.0) {}
    constexpr Vector2(double x_, double y_)
        : x(x_), y(y_) {}

    // =========================
    // Arithmetic Operators
    // =========================

    Vector2 operator+(const Vector2& v) const {
        return Vector2(x + v.x, y + v.y);
    }

    Vector2 operator-(const Vector2& v) const {
        return Vector2(x - v.x, y - v.y);
    }

    Vector2 operator*(double s) const {
        return Vector2(x * s, y * s);
    }

    Vector2 operator/(double s) const {
        return Vector2(x / s, y / s);
    }

    Vector2& operator+=(const Vector2& v) {
        x += v.x; y += v.y;
        return *this;
    }

    Vector2& operator-=(const Vector2& v) {
        x -= v.x; y -= v.y;
        return *this;
    }

    Vector2& operator*=(double s) {
        x *= s; y *= s;
        return *this;
    }

    Vector2& operator/=(double s) {
        x /= s; y /= s;
        return *this;
    }

    // =========================
    // Vector Math
    // =========================

    double dot2(const Vector2& v) const {
        return x*v.x + y*v.y;
    }

    double cross2(const Vector2& v) const {
        return x*v.y - y*v.x;
    }

    double lengthSquared2() const {
        return x*x + y*y;
    }

    double length2() const {
        return std::sqrt(lengthSquared2());
    }

    // =========================
    // Normalization
    // =========================

    Vector2 normalized2() const {
        double len = length2();
        if (len < 1e-12)
            return Vector2();
        return *this / len;
    }

    Vector2& normalize2() {
        double len = length2();
        if (len < 1e-12)
            return *this;
        return (*this /= len);
    }

    // =========================
    // Projection & Rejection
    // =========================

    Vector2 projectOnto2(const Vector2& b) const {
        double denom = b.lengthSquared2();
        if (denom < 1e-12)
            return Vector2();
        return (dot2(b) / denom) * b;
    }

    Vector2 rejectFrom2(const Vector2& b) const {
        return *this - projectOnto2(b);
    }

    // =========================
    // Reflection (normal must be normalized)
    // =========================

    Vector2 reflect2(const Vector2& normal) const {
        return *this - normal * (2.0 * dot2(normal));
    }

    // =========================
    // Angle Between
    // =========================

    double angleBetween2(const Vector2& v) const {
        double denom = length2() * v.length2();
        if (denom < 1e-12)
            return 0.0;

        double c = dot2(v) / denom;
        c = std::clamp(c, -1.0, 1.0);
        return std::acos(c);
    }
};

inline Vector2 operator*(double s, const Vector2& v) {
    return v * s;
}

Vector2 operator*(const Matrix2& m, const Vector2& v) {
    return Vector2(
        m.data[0]*v.x + m.data[1]*v.y,
        m.data[2]*v.x + m.data[3]*v.y
    );
}

//=======================================
//---------------3D Vector---------------
//=======================================

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

    double dot3(const Vector3& v) const {
        return x*v.x + y*v.y + z*v.z;
    }

    Vector3 cross3(const Vector3& v) const {
        return Vector3(
            y*v.z - z*v.y,
            z*v.x - x*v.z,
            x*v.y - y*v.x
        );
    }

    double lengthSquared3() const {
        return x*x + y*y + z*z;
    }

    double length3() const {
        return std::sqrt(lengthSquared3());
    }

    // =========================
    // Normalization
    // =========================

    Vector3 normalized3() const {
        double len = length3();
        if (len < 1e-12)
            return Vector3();
        return *this / len;
    }

    Vector3& normalize3() {
        double len = length3();
        if (len < 1e-12)
            return *this;
        return (*this /= len);
    }

    // =========================
    // Projection & Rejection
    // =========================

    Vector3 projectOnto3(const Vector3& b) const {
        double denom = b.lengthSquared3();
        if (denom < 1e-12)
            return Vector3();
        return (dot3(b) / denom) * b;
    }

    Vector3 rejectFrom(const Vector3& b) const {
        return *this - projectOnto3(b);
    }

    // =========================
    // Reflection
    // =========================

    Vector3 reflect(const Vector3& normal) const {
        return *this - normal * (2.0 * this->dot3(normal));
    }

    // =========================
    // Refraction (Snell's Law)
    // =========================

    Vector3 refract(const Vector3& normal, double eta) const {
        double cosi = std::clamp(this->dot3(normal), -1.0, 1.0);
        double k = 1.0 - eta*eta * (1.0 - cosi*cosi);

        if (k < 0.0)
            return Vector3();

        return (*this * eta) - normal * (eta * cosi + std::sqrt(k));
    }

    // =========================
    // Angle Between
    // =========================

    double angleBetween(const Vector3& v) const {
        double denom = length3() * v.length3();
        if (denom < 1e-12)
            return 0.0;

        double c = dot3(v) / denom;
        c = std::clamp(c, -1.0, 1.0);
        return std::acos(c);
    }

    // =========================
    // Scalar Triple Product
    // =========================

    double tripleProduct(const Vector3& b, const Vector3& c) const {
        return this->dot3(b.cross3(c));
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
        e1 = v.normalized3();

        Vector3 helper =
            (std::fabs(e1.x) < 0.9)
            ? Vector3(1.0, 0.0, 0.0)
            : Vector3(0.0, 1.0, 0.0);

        e2 = e1.cross3(helper).normalized3();
        e3 = e1.cross3(e2);
    }
};

inline Vector3 operator*(double s, const Vector3& v) {
    return v * s;
}

Vector3 operator*(const Matrix3& m, const Vector3& v) {
    return Vector3(
        m.data[0]*v.x + m.data[1]*v.y + m.data[2]*v.z,
        m.data[3]*v.x + m.data[4]*v.y + m.data[5]*v.z,
        m.data[6]*v.x + m.data[7]*v.y + m.data[8]*v.z
    );
}

//=======================================
//---------------4D Vector---------------
//=======================================

class Vector4 {
public:
    double x;
    double y;
    double z;
    double w;

    constexpr Vector4() : x(0), y(0), z(0), w(0) {}
    constexpr Vector4(double x_, double y_, double z_, double w_)
        : x(x_), y(y_), z(z_), w(w_) {}

    // =========================
    // Arithmetic Operators
    // =========================

    Vector4 operator+(const Vector4& v) const {
        return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
    }

    Vector4 operator-(const Vector4& v) const {
        return Vector4(x - v.x, y - v.y, z - v.z, w - v.w);
    }

    Vector4 operator*(double s) const {
        return Vector4(x * s, y * s, z * s, w * s);
    }

    Vector4 operator/(double s) const {
        return Vector4(x / s, y / s, z / s, w / s);
    }

    Vector4& operator+=(const Vector4& v) {
        x += v.x; y += v.y; z += v.z; w += v.w;
        return *this;
    }

    Vector4& operator-=(const Vector4& v) {
        x -= v.x; y -= v.y; z -= v.z; w -= v.w;
        return *this;
    }

    Vector4& operator*=(double s) {
        x *= s; y *= s; z *= s; w *= s;
        return *this;
    }

    Vector4& operator/=(double s) {
        x /= s; y /= s; z /= s; w /= s;
        return *this;
    }

    // =========================
    // Vector Math
    // =========================

    double dot4(const Vector4& v) const {
        return x*v.x + y*v.y + z*v.z + w*v.w;
    }

    double lengthSquared4() const {
        return dot4(*this);
    }

    double length4() const {
        return std::sqrt(lengthSquared4());
    }

    Vector4 normalized4() const {
        double len = length4();
        if (len < 1e-12)
            return Vector4();
        return *this / len;
    }

    Vector4& normalize4() {
        double len = length4();
        if (len < 1e-12)
            return *this;
        return (*this /= len);
    }

    double angleBetween4(const Vector4& v) const {
        double denom = length4() * v.length4();
        if (denom < 1e-12)
            return 0.0;

        double c = dot4(v) / denom;
        c = std::clamp(c, -1.0, 1.0);
        return std::acos(c);
    }
};

inline Vector4 operator*(double s, const Vector4& v) {
    return v * s;
}

Vector4 operator*(const Matrix4& m, const Vector4& v) {
    return Vector4(
        m.data[0]*v.x + m.data[1]*v.y + m.data[2]*v.z + m.data[3]*v.w,
        m.data[4]*v.x + m.data[5]*v.y + m.data[6]*v.z + m.data[7]*v.w,
        m.data[8]*v.x + m.data[9]*v.y + m.data[10]*v.z + m.data[11]*v.w,
        m.data[12]*v.x + m.data[13]*v.y + m.data[14]*v.z + m.data[15]*v.w
    );
}

#endif