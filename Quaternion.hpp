#ifndef QUATERNION_HPP
#define QUATERNION_HPP

// Depends on Vector3, Matrix3, and Matrix4 being fully defined.
// Therefore included last in Math.hpp

#include "Math.hpp"

#include <cmath>
#include <ostream>

// ============================================================
//  Quaternion
//  All rotation methods assume the quaternion is unit-length unless
//  otherwise noted. Call normalize() after accumulating many
//  small updates to prevent numerical drift.
// ============================================================

class Quaternion {
public:
    double w, x, y, z;

    constexpr Quaternion()
        : w(1.0), x(0.0), y(0.0), z(0.0) {}

    constexpr Quaternion(double w_, double x_, double y_, double z_)
        : w(w_), x(x_), y(y_), z(z_) {}

    // =========================
    // Static factories
    // =========================

    static Quaternion fromAxisAngle3(const Vector3& axis, double angle) {
        double half = angle * 0.5;
        double s    = std::sin(half);
        Vector3 u   = axis.normalized3();
        return Quaternion(std::cos(half), u.x*s, u.y*s, u.z*s);
    }

    static Quaternion fromMatrix3(const Matrix3& m) {
        double trace = m.data[0] + m.data[4] + m.data[8];
        Quaternion q;
        if (trace > 0.0) {
            double s = 0.5 / std::sqrt(trace + 1.0);
            q.w = 0.25 / s;
            q.x = (m.data[7] - m.data[5]) * s;
            q.y = (m.data[2] - m.data[6]) * s;
            q.z = (m.data[3] - m.data[1]) * s;
        } else if (m.data[0] > m.data[4] && m.data[0] > m.data[8]) {
            double s = 2.0 * std::sqrt(1.0 + m.data[0] - m.data[4] - m.data[8]);
            q.w = (m.data[7] - m.data[5]) / s;
            q.x = 0.25 * s;
            q.y = (m.data[1] + m.data[3]) / s;
            q.z = (m.data[2] + m.data[6]) / s;
        } else if (m.data[4] > m.data[8]) {
            double s = 2.0 * std::sqrt(1.0 + m.data[4] - m.data[0] - m.data[8]);
            q.w = (m.data[2] - m.data[6]) / s;
            q.x = (m.data[1] + m.data[3]) / s;
            q.y = 0.25 * s;
            q.z = (m.data[5] + m.data[7]) / s;
        } else {
            double s = 2.0 * std::sqrt(1.0 + m.data[8] - m.data[0] - m.data[4]);
            q.w = (m.data[3] - m.data[1]) / s;
            q.x = (m.data[2] + m.data[6]) / s;
            q.y = (m.data[5] + m.data[7]) / s;
            q.z = 0.25 * s;
        }
        return q;
    }

    // =========================
    // Unary
    // =========================

    Quaternion operator-() const { return Quaternion(-w, -x, -y, -z); }

    // =========================
    // Arithmetic
    // =========================

    Quaternion operator+(const Quaternion& q) const { return Quaternion(w+q.w, x+q.x, y+q.y, z+q.z); }
    Quaternion operator-(const Quaternion& q) const { return Quaternion(w-q.w, x-q.x, y-q.y, z-q.z); }
    Quaternion operator*(double s)            const { return Quaternion(w*s,   x*s,   y*s,   z*s);   }
    Quaternion operator/(double s)            const { return Quaternion(w/s,   x/s,   y/s,   z/s);   }

    Quaternion& operator+=(const Quaternion& q) { w+=q.w; x+=q.x; y+=q.y; z+=q.z; return *this; }
    Quaternion& operator-=(const Quaternion& q) { w-=q.w; x-=q.x; y-=q.y; z-=q.z; return *this; }
    Quaternion& operator*=(double s)            { w*=s; x*=s; y*=s; z*=s; return *this; }

    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y - x*q.z + y*q.w + z*q.x,
            w*q.z + x*q.y - y*q.x + z*q.w
        );
    }

    Quaternion& operator*=(const Quaternion& q) { *this = *this * q; return *this; }

    // =========================
    // Comparison
    // =========================

    bool operator==(const Quaternion& q) const { return w==q.w && x==q.x && y==q.y && z==q.z; }
    bool operator!=(const Quaternion& q) const { return !(*this == q); }

    bool equals(const Quaternion& q, double eps = 1e-9) const {
        return std::fabs(w-q.w) <= eps && std::fabs(x-q.x) <= eps
            && std::fabs(y-q.y) <= eps && std::fabs(z-q.z) <= eps;
    }

    // =========================
    // Properties
    // =========================

    double dot(const Quaternion& q) const { return w*q.w + x*q.x + y*q.y + z*q.z; }

    double normSquared() const { return w*w + x*x + y*y + z*z; }
    double norm()        const { return std::sqrt(normSquared()); }

    Quaternion conjugate() const { return Quaternion(w, -x, -y, -z); }

    Quaternion inverse() const {
        double ns = normSquared();
        if (ns < 1e-24) return Quaternion();
        return conjugate() / ns;
    }

    Quaternion normalized() const {
        double n = norm();
        if (n < 1e-12) return Quaternion();
        return *this / n;
    }

    Quaternion& normalize() {
        double n = norm();
        if (n >= 1e-12) { w/=n; x/=n; y/=n; z/=n; }
        return *this;
    }

    // =========================
    // Rotation
    // =========================

    Vector3 rotate3(const Vector3& v) const {
        Vector3 qv(x, y, z);
        Vector3 t = qv.cross3(v) * 2.0;
        return v + t * w + qv.cross3(t);
    }

    // =========================
    // Conversion
    // =========================

    Matrix3 toMatrix3() const {
        double x2=x*x, y2=y*y, z2=z*z;
        double xy=x*y, xz=x*z, yz=y*z;
        double wx=w*x, wy=w*y, wz=w*z;
        return Matrix3(
            1-2*(y2+z2),   2*(xy-wz),    2*(xz+wy),
              2*(xy+wz),   1-2*(x2+z2),  2*(yz-wx),
              2*(xz-wy),   2*(yz+wx),    1-2*(x2+y2)
        );
    }

    Matrix4 toMatrix4() const {
        return Matrix4::fromRotation3(toMatrix3());
    }

    // =========================
    // Physics integration
    // =========================

    Quaternion integrateAngularVelocity3(const Vector3& omega, double dt) const {
        Quaternion omegaQ(0.0, omega.x, omega.y, omega.z);
        Quaternion dq = (omegaQ * (*this)) * (0.5 * dt);
        return (*this + dq).normalized();
    }

    // =========================
    // Inertia Tensor Transform
    // =========================

    Matrix3 transformInertiaTensor3(const Matrix3& Ibody) const {
        Matrix3 R = toMatrix3();
        return Matrix3::transformInertiaTensor3(Ibody, R);
    }

    Matrix3 transformSymmetricInertiaTensor3(const Matrix3& Ibody) const {
        Matrix3 R = toMatrix3();
        return Matrix3::transformSymmetricInertiaTensor3(Ibody, R);
    }

    Matrix3 transformSymmetricInertiaTensor3(
        double ixx, double ixy, double ixz,
                    double iyy, double iyz,
                                double izz) const
    {
        Matrix3 R = toMatrix3();
        return Matrix3::transformSymmetricInertiaTensor3(ixx, ixy, ixz, iyy, iyz, izz, R);
    }

    // =========================
    // Interpolation
    // =========================

    static Quaternion slerp(const Quaternion& a, const Quaternion& b, double t) {
        double cosHalf = a.dot(b);

        Quaternion b2 = b;
        if (cosHalf < 0.0) { b2 = -b; cosHalf = -cosHalf; }

        if (cosHalf > 0.9995) {
            return Quaternion(
                a.w + t*(b2.w - a.w),
                a.x + t*(b2.x - a.x),
                a.y + t*(b2.y - a.y),
                a.z + t*(b2.z - a.z)
            ).normalized();
        }

        double halfAngle = std::acos(cosHalf);
        double sinHalf   = std::sqrt(1.0 - cosHalf*cosHalf);
        double wa = std::sin((1.0 - t) * halfAngle) / sinHalf;
        double wb = std::sin(t          * halfAngle) / sinHalf;
        return Quaternion(
            a.w*wa + b2.w*wb,
            a.x*wa + b2.x*wb,
            a.y*wa + b2.y*wb,
            a.z*wa + b2.z*wb
        );
    }

    static Quaternion nlerp(const Quaternion& a, const Quaternion& b, double t) {
        double cosHalf = a.dot(b);
        Quaternion b2  = (cosHalf < 0.0) ? -b : b;
        return Quaternion(
            a.w + t*(b2.w - a.w),
            a.x + t*(b2.x - a.x),
            a.y + t*(b2.y - a.y),
            a.z + t*(b2.z - a.z)
        ).normalized();
    }
};

inline Quaternion operator*(double s, const Quaternion& q) { return q * s; }

inline std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    return os << "Quaternion(w=" << q.w
              << ", x=" << q.x
              << ", y=" << q.y
              << ", z=" << q.z << ")";
}

#endif