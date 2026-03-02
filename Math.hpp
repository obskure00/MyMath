#ifndef MATH_HPP
#define MATH_HPP

// Math.hpp — the single header users include
// Internal include order:
//   1. Vector.hpp
//   2. Matrix.hpp
//   3. This file

#include "Vector.hpp"
#include "Matrix.hpp"

#include <cmath>
#include <stdexcept>
#include <ostream>

// ============================================================
//  Matrix * Vector  operator*  out-of-line bodies
// ============================================================

inline Vector2 Matrix2::operator*(const Vector2& v) const {
    return Vector2(
        data[0]*v.x + data[1]*v.y,
        data[2]*v.x + data[3]*v.y
    );
}

inline Vector3 Matrix3::operator*(const Vector3& v) const {
    return Vector3(
        data[0]*v.x + data[1]*v.y + data[2]*v.z,
        data[3]*v.x + data[4]*v.y + data[5]*v.z,
        data[6]*v.x + data[7]*v.y + data[8]*v.z
    );
}

inline Vector4 Matrix4::operator*(const Vector4& v) const {
    return Vector4(
        data[0]*v.x  + data[1]*v.y  + data[2]*v.z  + data[3]*v.w,
        data[4]*v.x  + data[5]*v.y  + data[6]*v.z  + data[7]*v.w,
        data[8]*v.x  + data[9]*v.y  + data[10]*v.z + data[11]*v.w,
        data[12]*v.x + data[13]*v.y + data[14]*v.z + data[15]*v.w
    );
}

// ============================================================
//  Matrix3  —  Vector3-dependent static factories
// ============================================================

inline Matrix3 Matrix3::fromAxisAngle3(const Vector3& axis, double angle) {
    return fromAxisAngle3(axis.x, axis.y, axis.z, angle);
}

// ============================================================
//  Matrix4  —  Vector3-dependent static factories & helpers
// ============================================================

inline Matrix4 Matrix4::fromTranslation3(const Vector3& t) {
    Matrix4 m = Identity();
    m(0,3) = t.x;
    m(1,3) = t.y;
    m(2,3) = t.z;
    return m;
}

inline Matrix4 Matrix4::fromScale3(const Vector3& s) {
    Matrix4 m;
    m(0,0) = s.x;
    m(1,1) = s.y;
    m(2,2) = s.z;
    m(3,3) = 1.0;
    return m;
}

inline Matrix4 Matrix4::composeTRS3(const Vector3& translation, const Matrix3& rotation, const Vector3& scale) {
    Matrix4 m;
    for (int row = 0; row < 3; ++row) {
        m(row, 0) = rotation(row, 0) * scale.x;
        m(row, 1) = rotation(row, 1) * scale.y;
        m(row, 2) = rotation(row, 2) * scale.z;
    }
    m(0,3) = translation.x;
    m(1,3) = translation.y;
    m(2,3) = translation.z;
    m(3,3) = 1.0;
    return m;
}

inline Vector3 Matrix4::transformPoint3(const Vector3& p) const {
    return Vector3(
        data[0]*p.x + data[1]*p.y + data[2]*p.z  + data[3],
        data[4]*p.x + data[5]*p.y + data[6]*p.z  + data[7],
        data[8]*p.x + data[9]*p.y + data[10]*p.z + data[11]
    );
}

inline Vector3 Matrix4::transformDirection3(const Vector3& d) const {
    return Vector3(
        data[0]*d.x + data[1]*d.y + data[2]*d.z,
        data[4]*d.x + data[5]*d.y + data[6]*d.z,
        data[8]*d.x + data[9]*d.y + data[10]*d.z
    );
}

// ============================================================
//  Axis-aligned rotation applied directly to a Vector3
// ============================================================

inline Vector3 rotationX3(double rad, const Vector3& v) { return rotationX3(rad) * v; }
inline Vector3 rotationY3(double rad, const Vector3& v) { return rotationY3(rad) * v; }
inline Vector3 rotationZ3(double rad, const Vector3& v) { return rotationZ3(rad) * v; }

// ============================================================
//  Quaternion  —  see Quaternion.hpp
// ============================================================

#include "Quaternion.hpp"

// ============================================================
//  Transform3  —  see Transformation.hpp
// ============================================================

#include "Transformation.hpp"

// ============================================================
//  MassProperties + inertia utilities  —  see Inertia.hpp
// ============================================================

#include "Inertia.hpp"

// ============================================================
//  Time integration helpers  —  see Integration.hpp
// ============================================================

#include "Integration.hpp"

// ============================================================
//  Collision math primitives  —  see Collision.hpp
// ============================================================

#include "Collision.hpp"

// ============================================================
//  Advanced linear algebra utilities  —  see LinAlg.hpp
// ============================================================

#include "LinAlg.hpp"

#endif