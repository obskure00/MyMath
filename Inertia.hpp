#ifndef INERTIA_HPP
#define INERTIA_HPP

// Inertia.hpp — mass properties and inertia tensor utilities.
// All tensor formulas assume uniform density.
// Symmetry of physical inertia tensors is exploited throughout — the optimised
// Matrix3::transformSymmetricInertiaTensor3 path is used everywhere.
// Include via Math.hpp pulled in directly.

#include "Math.hpp"

#include <cmath>
#include <stdexcept>

// ============================================================
//  MassProperties
// ============================================================

struct MassProperties {
    double  mass;
    double  inverseMass;
    Matrix3 inertiaTensorBody;
    Matrix3 inverseInertiaTensorBody;

    MassProperties()
        : mass(0.0)
        , inverseMass(0.0)
        , inertiaTensorBody()
        , inverseInertiaTensorBody()
    {}

    MassProperties(double m, const Matrix3& I)
        : mass(m)
        , inverseMass(m > 0.0 ? 1.0 / m : 0.0)
        , inertiaTensorBody(I)
        , inverseInertiaTensorBody(
            (m > 0.0 && std::abs(I.determinant3()) > 1e-9)
                ? I.inverse3()
                : Matrix3())
    {}
};

// ============================================================
//  computeSphereInertia
// ============================================================

inline MassProperties computeSphereInertia(double radius, double mass) {
    if (radius < 0.0)
        throw std::invalid_argument("computeSphereInertia: radius must be >= 0.");
    if (mass < 0.0)
        throw std::invalid_argument("computeSphereInertia: mass must be >= 0.");

    const double I = (2.0 / 5.0) * mass * radius * radius;

    return MassProperties(mass, Matrix3(
        I,   0.0, 0.0,
        0.0, I,   0.0,
        0.0, 0.0, I
    ));
}

// ============================================================
//  computeBoxInertia
// ============================================================

inline MassProperties computeBoxInertia(double width, double height, double depth, double mass) {
    if (width < 0.0 || height < 0.0 || depth < 0.0)
        throw std::invalid_argument("computeBoxInertia: dimensions must be >= 0.");
    if (mass < 0.0)
        throw std::invalid_argument("computeBoxInertia: mass must be >= 0.");

    const double w2 = width  * width;
    const double h2 = height * height;
    const double d2 = depth  * depth;
    const double k  = mass / 12.0;

    return MassProperties(mass, Matrix3(
        k * (h2 + d2), 0.0,           0.0,
        0.0,           k * (w2 + d2), 0.0,
        0.0,           0.0,           k * (w2 + h2)
    ));
}

// ============================================================
//  computeCapsuleInertia
// ============================================================

inline MassProperties computeCapsuleInertia(double radius, double halfHeight, double mass) {
    if (radius < 0.0)
        throw std::invalid_argument("computeCapsuleInertia: radius must be >= 0.");
    if (halfHeight < 0.0)
        throw std::invalid_argument("computeCapsuleInertia: halfHeight must be >= 0.");
    if (mass < 0.0)
        throw std::invalid_argument("computeCapsuleInertia: mass must be >= 0.");

    const double r  = radius;
    const double hl = halfHeight;
    const double r2 = r * r;

    const double vCyl    = 2.0 * hl;
    const double vSphere = (4.0 / 3.0) * r;
    const double vTotal  = vCyl + vSphere;

    double mCyl, mSphere;
    if (vTotal < 1e-30) {
        mCyl    = 0.0;
        mSphere = mass;
    } else {
        mCyl    = mass * vCyl    / vTotal;
        mSphere = mass * vSphere / vTotal;
    }

    const double Iyy = 0.5  * mCyl    * r2
                     + 0.4  * mSphere * r2;

    const double Ixx = mCyl    * (r2 / 4.0 + hl * hl / 3.0)
                     + mSphere * (0.4 * r2 + hl * hl + 0.75 * hl * r);

    return MassProperties(mass, Matrix3(
        Ixx, 0.0, 0.0,
        0.0, Iyy, 0.0,
        0.0, 0.0, Ixx
    ));
}

// ============================================================
//  transformInertiaToWorld
//  Called once per body per physics step to cache the result
// ============================================================

inline Matrix3 transformInertiaToWorld(const Matrix3& R, const Matrix3& Ibody) {
    return Matrix3::transformSymmetricInertiaTensor3(Ibody, R);
}

// ============================================================
//  computeInverseInertiaWorld
// ============================================================

inline Matrix3 computeInverseInertiaWorld(const Matrix3& R, const Matrix3& IbodyInv) {
    return Matrix3::transformSymmetricInertiaTensor3(IbodyInv, R);
}

#endif