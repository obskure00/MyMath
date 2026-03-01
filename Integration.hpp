#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

// Integration.hpp — time integration helpers for rigid-body physics.
// Include via Math.hpp (pulled in automatically) or directly.

#include "Math.hpp"

namespace Integrator {

// ============================================================
//  Quaternion derivative
// ============================================================

inline Quaternion quaternionDerivative(const Quaternion& q, const Vector3& omega) {
    Quaternion omegaQ(0.0, omega.x, omega.y, omega.z);
    return (omegaQ * q) * 0.5;
}

// ============================================================
//  Explicit (Forward) Euler
// ============================================================

inline void integrateEulerLinear(Vector3& pos, Vector3& vel,
                                  const Vector3& acc, double dt)
{
    pos += vel * dt;
    vel += acc * dt;
}

inline Quaternion integrateEulerOrientation(const Quaternion& q,
                                             const Vector3&    omega,
                                             double dt)
{
    return q.integrateAngularVelocity3(omega, dt);
}

// ============================================================
//  Semi-Implicit (Symplectic) Euler
// ============================================================

inline void integrateSemiImplicitEulerLinear(Vector3& pos, Vector3& vel,
                                              const Vector3& acc, double dt)
{
    vel += acc * dt;
    pos += vel * dt;
}

inline void integrateSemiImplicitEulerOrientation(Quaternion& q,
                                                   Vector3&    omega,
                                                   const Vector3& alpha,
                                                   double dt)
{
    omega += alpha * dt;
    q      = integrateEulerOrientation(q, omega, dt);
}

// ============================================================
//  Velocity Verlet
//  Mandatory call order:
//    Vector3 newPos = verletIntegratePosition(pos, vel, acc, dt);
//    Vector3 newAcc = evaluateForces(newPos);     // YOUR PHYSICS CODE
//    Vector3 newVel = verletIntegrateVelocity(vel, acc, newAcc, dt);
// ============================================================

inline Vector3 verletIntegratePosition(const Vector3& pos,
                                       const Vector3& vel,
                                       const Vector3& acc,
                                       double dt)
{
    return pos + vel * dt + acc * (0.5 * dt * dt);
}

inline Vector3 verletIntegrateVelocity(const Vector3& vel,
                                        const Vector3& accOld,
                                        const Vector3& accNew,
                                        double dt)
{
    return vel + (accOld + accNew) * (0.5 * dt);
}

// ============================================================
//  Vec3Pair - state type for (position, velocity) RK4
// ============================================================

struct Vec3Pair {
    Vector3 position;
    Vector3 velocity;

    Vec3Pair operator+(const Vec3Pair& o) const {
        return { position + o.position, velocity + o.velocity };
    }
    Vec3Pair operator*(double s) const {
        return { position * s, velocity * s };
    }
};

// ============================================================
//  Generic RK4
//    Use integrateRK4Orientation / integrateRK4OrientationVarying
//    below rather than calling this directly.
// ============================================================

template<typename State, typename DerivFn>
State integrateRK4(const State& y, double t, double dt, DerivFn f)
{
    const double half = dt * 0.5;

    State k1 = f(t,        y);
    State k2 = f(t + half, y + k1 * half);
    State k3 = f(t + half, y + k2 * half);
    State k4 = f(t + dt,   y + k3 * dt);

    return y + (k1 + (k2 + k3) * 2.0 + k4) * (dt / 6.0);
}

// ============================================================
//  RK4 - orientation under constant angular velocity
// ============================================================

inline Quaternion integrateRK4Orientation(const Quaternion& q,
                                           const Vector3&    omega,
                                           double dt)
{
    auto f = [&](double /*t*/, const Quaternion& q_) -> Quaternion {
        return quaternionDerivative(q_, omega);
    };
    return integrateRK4(q, 0.0, dt, f).normalized();
}

// ============================================================
//  RK4 - orientation under time-varying angular velocity
// ============================================================

template<typename OmegaFn>
Quaternion integrateRK4OrientationVarying(const Quaternion& q,
                                           double t, double dt,
                                           OmegaFn omegaFn)
{
    auto f = [&](double ti, const Quaternion& q_) -> Quaternion {
        return quaternionDerivative(q_, omegaFn(ti));
    };
    return integrateRK4(q, t, dt, f).normalized();
}

}

#endif