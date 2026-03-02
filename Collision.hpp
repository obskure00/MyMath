#ifndef COLLISION_HPP
#define COLLISION_HPP

// Collision.hpp — geometric primitives and collision math for rigid-body physics.
// Include via Math.hpp pulled in directly.

#include "Math.hpp"

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace Collision {

// ============================================================
//  GEOMETRIC PRIMITIVES
// ============================================================

struct Ray {
    Vector3 origin;
    Vector3 direction;

    Ray() = default;
    Ray(const Vector3& o, const Vector3& d) : origin(o), direction(d) {}

    Vector3 at(double t) const { return origin + direction * t; }
};

struct Plane {
    Vector3 normal;
    double  d;

    Plane() : normal(0,1,0), d(0.0) {}
    Plane(const Vector3& n, double d_) : normal(n), d(d_) {}

    static Plane fromNormalPoint(const Vector3& n, const Vector3& p) {
        Vector3 unit = n.normalized3();
        return Plane(unit, unit.dot3(p));
    }

    static Plane fromTriangle(const Vector3& a, const Vector3& b, const Vector3& c) {
        return fromNormalPoint((b - a).cross3(c - a), a);
    }

    double signedDistance(const Vector3& p) const { return normal.dot3(p) - d; }
};

struct AABB {
    Vector3 min;
    Vector3 max;

    AABB() : min(0,0,0), max(0,0,0) {}
    AABB(const Vector3& mn, const Vector3& mx) : min(mn), max(mx) {}

    static AABB fromCenterHalfExtents(const Vector3& c, const Vector3& h) {
        return AABB(c - h, c + h);
    }

    Vector3 center()      const { return (min + max) * 0.5; }
    Vector3 halfExtents() const { return (max - min) * 0.5; }
    Vector3 size()        const { return max - min; }

    bool contains(const Vector3& p) const {
        return p.x >= min.x && p.x <= max.x
            && p.y >= min.y && p.y <= max.y
            && p.z >= min.z && p.z <= max.z;
    }

    bool overlaps(const AABB& o) const {
        return min.x <= o.max.x && max.x >= o.min.x
            && min.y <= o.max.y && max.y >= o.min.y
            && min.z <= o.max.z && max.z >= o.min.z;
    }
};

struct OBB {
    Vector3 center;
    Vector3 halfExtents;
    Matrix3 rotation;

    OBB() : center(0,0,0), halfExtents(1,1,1), rotation(Matrix3::Identity()) {}

    OBB(const Vector3& c, const Vector3& h, const Matrix3& r)
        : center(c), halfExtents(h), rotation(r) {}

    static OBB fromAABB(const AABB& aabb) {
        return OBB(aabb.center(), aabb.halfExtents(), Matrix3::Identity());
    }
};

struct Sphere {
    Vector3 center;
    double  radius;

    Sphere() : center(0,0,0), radius(0.0) {}
    Sphere(const Vector3& c, double r) : center(c), radius(r) {}
};

struct Capsule {
    Vector3 a;
    Vector3 b;
    double  radius;

    Capsule() : a(0,0,0), b(0,1,0), radius(0.5) {}
    Capsule(const Vector3& a_, const Vector3& b_, double r) : a(a_), b(b_), radius(r) {}
};

struct Triangle {
    Vector3 a, b, c;

    Triangle() = default;
    Triangle(const Vector3& a_, const Vector3& b_, const Vector3& c_)
        : a(a_), b(b_), c(c_) {}

    Vector3 normal() const { return (b - a).cross3(c - a); }

    Vector3 unitNormal() const { return normal().normalized3(); }
};

struct ConvexPolytope {
    std::vector<Vector3> vertices;

    ConvexPolytope() = default;
    explicit ConvexPolytope(std::vector<Vector3> verts)
        : vertices(std::move(verts)) {}
};

// ============================================================
//  Internal helpers
// ============================================================

inline Vector3 obbAxis(const Matrix3& R, int j) {
    return Vector3(R.data[j], R.data[3+j], R.data[6+j]);
}

// ============================================================
//  MATH OPERATIONS
// ============================================================

// ============================================================
//  Closest point on segment
// ============================================================

inline Vector3 closestPointOnSegment(const Vector3& p,
                                      const Vector3& a, const Vector3& b,
                                      double* t_out = nullptr)
{
    Vector3 ab = b - a;
    double denom = ab.lengthSquared3();
    if (denom < 1e-24) {
        if (t_out) *t_out = 0.0;
        return a;
    }
    double t = std::clamp((p - a).dot3(ab) / denom, 0.0, 1.0);
    if (t_out) *t_out = t;
    return a + ab * t;
}

// ============================================================
//  Closest point on triangle
// ============================================================

inline Vector3 closestPointOnTriangle(const Vector3& p, const Triangle& tri)
{
    const Vector3& A = tri.a;
    const Vector3& B = tri.b;
    const Vector3& C = tri.c;

    Vector3 ab = B - A;
    Vector3 ac = C - A;
    Vector3 ap = p - A;

    double d1 = ab.dot3(ap), d2 = ac.dot3(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return A;
    Vector3 bp = p - B;
    double d3 = ab.dot3(bp), d4 = ac.dot3(bp);
    if (d3 >= 0.0 && d4 <= d3) return B;
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        return A + ab * v;
    }

    Vector3 cp = p - C;
    double d5 = ab.dot3(cp), d6 = ac.dot3(cp);
    if (d6 >= 0.0 && d5 <= d6) return C;

    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        return A + ac * w;
    }

    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return B + (C - B) * w;
    }

    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    return A + ab * v + ac * w;
}

// ============================================================
//  Closest point on AABB
//
//  Component-wise clamp of p to [min, max].
// ============================================================

inline Vector3 closestPointOnAABB(const Vector3& p, const AABB& aabb) {
    return Vector3(
        std::clamp(p.x, aabb.min.x, aabb.max.x),
        std::clamp(p.y, aabb.min.y, aabb.max.y),
        std::clamp(p.z, aabb.min.z, aabb.max.z)
    );
}

// ============================================================
//  Closest point on OBB
//
//  Transform p to OBB local space, clamp each axis independently,
//  transform back to world space.
// ============================================================

inline Vector3 closestPointOnOBB(const Vector3& p, const OBB& obb)
{
    Vector3 d = p - obb.center;
    Vector3 q = obb.center;

    const double he[3] = { obb.halfExtents.x, obb.halfExtents.y, obb.halfExtents.z };

    for (int i = 0; i < 3; ++i) {
        Vector3 axis = obbAxis(obb.rotation, i);
        double dist = std::clamp(d.dot3(axis), -he[i], he[i]);
        q = q + axis * dist;
    }
    return q;
}

// ============================================================
//  Barycentric coordinates
// ============================================================

struct BarycentricCoords {
    double u, v, w;
};

inline BarycentricCoords barycentricCoords(const Vector3& p, const Triangle& tri)
{
    Vector3 ab = tri.b - tri.a;
    Vector3 ac = tri.c - tri.a;
    Vector3 ap = p     - tri.a;

    double d00 = ab.dot3(ab);
    double d01 = ab.dot3(ac);
    double d11 = ac.dot3(ac);
    double d20 = ap.dot3(ab);
    double d21 = ap.dot3(ac);

    double denom = d00 * d11 - d01 * d01;
    if (std::fabs(denom) < 1e-24) {
        return { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
    }

    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    return { 1.0 - v - w, v, w };
}

inline bool pointInTriangle(const Vector3& p, const Triangle& tri) {
    auto bc = barycentricCoords(p, tri);
    return bc.u >= 0.0 && bc.v >= 0.0 && bc.w >= 0.0;
}

// ============================================================
//  Distance — point to plane (signed)
// ============================================================

inline double distancePointPlane(const Vector3& p, const Plane& plane) {
    return plane.signedDistance(p);
}

// ============================================================
//  Distance — point to segment
// ============================================================

inline double distancePointSegment(const Vector3& p,
                                    const Vector3& a, const Vector3& b)
{
    return (p - closestPointOnSegment(p, a, b)).length3();
}

// ============================================================
//  Distance — point to AABB
// ============================================================

inline double distancePointAABB(const Vector3& p, const AABB& aabb) {
    return (p - closestPointOnAABB(p, aabb)).length3();
}

// ============================================================
//  Distance — point to OBB
// ============================================================

inline double distancePointOBB(const Vector3& p, const OBB& obb) {
    return (p - closestPointOnOBB(p, obb)).length3();
}

// ============================================================
//  Closest points between two segments
// ============================================================

struct SegmentSegmentResult {
    double  t1;
    double  t2;
    double  distSquared;
    Vector3 p1;
    Vector3 p2;
};

inline SegmentSegmentResult closestPointsSegmentSegment(
    const Vector3& a1, const Vector3& b1,
    const Vector3& a2, const Vector3& b2)
{
    Vector3 d1 = b1 - a1;
    Vector3 d2 = b2 - a2;
    Vector3 r  = a1 - a2;

    double a = d1.lengthSquared3();
    double e = d2.lengthSquared3();
    double f = d2.dot3(r);

    double s, t;

    if (a < 1e-24 && e < 1e-24) {
        s = 0.0; t = 0.0;
    } else if (a < 1e-24) {
        s = 0.0;
        t = std::clamp(f / e, 0.0, 1.0);
    } else {
        double c = d1.dot3(r);
        if (e < 1e-24) {
            t = 0.0;
            s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            double b    = d1.dot3(d2);
            double denom = a * e - b * b;

            if (denom > 1e-24) {
                s = std::clamp((b * f - c * e) / denom, 0.0, 1.0);
            } else {
                s = 0.0;
            }

            t = (b * s + f) / e;

            if (t < 0.0) {
                t = 0.0;
                s = std::clamp(-c / a, 0.0, 1.0);
            } else if (t > 1.0) {
                t = 1.0;
                s = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }

    Vector3 p1 = a1 + d1 * s;
    Vector3 p2 = a2 + d2 * t;
    double dist2 = (p1 - p2).lengthSquared3();
    return { s, t, dist2, p1, p2 };
}

// ============================================================
//  Ray intersection results
// ============================================================

struct RayHit {
    bool    hit    = false;
    double  t      = 0.0;
    Vector3 point;
    Vector3 normal;
};

// ============================================================
//  Ray vs Plane
// ============================================================

inline RayHit rayVsPlane(const Ray& ray, const Plane& plane) {
    double denom = plane.normal.dot3(ray.direction);
    if (std::fabs(denom) < 1e-12) return {};

    double t = (plane.d - plane.normal.dot3(ray.origin)) / denom;
    if (t < 0.0) return {};

    RayHit h;
    h.hit    = true;
    h.t      = t;
    h.point  = ray.at(t);
    h.normal = (denom < 0.0) ? plane.normal : -plane.normal;
    return h;
}

// ============================================================
//  Ray vs Sphere
// ============================================================

inline RayHit rayVsSphere(const Ray& ray, const Sphere& sphere) {
    Vector3 oc = ray.origin - sphere.center;
    double a = ray.direction.lengthSquared3();
    double b = 2.0 * oc.dot3(ray.direction);
    double c = oc.lengthSquared3() - sphere.radius * sphere.radius;
    double disc = b * b - 4.0 * a * c;

    if (disc < 0.0) return {};

    double sqrtDisc = std::sqrt(disc);
    double t = (-b - sqrtDisc) / (2.0 * a);
    if (t < 0.0) t = (-b + sqrtDisc) / (2.0 * a);
    if (t < 0.0) return {};

    RayHit h;
    h.hit    = true;
    h.t      = t;
    h.point  = ray.at(t);
    h.normal = (h.point - sphere.center).normalized3();
    return h;
}

// ============================================================
//  Ray vs AABB  (slab method)
// ============================================================

inline RayHit rayVsAABB(const Ray& ray, const AABB& aabb) {
    double tmin = -std::numeric_limits<double>::infinity();
    double tmax =  std::numeric_limits<double>::infinity();
    int    normalAxis = 0;
    bool   normalFlip = false;

    const double* orig = &ray.origin.x;
    const double* dir  = &ray.direction.x;
    const double* mn   = &aabb.min.x;
    const double* mx   = &aabb.max.x;

    for (int i = 0; i < 3; ++i) {
        if (std::fabs(dir[i]) < 1e-12) {
            if (orig[i] < mn[i] || orig[i] > mx[i]) return {};
        } else {
            double invD = 1.0 / dir[i];
            double t1 = (mn[i] - orig[i]) * invD;
            double t2 = (mx[i] - orig[i]) * invD;
            bool   flipped = t1 > t2;
            if (flipped) std::swap(t1, t2);
            if (t1 > tmin) {
                tmin       = t1;
                normalAxis = i;
                normalFlip = flipped;
            }
            tmax = std::min(tmax, t2);
            if (tmin > tmax) return {};
        }
    }

    if (tmax < 0.0) return {};

    double t = (tmin >= 0.0) ? tmin : tmax;
    Vector3 n(0,0,0);
    (&n.x)[normalAxis] = normalFlip ? 1.0 : -1.0;

    RayHit h;
    h.hit    = true;
    h.t      = t;
    h.point  = ray.at(t);
    h.normal = n;
    return h;
}

// ============================================================
//  Ray vs OBB
// ============================================================

inline RayHit rayVsOBB(const Ray& ray, const OBB& obb) {
    Vector3 diff = ray.origin - obb.center;
    Matrix3 RT   = obb.rotation.transpose3();

    Vector3 localOrigin    = RT * diff;
    Vector3 localDirection = RT * ray.direction;

    double tmin = -std::numeric_limits<double>::infinity();
    double tmax =  std::numeric_limits<double>::infinity();
    int    normalAxis = 0;
    bool   normalFlip = false;

    const double he[3] = { obb.halfExtents.x, obb.halfExtents.y, obb.halfExtents.z };
    const double* lo   = &localOrigin.x;
    const double* ld   = &localDirection.x;

    for (int i = 0; i < 3; ++i) {
        if (std::fabs(ld[i]) < 1e-12) {
            if (lo[i] < -he[i] || lo[i] > he[i]) return {};
        } else {
            double invD = 1.0 / ld[i];
            double t1   = (-he[i] - lo[i]) * invD;
            double t2   = ( he[i] - lo[i]) * invD;
            bool   flipped = t1 > t2;
            if (flipped) std::swap(t1, t2);
            if (t1 > tmin) {
                tmin       = t1;
                normalAxis = i;
                normalFlip = flipped;
            }
            tmax = std::min(tmax, t2);
            if (tmin > tmax) return {};
        }
    }

    if (tmax < 0.0) return {};

    double t = (tmin >= 0.0) ? tmin : tmax;

    Vector3 localNormal(0,0,0);
    (&localNormal.x)[normalAxis] = normalFlip ? 1.0 : -1.0;
    Vector3 worldNormal = obb.rotation * localNormal;

    RayHit h;
    h.hit    = true;
    h.t      = t;
    h.point  = ray.at(t);
    h.normal = worldNormal;
    return h;
}

// ============================================================
//  Ray vs Triangle  (Moller-Trumbore, two-sided)
// ============================================================

inline RayHit rayVsTriangle(const Ray& ray, const Triangle& tri) {
    constexpr double EPS = 1e-10;

    Vector3 e1 = tri.b - tri.a;
    Vector3 e2 = tri.c - tri.a;
    Vector3 h  = ray.direction.cross3(e2);
    double  a  = e1.dot3(h);

    if (std::fabs(a) < EPS) return {};

    double  f = 1.0 / a;
    Vector3 s = ray.origin - tri.a;
    double  u = f * s.dot3(h);
    if (u < 0.0 || u > 1.0) return {};

    Vector3 q = s.cross3(e1);
    double  v = f * ray.direction.dot3(q);
    if (v < 0.0 || u + v > 1.0) return {};

    double t = f * e2.dot3(q);
    if (t < EPS) return {};

    Vector3 geomNormal = e1.cross3(e2).normalized3();
    Vector3 normal = (a < 0.0) ? -geomNormal : geomNormal;

    RayHit h2;
    h2.hit    = true;
    h2.t      = t;
    h2.point  = ray.at(t);
    h2.normal = normal;
    return h2;
}

// ============================================================
//  Ray vs Capsule
// ============================================================

inline RayHit rayVsCapsule(const Ray& ray, const Capsule& cap) {
    Vector3 seg = cap.b - cap.a;
    Vector3 d   = ray.direction;
    Vector3 w   = ray.origin - cap.a;

    double dd = d.dot3(d);
    double ds = d.dot3(seg);
    double ss = seg.dot3(seg);
    double dw = d.dot3(w);
    double sw = seg.dot3(w);
    double ww = w.dot3(w);

    if (ss < 1e-24) {
        Sphere sph(cap.a, cap.radius);
        return rayVsSphere(ray, sph);
    }

    double denom = dd * ss - ds * ds;

    RayHit best;

    if (std::fabs(denom) > 1e-12) {
        double a2 = denom;
        double a1 = ss * (2.0 * dw) - ds * (2.0 * sw);
        double a0 = ss * (ww - cap.radius * cap.radius) - sw * sw;
        double disc = a1 * a1 - 4.0 * a2 * a0;

        if (disc >= 0.0) {
            double sqrtDisc = std::sqrt(disc);
            for (int sign : {-1, 1}) {
                double t = (-a1 + sign * sqrtDisc) / (2.0 * a2);
                if (t < 0.0) continue;
                double u = (sw + ds * t) / ss;
                if (u < 0.0 || u > 1.0) continue;
                if (!best.hit || t < best.t) {
                    best.hit    = true;
                    best.t      = t;
                    best.point  = ray.at(t);
                    Vector3 segPt = cap.a + seg * u;
                    best.normal = (best.point - segPt).normalized3();
                }
                break;
            }
        }
    }

    for (const Vector3* ep : { &cap.a, &cap.b }) {
        RayHit sh = rayVsSphere(ray, Sphere(*ep, cap.radius));
        if (sh.hit && (!best.hit || sh.t < best.t)) {
            Vector3 cp = closestPointOnSegment(sh.point, cap.a, cap.b);
            if ((sh.point - cp).lengthSquared3() >= cap.radius * cap.radius - 1e-9)
                best = sh;
        }
    }

    return best;
}

// ============================================================
//  SAT — interval helpers
// ============================================================

struct Interval {
    double min, max;
};

inline bool intervalsOverlap(const Interval& a, const Interval& b) {
    return a.max >= b.min && b.max >= a.min;
}

inline double intervalPenetration(const Interval& a, const Interval& b) {
    return std::min(a.max, b.max) - std::max(a.min, b.min);
}

inline Interval projectAABBOnAxis(const AABB& aabb, const Vector3& axis) {
    Vector3 c = aabb.center();
    Vector3 h = aabb.halfExtents();
    double  r = std::fabs(h.x * axis.x)
              + std::fabs(h.y * axis.y)
              + std::fabs(h.z * axis.z);
    double  cProj = c.dot3(axis);
    return { cProj - r, cProj + r };
}

inline Interval projectOBBOnAxis(const OBB& obb, const Vector3& axis) {
    double r = 0.0;
    const double he[3] = { obb.halfExtents.x, obb.halfExtents.y, obb.halfExtents.z };
    for (int i = 0; i < 3; ++i)
        r += he[i] * std::fabs(obbAxis(obb.rotation, i).dot3(axis));
    double cProj = obb.center.dot3(axis);
    return { cProj - r, cProj + r };
}

// ============================================================
//  SAT result
// ============================================================

struct SATResult {
    bool    separated  = false;
    double  penetration = 0.0;
    Vector3 axis;
};

// ============================================================
//  SAT — AABB vs AABB
// ============================================================

inline SATResult satAABBvsAABB(const AABB& a, const AABB& b) {
    SATResult result;
    result.penetration = std::numeric_limits<double>::infinity();

    const Vector3 axes[3] = { {1,0,0}, {0,1,0}, {0,0,1} };
    for (const Vector3& axis : axes) {
        Interval ia = projectAABBOnAxis(a, axis);
        Interval ib = projectAABBOnAxis(b, axis);
        if (!intervalsOverlap(ia, ib)) {
            result.separated   = true;
            result.penetration = 0.0;
            result.axis        = axis;
            return result;
        }
        double pen = intervalPenetration(ia, ib);
        if (pen < result.penetration) {
            result.penetration = pen;
            result.axis        = axis;
        }
    }
    return result;
}

// ============================================================
//  SAT — OBB vs OBB  (15 separating axes)
//  Implements the Gottschalk/Ericson formulation
//  (Real-Time Collision Detection).
// ============================================================

inline SATResult satOBBvsOBB(const OBB& A, const OBB& B) {
    constexpr double EPS = 1e-6;

    const Vector3 a[3] = { obbAxis(A.rotation, 0),
                            obbAxis(A.rotation, 1),
                            obbAxis(A.rotation, 2) };
    const Vector3 b_ax[3] = { obbAxis(B.rotation, 0),
                               obbAxis(B.rotation, 1),
                               obbAxis(B.rotation, 2) };

    const double ha[3] = { A.halfExtents.x, A.halfExtents.y, A.halfExtents.z };
    const double hb[3] = { B.halfExtents.x, B.halfExtents.y, B.halfExtents.z };

    double R[3][3], AbsR[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            R[i][j]    = a[i].dot3(b_ax[j]);
            AbsR[i][j] = std::fabs(R[i][j]) + EPS;
        }

    Vector3 T_world = B.center - A.center;
    double t[3] = { T_world.dot3(a[0]),
                    T_world.dot3(a[1]),
                    T_world.dot3(a[2]) };

    SATResult result;
    result.penetration = std::numeric_limits<double>::infinity();

    auto testAxis = [&](const Vector3& axis, double sepDist, double rA, double rB) -> bool {
        double lenSq = axis.lengthSquared3();
        if (lenSq < EPS * EPS) return false;
        double pen  = rA + rB - sepDist;
        double pen_normalised = pen / std::sqrt(lenSq);
        if (pen_normalised < 0.0) {
            result.separated   = true;
            result.penetration = 0.0;
            result.axis        = axis;
            return true;
        }
        if (pen_normalised < result.penetration) {
            result.penetration = pen_normalised;
            result.axis        = axis / std::sqrt(lenSq);
        }
        return false;
    };

    for (int i = 0; i < 3; ++i) {
        double rB = hb[0]*AbsR[i][0] + hb[1]*AbsR[i][1] + hb[2]*AbsR[i][2];
        if (testAxis(a[i], std::fabs(t[i]), ha[i], rB)) return result;
    }

    for (int j = 0; j < 3; ++j) {
        double tProj = t[0]*R[0][j] + t[1]*R[1][j] + t[2]*R[2][j];
        double rA    = ha[0]*AbsR[0][j] + ha[1]*AbsR[1][j] + ha[2]*AbsR[2][j];
        if (testAxis(b_ax[j], std::fabs(tProj), rA, hb[j])) return result;
    }

    {
        // A0 x B0
        double sep = std::fabs(t[2]*R[1][0] - t[1]*R[2][0]);
        double rA  = ha[1]*AbsR[2][0] + ha[2]*AbsR[1][0];
        double rB  = hb[1]*AbsR[0][2] + hb[2]*AbsR[0][1];
        Vector3 axis = a[0].cross3(b_ax[0]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A0 x B1
        double sep = std::fabs(t[2]*R[1][1] - t[1]*R[2][1]);
        double rA  = ha[1]*AbsR[2][1] + ha[2]*AbsR[1][1];
        double rB  = hb[0]*AbsR[0][2] + hb[2]*AbsR[0][0];
        Vector3 axis = a[0].cross3(b_ax[1]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A0 x B2
        double sep = std::fabs(t[2]*R[1][2] - t[1]*R[2][2]);
        double rA  = ha[1]*AbsR[2][2] + ha[2]*AbsR[1][2];
        double rB  = hb[0]*AbsR[0][1] + hb[1]*AbsR[0][0];
        Vector3 axis = a[0].cross3(b_ax[2]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A1 x B0
        double sep = std::fabs(t[0]*R[2][0] - t[2]*R[0][0]);
        double rA  = ha[0]*AbsR[2][0] + ha[2]*AbsR[0][0];
        double rB  = hb[1]*AbsR[1][2] + hb[2]*AbsR[1][1];
        Vector3 axis = a[1].cross3(b_ax[0]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A1 x B1
        double sep = std::fabs(t[0]*R[2][1] - t[2]*R[0][1]);
        double rA  = ha[0]*AbsR[2][1] + ha[2]*AbsR[0][1];
        double rB  = hb[0]*AbsR[1][2] + hb[2]*AbsR[1][0];
        Vector3 axis = a[1].cross3(b_ax[1]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A1 x B2
        double sep = std::fabs(t[0]*R[2][2] - t[2]*R[0][2]);
        double rA  = ha[0]*AbsR[2][2] + ha[2]*AbsR[0][2];
        double rB  = hb[0]*AbsR[1][1] + hb[1]*AbsR[1][0];
        Vector3 axis = a[1].cross3(b_ax[2]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A2 x B0
        double sep = std::fabs(t[1]*R[0][0] - t[0]*R[1][0]);
        double rA  = ha[0]*AbsR[1][0] + ha[1]*AbsR[0][0];
        double rB  = hb[1]*AbsR[2][2] + hb[2]*AbsR[2][1];
        Vector3 axis = a[2].cross3(b_ax[0]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A2 x B1
        double sep = std::fabs(t[1]*R[0][1] - t[0]*R[1][1]);
        double rA  = ha[0]*AbsR[1][1] + ha[1]*AbsR[0][1];
        double rB  = hb[0]*AbsR[2][2] + hb[2]*AbsR[2][0];
        Vector3 axis = a[2].cross3(b_ax[1]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }
    {
        // A2 x B2
        double sep = std::fabs(t[1]*R[0][2] - t[0]*R[1][2]);
        double rA  = ha[0]*AbsR[1][2] + ha[1]*AbsR[0][2];
        double rB  = hb[0]*AbsR[2][1] + hb[1]*AbsR[2][0];
        Vector3 axis = a[2].cross3(b_ax[2]);
        if (testAxis(axis, sep, rA, rB)) return result;
    }

    return result;
}

// ============================================================
//  GJK Support functions
// ============================================================

inline Vector3 supportSphere(const Sphere& s, const Vector3& dir) {
    return s.center + dir.normalized3() * s.radius;
}

inline Vector3 supportAABB(const AABB& aabb, const Vector3& dir) {
    return Vector3(
        dir.x >= 0.0 ? aabb.max.x : aabb.min.x,
        dir.y >= 0.0 ? aabb.max.y : aabb.min.y,
        dir.z >= 0.0 ? aabb.max.z : aabb.min.z
    );
}

inline Vector3 supportOBB(const OBB& obb, const Vector3& dir) {
    Vector3 result = obb.center;
    const double he[3] = { obb.halfExtents.x, obb.halfExtents.y, obb.halfExtents.z };
    for (int i = 0; i < 3; ++i) {
        Vector3 axis = obbAxis(obb.rotation, i);
        double  proj = dir.dot3(axis);
        result = result + axis * (proj >= 0.0 ? he[i] : -he[i]);
    }
    return result;
}

inline Vector3 supportCapsule(const Capsule& cap, const Vector3& dir) {
    Vector3 tip = (dir.dot3(cap.b - cap.a) >= 0.0) ? cap.b : cap.a;
    return tip + dir.normalized3() * cap.radius;
}

inline Vector3 supportTriangle(const Triangle& tri, const Vector3& dir) {
    double da = dir.dot3(tri.a);
    double db = dir.dot3(tri.b);
    double dc = dir.dot3(tri.c);
    if (da >= db && da >= dc) return tri.a;
    if (db >= dc)             return tri.b;
    return tri.c;
}

inline Vector3 supportConvex(const ConvexPolytope& poly, const Vector3& dir) {
    if (poly.vertices.empty())
        throw std::runtime_error("supportConvex: empty vertex set");
    double  best = -std::numeric_limits<double>::infinity();
    Vector3 result;
    for (const Vector3& v : poly.vertices) {
        double d = v.dot3(dir);
        if (d > best) { best = d; result = v; }
    }
    return result;
}

inline Vector3 supportMinkowskiDiff(const Vector3& sA, const Vector3& sB) {
    return sA - sB;
}

// ============================================================
//  EPA helpers
//  Usage:
//    1. Build initial polytope from GJK simplex:
//         EPAPolytope epa = makeEPAPolytope(v0, v1, v2, v3);
//    2. Loop:
//         const EPAFace& f = findClosestFace(epa);
//         Vector3 sup = supportMinkowskiDiff(
//                           supportA(dir),
//                           supportB(-dir))    // dir = f.normal
//         if (sup.dot3(f.normal) - f.dist < EPA_EPSILON) break;
//         expandPolytope(epa, sup);
//    3. Penetration depth = f.dist, direction = f.normal.
// ============================================================

struct EPAFace {
    int     i0, i1, i2;
    Vector3 normal;
    double  dist;
};

struct EPAPolytope {
    std::vector<Vector3> vertices;
    std::vector<EPAFace> faces;
};

inline EPAFace makeEPAFace(const std::vector<Vector3>& verts, int i0, int i1, int i2, const Vector3& interior) {
    const Vector3& v0 = verts[i0];
    const Vector3& v1 = verts[i1];
    const Vector3& v2 = verts[i2];

    Vector3 n = (v1 - v0).cross3(v2 - v0);
    double  lenSq = n.lengthSquared3();

    EPAFace f;
    if (lenSq < 1e-24) {
        f.i0     = i0; f.i1 = i1; f.i2 = i2;
        f.normal = Vector3(0,1,0);
        f.dist   = std::numeric_limits<double>::infinity();
        return f;
    }

    n = n / std::sqrt(lenSq);

    if (n.dot3(interior - v0) > 0.0) {
        n = -n;
        std::swap(f.i1, f.i2);
    }
    f.i0     = i0;
    f.i1     = i1;
    f.i2     = i2;
    if (n.dot3(interior - v0) > 0.0) n = -n;
    f.normal = n;
    f.dist   = std::fabs(n.dot3(v0));

    return f;
}

inline EPAPolytope makeEPAPolytope(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    EPAPolytope poly;
    poly.vertices = { a, b, c, d };

    poly.faces.push_back(makeEPAFace(poly.vertices, 0, 1, 2, d));
    poly.faces.push_back(makeEPAFace(poly.vertices, 0, 1, 3, c));
    poly.faces.push_back(makeEPAFace(poly.vertices, 0, 2, 3, b));
    poly.faces.push_back(makeEPAFace(poly.vertices, 1, 2, 3, a));

    return poly;
}

inline const EPAFace& findClosestFace(const EPAPolytope& poly) {
    if (poly.faces.empty())
        throw std::runtime_error("findClosestFace: empty polytope");

    size_t bestIdx = 0;
    for (size_t i = 1; i < poly.faces.size(); ++i)
        if (poly.faces[i].dist < poly.faces[bestIdx].dist)
            bestIdx = i;

    return poly.faces[bestIdx];
}

inline bool expandPolytope(EPAPolytope& poly, const Vector3& newVert) {
    const double EPS = 1e-10;

    const EPAFace& closest = findClosestFace(poly);
    if (newVert.dot3(closest.normal) - closest.dist < EPS)
        return false;

    int newIdx = static_cast<int>(poly.vertices.size());
    poly.vertices.push_back(newVert);

    using Edge = std::pair<int,int>;
    std::vector<Edge> visibleEdges;
    std::vector<bool> faceVisible(poly.faces.size(), false);

    for (size_t k = 0; k < poly.faces.size(); ++k) {
        const EPAFace& f = poly.faces[k];
        const Vector3& v = poly.vertices[f.i0];
        if (f.normal.dot3(newVert - v) > 0.0) {
            faceVisible[k] = true;
            visibleEdges.push_back({ f.i0, f.i1 });
            visibleEdges.push_back({ f.i1, f.i2 });
            visibleEdges.push_back({ f.i2, f.i0 });
        }
    }

    std::vector<Edge> horizon;
    for (const Edge& e : visibleEdges) {
        bool isBoundary = true;
        for (const Edge& other : visibleEdges) {
            if (other.first == e.second && other.second == e.first) {
                isBoundary = false;
                break;
            }
        }
        if (isBoundary) horizon.push_back(e);
    }

    std::vector<EPAFace> kept;
    kept.reserve(poly.faces.size());
    for (size_t k = 0; k < poly.faces.size(); ++k)
        if (!faceVisible[k])
            kept.push_back(poly.faces[k]);
    poly.faces = std::move(kept);

    Vector3 interior = poly.vertices[0];
    if (interior.equals3(newVert, 1e-15) && poly.vertices.size() > 1)
        interior = poly.vertices[1];

    for (const Edge& e : horizon) {
        EPAFace nf = makeEPAFace(poly.vertices, e.first, e.second, newIdx, interior);
        poly.faces.push_back(nf);
    }

    return true;
}

}

#endif