#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

// Transformation.hpp — rigid-body transform for the physics simulation loop.
// Depends on Quaternion (and transitively on Vector3, Matrix3, Matrix4).
// Include via Math.hpp, which pulls this in directly.

#include "Quaternion.hpp"

// ============================================================
//  Transform3
// ============================================================

class Transform3 {
    public:
        Vector3    position;
        Quaternion orientation;

        // =========================
        // Constructors
        // =========================

        Transform3()
            : position(0.0, 0.0, 0.0)
            , orientation() {}

        Transform3(const Vector3& pos, const Quaternion& orient)
            : position(pos)
            , orientation(orient) {}

        static Transform3 Identity() {
            return Transform3();
        }

        // =========================
        // Transform a point
        // =========================

        Vector3 transformPoint3(const Vector3& p) const {
            return orientation.rotate3(p) + position;
        }

        // =========================
        // Transform a direction
        // =========================

        Vector3 transformDirection3(const Vector3& d) const {
            return orientation.rotate3(d);
        }

        // =========================
        // Inverse
        // =========================

        Transform3 inverse() const {
            Quaternion invQ = orientation.conjugate();
            Vector3    invT = invQ.rotate3(-position);
            return Transform3(invT, invQ);
        }

        // =========================
        // Combine  (operator*)
        // =========================

        Transform3 operator*(const Transform3& other) const {
            return Transform3(
                orientation.rotate3(other.position) + position,
                orientation * other.orientation
            );
        }

        Transform3& operator*=(const Transform3& other) {
            *this = *this * other;
            return *this;
        }

        Transform3 combine(const Transform3& other) const {
            return *this * other;
        }

        // =========================
        // Interpolation
        // =========================

        Transform3 interpolate(const Transform3& other, double t) const {
            return Transform3(
                lerp3(position, other.position, t),
                Quaternion::slerp(orientation, other.orientation, t)
            );
        }

        Transform3 interpolateFast(const Transform3& other, double t) const {
            return Transform3(
                lerp3(position, other.position, t),
                Quaternion::nlerp(orientation, other.orientation, t)
            );
        }

        // =========================
        // Convert to Matrix4  (for rendering / GPU upload only)
        // =========================

        Matrix4 toMatrix4() const {
            Matrix4 m = orientation.toMatrix4();
            m(0, 3) = position.x;
            m(1, 3) = position.y;
            m(2, 3) = position.z;
            return m;
        }

        // =========================
        // Comparison
        // =========================

        bool operator==(const Transform3& other) const {
            return position == other.position && orientation == other.orientation;
        }

        bool operator!=(const Transform3& other) const {
            return !(*this == other);
        }

        bool equals(const Transform3& other, double eps = 1e-9) const {
            return position.equals3(other.position, eps)
                && orientation.equals(other.orientation, eps);
        }

        // =========================
        // Normalise orientation
        // =========================

        Transform3& normalizeOrientation() {
            orientation.normalize();
            return *this;
        }
};

inline std::ostream& operator<<(std::ostream& os, const Transform3& t) {
    return os << "Transform3(position=" << t.position
              << ", orientation=" << t.orientation << ")";
}

#endif