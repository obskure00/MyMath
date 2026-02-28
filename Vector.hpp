#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <algorithm>
#include <ostream>

class Matrix2;
class Matrix3;
class Matrix4;

//=======================================
//---------------- Vector2 --------------
//=======================================

class Vector2 {
    public:
        double x, y;

        constexpr Vector2()                      : x(0.0), y(0.0) {}
        constexpr Vector2(double x_, double y_)  : x(x_),  y(y_)  {}

        // =========================
        // Unary
        // =========================

        Vector2 operator-() const { return Vector2(-x, -y); }

        // =========================
        // Arithmetic
        // =========================

        Vector2  operator+ (const Vector2& v) const { return Vector2(x+v.x, y+v.y); }
        Vector2  operator- (const Vector2& v) const { return Vector2(x-v.x, y-v.y); }
        Vector2  operator* (double s)         const { return Vector2(x*s,   y*s);   }
        Vector2  operator/ (double s)         const { return Vector2(x/s,   y/s);   }

        Vector2& operator+=(const Vector2& v) { x+=v.x; y+=v.y; return *this; }
        Vector2& operator-=(const Vector2& v) { x-=v.x; y-=v.y; return *this; }
        Vector2& operator*=(double s)         { x*=s;   y*=s;   return *this; }
        Vector2& operator/=(double s)         { x/=s;   y/=s;   return *this; }

        // =========================
        // Comparison
        // =========================

        bool operator==(const Vector2& v) const { return x == v.x && y == v.y; }
        bool operator!=(const Vector2& v) const { return !(*this == v); }

        bool equals2(const Vector2& v, double eps = 1e-9) const {
            return std::fabs(x-v.x) <= eps && std::fabs(y-v.y) <= eps;
        }

        // =========================
        // Vector Math
        // =========================

        double dot2(const Vector2& v)   const { return x*v.x + y*v.y; }
        double cross2(const Vector2& v) const { return x*v.y - y*v.x; }

        double lengthSquared2() const { return x*x + y*y; }
        double length2()        const { return std::sqrt(lengthSquared2()); }

        // =========================
        // Normalization
        // =========================

        Vector2 normalized2() const {
            double len = length2();
            return (len < 1e-12) ? Vector2() : *this / len;
        }

        Vector2& normalize2() {
            double len = length2();
            if (len >= 1e-12) *this /= len;
            return *this;
        }

        // =========================
        // Projection & Rejection
        // =========================

        Vector2 projectOnto2(const Vector2& b) const {
            double denom = b.lengthSquared2();
            return (denom < 1e-12) ? Vector2() : b * (dot2(b) / denom);
        }

        Vector2 rejectFrom2(const Vector2& b) const {
            return *this - projectOnto2(b);
        }

        // =========================
        // Reflection
        // =========================

        Vector2 reflect2(const Vector2& normal) const {
            return *this - normal * (2.0 * dot2(normal));
        }

        // =========================
        // Angle Between
        // =========================

        double angleBetween2(const Vector2& v) const {
            double denom = length2() * v.length2();
            if (denom < 1e-12) return 0.0;
            return std::acos(std::clamp(dot2(v) / denom, -1.0, 1.0));
        }
};

inline Vector2  operator*(double s, const Vector2& v) { return v * s; }

inline Vector2  lerp2(const Vector2& a, const Vector2& b, double t) {
    return a + (b - a) * t;
}

inline std::ostream& operator<<(std::ostream& os, const Vector2& v) {
    return os << "Vector2(" << v.x << ", " << v.y << ")";
}

//=======================================
//---------------- Vector3 --------------
//=======================================

class Vector3 {
    public:
        double x, y, z;

        constexpr Vector3()                               : x(0.0), y(0.0), z(0.0) {}
        constexpr Vector3(double x_, double y_, double z_) : x(x_),  y(y_),  z(z_)  {}

        // =========================
        // Unary
        // =========================

        Vector3 operator-() const { return Vector3(-x, -y, -z); }

        // =========================
        // Arithmetic
        // =========================

        Vector3  operator+ (const Vector3& v) const { return Vector3(x+v.x, y+v.y, z+v.z); }
        Vector3  operator- (const Vector3& v) const { return Vector3(x-v.x, y-v.y, z-v.z); }
        Vector3  operator* (double s)         const { return Vector3(x*s,   y*s,   z*s);   }
        Vector3  operator/ (double s)         const { return Vector3(x/s,   y/s,   z/s);   }

        Vector3& operator+=(const Vector3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
        Vector3& operator-=(const Vector3& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
        Vector3& operator*=(double s)         { x*=s;   y*=s;   z*=s;   return *this; }
        Vector3& operator/=(double s)         { x/=s;   y/=s;   z/=s;   return *this; }

        // =========================
        // Comparison
        // =========================

        bool operator==(const Vector3& v) const { return x==v.x && y==v.y && z==v.z; }
        bool operator!=(const Vector3& v) const { return !(*this == v); }

        bool equals3(const Vector3& v, double eps = 1e-9) const {
            return std::fabs(x-v.x) <= eps
                && std::fabs(y-v.y) <= eps
                && std::fabs(z-v.z) <= eps;
        }

        // =========================
        // Vector Math
        // =========================

        double dot3(const Vector3& v) const { return x*v.x + y*v.y + z*v.z; }

        Vector3 cross3(const Vector3& v) const {
            return Vector3(
                y*v.z - z*v.y,
                z*v.x - x*v.z,
                x*v.y - y*v.x
            );
        }

        double lengthSquared3() const { return x*x + y*y + z*z; }
        double length3()        const { return std::sqrt(lengthSquared3()); }

        // =========================
        // Normalization
        // =========================

        Vector3 normalized3() const {
            double len = length3();
            return (len < 1e-12) ? Vector3() : *this / len;
        }

        Vector3& normalize3() {
            double len = length3();
            if (len >= 1e-12) *this /= len;
            return *this;
        }

        // =========================
        // Projection & Rejection
        // =========================

        Vector3 projectOnto3(const Vector3& b) const {
            double denom = b.lengthSquared3();
            return (denom < 1e-12) ? Vector3() : b * (dot3(b) / denom);
        }

        Vector3 rejectFrom3(const Vector3& b) const {
            return *this - projectOnto3(b);
        }

        // =========================
        // Reflection
        // =========================

        Vector3 reflect3(const Vector3& normal) const {
            return *this - normal * (2.0 * dot3(normal));
        }

        // =========================
        // Refraction  (Snell's Law)
        // =========================

        Vector3 refract3(const Vector3& normal, double eta) const {
            double cosTheta = std::clamp(dot3(normal), -1.0, 1.0);
            double k = 1.0 - eta * eta * (1.0 - cosTheta * cosTheta);
            if (k < 0.0) return Vector3();
            return (*this * eta) - normal * (eta * cosTheta + std::sqrt(k));
        }

        // =========================
        // Angle Between
        // =========================

        double angleBetween3(const Vector3& v) const {
            double denom = length3() * v.length3();
            if (denom < 1e-12) return 0.0;
            return std::acos(std::clamp(dot3(v) / denom, -1.0, 1.0));
        }

        // =========================
        // Scalar Triple Product
        // =========================

        double tripleProduct3(const Vector3& b, const Vector3& c) const {
            return dot3(b.cross3(c));
        }

        // =========================
        // Orthonormal Basis from a single non-zero vector
        // =========================

        static void orthonormalBasis3(
            const Vector3& v,
            Vector3& e1,
            Vector3& e2,
            Vector3& e3)
        {
            e1 = v.normalized3();
            Vector3 helper = (std::fabs(e1.x) < 0.9)
                        ? Vector3(1.0, 0.0, 0.0)
                        : Vector3(0.0, 1.0, 0.0);
            e2 = e1.cross3(helper).normalized3();
            e3 = e1.cross3(e2);
        }
};

inline Vector3 operator*(double s, const Vector3& v) { return v * s; }

inline Vector3 lerp3(const Vector3& a, const Vector3& b, double t) {
    return a + (b - a) * t;
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    return os << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
}

//=======================================
//---------------- Vector4 --------------
//=======================================

class Vector4 {
    public:
        double x, y, z, w;

        constexpr Vector4()                                           : x(0), y(0), z(0), w(0) {}
        constexpr Vector4(double x_, double y_, double z_, double w_) : x(x_), y(y_), z(z_), w(w_) {}

        // =========================
        // Unary
        // =========================

        Vector4 operator-() const { return Vector4(-x, -y, -z, -w); }

        // =========================
        // Arithmetic
        // =========================

        Vector4  operator+ (const Vector4& v) const { return Vector4(x+v.x, y+v.y, z+v.z, w+v.w); }
        Vector4  operator- (const Vector4& v) const { return Vector4(x-v.x, y-v.y, z-v.z, w-v.w); }
        Vector4  operator* (double s)         const { return Vector4(x*s,   y*s,   z*s,   w*s);   }
        Vector4  operator/ (double s)         const { return Vector4(x/s,   y/s,   z/s,   w/s);   }

        Vector4& operator+=(const Vector4& v) { x+=v.x; y+=v.y; z+=v.z; w+=v.w; return *this; }
        Vector4& operator-=(const Vector4& v) { x-=v.x; y-=v.y; z-=v.z; w-=v.w; return *this; }
        Vector4& operator*=(double s)         { x*=s;   y*=s;   z*=s;   w*=s;   return *this; }
        Vector4& operator/=(double s)         { x/=s;   y/=s;   z/=s;   w/=s;   return *this; }

        // =========================
        // Comparison
        // =========================

        bool operator==(const Vector4& v) const { return x==v.x && y==v.y && z==v.z && w==v.w; }
        bool operator!=(const Vector4& v) const { return !(*this == v); }

        bool equals4(const Vector4& v, double eps = 1e-9) const {
            return std::fabs(x-v.x) <= eps && std::fabs(y-v.y) <= eps
                && std::fabs(z-v.z) <= eps && std::fabs(w-v.w) <= eps;
        }

        // =========================
        // Vector Math
        // =========================

        double dot4(const Vector4& v) const { return x*v.x + y*v.y + z*v.z + w*v.w; }

        double lengthSquared4() const { return dot4(*this); }
        double length4()        const { return std::sqrt(lengthSquared4()); }

        Vector4 normalized4() const {
            double len = length4();
            return (len < 1e-12) ? Vector4() : *this / len;
        }

        Vector4& normalize4() {
            double len = length4();
            if (len >= 1e-12) *this /= len;
            return *this;
        }

        double angleBetween4(const Vector4& v) const {
            double denom = length4() * v.length4();
            if (denom < 1e-12) return 0.0;
            return std::acos(std::clamp(dot4(v) / denom, -1.0, 1.0));
        }
};

inline Vector4 operator*(double s, const Vector4& v) { return v * s; }

inline Vector4 lerp4(const Vector4& a, const Vector4& b, double t) {
    return a + (b - a) * t;
}

inline std::ostream& operator<<(std::ostream& os, const Vector4& v) {
    return os << "Vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
}

#endif