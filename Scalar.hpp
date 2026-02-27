#ifndef SCALAR_HPP
#define SCALAR_HPP

#include <limits>
#include <cstdint>
#include <cmath>

namespace math {

    inline constexpr double pi = 3.141592653589793;
    inline constexpr double tau = 6.283185307179586;
    inline constexpr double root2 = 1.41421356237309504;
    constexpr double machineEpsilon = std::numeric_limits<double>::epsilon();


    constexpr double degRad(double deg) noexcept {
        return deg * (pi / 180.0);
    }

    constexpr double radDeg(double rad) noexcept {
        return rad * (180.0 / pi);
    }

    template<typename T>
    constexpr T abs(T x) noexcept {
        static_assert(std::is_arithmetic_v<T>);
        return x < T(0) ? -x : x;
    }

    template<typename T>
    constexpr T min(T a, T b) noexcept {
        return (b < a) ? b : a;
    }

    template<typename T>
    constexpr T max(T a, T b) noexcept {
        return (a < b) ? b : a;
    }

    template<typename T>
    constexpr T clamp(T v, T lo, T hi) noexcept {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }

    template<typename T>
    constexpr T saturate(T v) noexcept {
        return clamp(v, T(0), T(1));
    }

    template<typename T, typename U>
    constexpr T lerp(T a, T b, U t) noexcept {
        return a + (b - a) * t;
    }

    template<typename T>
    constexpr T inverseLerp(T a, T b, T x) noexcept {
        return (x - a) / (b - a);
    }

    template<typename T>
    constexpr T smoothstep(T edge0, T edge1, T x) noexcept {
        T t = saturate((x - edge0) / (edge1 - edge0));
        return t * t * (T(3) - T(2) * t);
    }

    template<typename T>
    constexpr T sign(T x) noexcept {
        return (T(0) < x) - (x < T(0));
    }

    template<typename T>
    constexpr T square(T x) noexcept {
        return x * x;
    }

    template<typename T>
    constexpr T cube(T x) noexcept {
        return x * x * x;
    }

    template<typename T>
    constexpr T reciprocal(T x) noexcept {
        return T(1) / x;
    }

    inline float reciprocalFast(float x) noexcept {
        return 1.0f / x;
    }

    template<typename T>
    inline T invSqrt(T x) {
        return T(1) / std::sqrt(x);
    }

}

#endif