#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <memory>
#include <vector>
#include <tuple>

class Vector3 {
    private:
        double x;
        double y;
        double z;

    public:
        Vector3 (double x, double y, double z)
            : x(x), y(y), z(z) {}

        std::tuple<double, double, double> getVector3 () {
            return std::make_tuple(x, y, z);
        }

        void setVector3 (double Newx, double Newy, double Newz) {
            x = Newx;
            y = Newy;
            z = Newz;
            return;
        }

        Vector3 operator+ (const Vector3& other) {
            return Vector3(x + other.x, y + other.y, z + other.z);
        }

        Vector3 operator- (const Vector3& other) {
            return Vector3(x - other.x, y - other.y, z - other.z);
        }

        Vector3 crossProd (const Vector3& other) {
            Vector3 vec(0, 0, 0);
            vec.x = (y * other.z) - (z * other.y);
            vec.y = (z * other.x) - (x * other.z);
            vec.z = (x * other.y) - (y * other.x);
            return vec;
        }

        void scalarProd (double scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return;
        }
};

#endif