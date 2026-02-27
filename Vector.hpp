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

        
};

#endif