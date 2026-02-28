#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <array>
#include <initializer_list>
#include <algorithm>
#include <stdexcept>

class Matrix3x3 {
    private:
        std::array<double, 9> data;

    public:

        Matrix3x3() {
            data.fill(0.0);
        }

        explicit Matrix3x3(double initial) {
            data.fill(initial);
        }

        Matrix3x3(double m00, double m01, double m02,
                double m10, double m11, double m12,
                double m20, double m21, double m22)
                : data{m00, m01, m02, m10, m11, m12, m20, m21, m22} {}

        static Matrix3x3 Identity() {
            Matrix3x3 m;
            m.data[0] = 1.0; m.data[4] = 1.0; m.data[8] = 1.0;
            return m;
        }
};

#endif