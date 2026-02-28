#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <array>

class Matrix3 {
    private:
        std::array<double, 9> data;

    public:

        Matrix3() {
            data.fill(0.0);
        }

        explicit Matrix3(double initial) {
            data.fill(initial);
        }

        Matrix3(double m00, double m01, double m02,
                double m10, double m11, double m12,
                double m20, double m21, double m22)
                : data{m00, m01, m02, m10, m11, m12, m20, m21, m22} {}

        static Matrix3 Identity() {
            Matrix3 m;
            m.data[0] = 1.0; m.data[4] = 1.0; m.data[8] = 1.0;
            return m;
        }

        Matrix3 operator+(const Matrix3& m) const {
            Matrix3 result;
            for (int i = 0; i < 9; ++i)
                result.data[i] = data[i] + m.data[i];
            return result;
        }

        Matrix3 operator-(const Matrix3& m) const {
            Matrix3 result;
            for (int i = 0; i < 9; ++i)
                result.data[i] = data[i] - m.data[i];
            return result;
        }

        Matrix3 operator*(const Matrix3& m) const {
            Matrix3 result;

            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col) {
                    result.data[row * 3 + col] =
                        data[row * 3 + 0] * m.data[0 * 3 + col] +
                        data[row * 3 + 1] * m.data[1 * 3 + col] +
                        data[row * 3 + 2] * m.data[2 * 3 + col];
                }
            }

            return result;
        }

        Matrix3& operator+=(const Matrix3& m) {
            for (std::size_t i = 0; i < 9; ++i)
                data[i] += m.data[i];
            return *this;
        }

        Matrix3& operator-=(const Matrix3& m) {
            for (std::size_t i = 0; i < 9; ++i)
                data[i] -= m.data[i];
            return *this;
        }

        Matrix3& operator*= (const Matrix3& m) {
            *this = *this * m;
            return *this;
        }
};

#endif