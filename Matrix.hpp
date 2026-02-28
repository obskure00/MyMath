#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <array>
#include <stdexcept>

class Matrix3 {
    public:
        std::array<double, 9> data;

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

    // =========================
    // Arithmetic Operators
    // =========================

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

        Matrix3 operator*(const double scalar) const {
            Matrix3 result;

            for (int i = 0; i < 9; ++i) {
                result.data[i] = data[i] * scalar;
            }

            return result;
        }
    
        Matrix3 operator/(const double scalar) const {
            Matrix3 result;

            for (int i = 0; i < 9; ++i) {
                result.data[i] = data[i] / scalar;
            }

            return result;
        }

        Matrix3& operator*=(const double scalar) {
            for (int i = 0; i < 9; ++i) {
                data[i] *= scalar;
            }
            return *this;
        }

        Matrix3& operator+=(const Matrix3& m) {
            for (std::size_t i = 0; i < 9; ++i) {
                data[i] += m.data[i];
            }

            return *this;
        }

        Matrix3& operator-=(const Matrix3& m) {
            for (std::size_t i = 0; i < 9; ++i) {
                data[i] -= m.data[i];
            }

            return *this;
        }

        Matrix3& operator*= (const Matrix3& m) {
            *this = *this * m;
            return *this;
        }

    // =========================
    // Determinant
    // =========================

        double determinant3() const {
            return data[0]*data[4]*data[8] + data[1]*data[5]*data[6] + data[2]*data[3]*data[7] -
                   data[2]*data[4]*data[6] - data[1]*data[3]*data[8] - data[0]*data[5]*data[7];
        }

    // =========================
    // Transpose
    // =========================

        Matrix3 transpose3() const {
            return Matrix3(
                data[0], data[3], data[6],
                data[1], data[4], data[7],
                data[2], data[5], data[8]
            );
        }

        Matrix3& transposeInPlace3() {
            std::swap(data[1], data[3]);
            std::swap(data[2], data[6]);
            std::swap(data[5], data[7]);
            return *this;
        }

    // =========================
    // Cofactor
    // =========================

        Matrix3 cofactorMatrix3() const {
            return Matrix3(
                data[4]*data[8] - data[5]*data[7], -(data[3]*data[8] - data[5]*data[6]), data[3]*data[7] - data[4]*data[6],
                -(data[1]*data[8] - data[2]*data[7]), data[0]*data[8] - data[2]*data[6], -(data[0]*data[7] - data[1]*data[6]),
                data[1]*data[5] - data[2]*data[4], -(data[0]*data[5] - data[2]*data[3]), data[0]*data[4] - data[1]*data[3]
            );
        }

        Matrix3& cofactorInPlace3() {
            *this = cofactorMatrix3();
            return *this;
        }
    // =========================
    // Inverse (general)
    // =========================

        Matrix3 inverse3() const {
            double det = determinant3();
            
            if (std::abs(det) < 1e-9) { 
                throw std::runtime_error("Matrix is singular.");
            }

            Matrix3 adjugate = cofactorMatrix3();
            adjugate.transposeInPlace3();

            adjugate *= (1.0 / det);
            return adjugate;
        }

        Matrix3& invertInPlace3() {
            double det = determinant3();

            if (std::abs(det) < 1e-9) {
                throw std::runtime_error("Matrix is singular and cannot be inverted.");
            }

            cofactorInPlace3();
            transposeInPlace3();
            *this *= (1.0 / det);

            return *this;
        }
};

inline Matrix3 operator*(double scalar, const Matrix3& m) {
    return m * scalar;
}

#endif