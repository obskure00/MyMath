#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <array>
#include <stdexcept>
#include <cmath>
#include <ostream>

class Vector2;
class Vector3;
class Vector4;

//=======================================
//-------------- Matrix2 ----------------
//=======================================

class Matrix2 {
public:
    std::array<double, 4> data;

    Matrix2() { data.fill(0.0); }
    explicit Matrix2(double fill) { data.fill(fill); }

    Matrix2(double m00, double m01,
            double m10, double m11)
        : data{m00, m01,
               m10, m11} {}

    static Matrix2 Identity() {
        return Matrix2(1.0, 0.0,
                       0.0, 1.0);
    }

    // =========================
    // Element Access
    // =========================

    double&       operator()(int row, int col)       { return data[row*2 + col]; }
    const double& operator()(int row, int col) const { return data[row*2 + col]; }

    // =========================
    // Comparison
    // =========================

    bool operator==(const Matrix2& m) const { return data == m.data; }
    bool operator!=(const Matrix2& m) const { return data != m.data; }

    // =========================
    // Unary
    // =========================

    Matrix2 operator-() const {
        Matrix2 r;
        for (int i = 0; i < 4; ++i) r.data[i] = -data[i];
        return r;
    }

    // =========================
    // Arithmetic
    // =========================

    Matrix2 operator+(const Matrix2& m) const {
        Matrix2 r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] + m.data[i];
        return r;
    }

    Matrix2 operator-(const Matrix2& m) const {
        Matrix2 r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] - m.data[i];
        return r;
    }

    Matrix2 operator*(const Matrix2& m) const {
        return Matrix2(
            data[0]*m.data[0] + data[1]*m.data[2],
            data[0]*m.data[1] + data[1]*m.data[3],
            data[2]*m.data[0] + data[3]*m.data[2],
            data[2]*m.data[1] + data[3]*m.data[3]
        );
    }

    Matrix2 operator*(double s) const {
        Matrix2 r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] * s;
        return r;
    }

    Matrix2 operator/(double s) const {
        Matrix2 r;
        for (int i = 0; i < 4; ++i) r.data[i] = data[i] / s;
        return r;
    }

    Matrix2& operator+=(const Matrix2& m)  { for (int i=0;i<4;++i) data[i]+=m.data[i]; return *this; }
    Matrix2& operator-=(const Matrix2& m)  { for (int i=0;i<4;++i) data[i]-=m.data[i]; return *this; }
    Matrix2& operator*=(const Matrix2& m)  { *this = *this * m; return *this; }
    Matrix2& operator*=(double s)          { for (int i=0;i<4;++i) data[i]*=s; return *this; }
    Matrix2& operator/=(double s)          { for (int i=0;i<4;++i) data[i]/=s; return *this; }

    // =========================
    // Determinant
    // =========================

    double determinant2() const {
        return data[0]*data[3] - data[1]*data[2];
    }

    // =========================
    // Transpose
    // =========================

    Matrix2 transpose2() const {
        return Matrix2(data[0], data[2],
                       data[1], data[3]);
    }

    Matrix2& transposeInPlace2() {
        std::swap(data[1], data[2]);
        return *this;
    }

    // =========================
    // Inverse
    // =========================

    Matrix2 inverse2() const {
        double det = determinant2();
        if (std::abs(det) < 1e-9)
            throw std::runtime_error("Matrix2 is singular.");
        return Matrix2( data[3], -data[1],
                       -data[2],  data[0]) / det;
    }

    Matrix2& invertInPlace2() {
        *this = inverse2();
        return *this;
    }

    // =========================
    // Matrix-Vector multiplication — body defined in Math.hpp
    // =========================

    inline Vector2 operator*(const Vector2& v) const;
};

inline Matrix2 operator*(double s, const Matrix2& m) { return m * s; }

inline std::ostream& operator<<(std::ostream& os, const Matrix2& m) {
    return os << "Matrix2([" << m.data[0] << ", " << m.data[1] << "], ["
                             << m.data[2] << ", " << m.data[3] << "])";
}

//=======================================
//-------------- Matrix3 ----------------
//=======================================

class Matrix3 {
public:
    std::array<double, 9> data;

    // =========================
    // Constructors
    // =========================

    Matrix3() { data.fill(0.0); }
    explicit Matrix3(double fill) { data.fill(fill); }

    Matrix3(double m00, double m01, double m02,
            double m10, double m11, double m12,
            double m20, double m21, double m22)
        : data{m00, m01, m02,
               m10, m11, m12,
               m20, m21, m22} {}

    static Matrix3 Identity() {
        Matrix3 m;
        m.data[0] = 1.0; m.data[4] = 1.0; m.data[8] = 1.0;
        return m;
    }

    static Matrix3 fromAxisAngle3(double ax, double ay, double az, double angle) {
        double len = std::sqrt(ax*ax + ay*ay + az*az);
        if (len < 1e-12) return Identity();
        ax /= len; ay /= len; az /= len;
        double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
        return Matrix3(
            t*ax*ax + c,     t*ax*ay - s*az,  t*ax*az + s*ay,
            t*ax*ay + s*az,  t*ay*ay + c,     t*ay*az - s*ax,
            t*ax*az - s*ay,  t*ay*az + s*ax,  t*az*az + c
        );
    }

    inline static Matrix3 fromAxisAngle3(const Vector3& axis, double angle);

    // =========================
    // Element Access
    // =========================

    double&       operator()(int row, int col)       { return data[row*3 + col]; }
    const double& operator()(int row, int col) const { return data[row*3 + col]; }

    // =========================
    // Comparison
    // =========================

    bool operator==(const Matrix3& m) const { return data == m.data; }
    bool operator!=(const Matrix3& m) const { return data != m.data; }

    // =========================
    // Unary
    // =========================

    Matrix3 operator-() const {
        Matrix3 r;
        for (int i = 0; i < 9; ++i) r.data[i] = -data[i];
        return r;
    }

    // =========================
    // Arithmetic
    // =========================

    Matrix3 operator+(const Matrix3& m) const {
        Matrix3 r;
        for (int i = 0; i < 9; ++i) r.data[i] = data[i] + m.data[i];
        return r;
    }

    Matrix3 operator-(const Matrix3& m) const {
        Matrix3 r;
        for (int i = 0; i < 9; ++i) r.data[i] = data[i] - m.data[i];
        return r;
    }

    Matrix3 operator*(const Matrix3& m) const {
        Matrix3 r;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                r.data[row*3+col] =
                    data[row*3+0] * m.data[0*3+col] +
                    data[row*3+1] * m.data[1*3+col] +
                    data[row*3+2] * m.data[2*3+col];
        return r;
    }

    Matrix3 operator*(double s) const {
        Matrix3 r;
        for (int i = 0; i < 9; ++i) r.data[i] = data[i] * s;
        return r;
    }

    Matrix3 operator/(double s) const {
        Matrix3 r;
        for (int i = 0; i < 9; ++i) r.data[i] = data[i] / s;
        return r;
    }

    Matrix3& operator+=(const Matrix3& m)  { for (int i=0;i<9;++i) data[i]+=m.data[i]; return *this; }
    Matrix3& operator-=(const Matrix3& m)  { for (int i=0;i<9;++i) data[i]-=m.data[i]; return *this; }
    Matrix3& operator*=(const Matrix3& m)  { *this = *this * m; return *this; }
    Matrix3& operator*=(double s)          { for (int i=0;i<9;++i) data[i]*=s; return *this; }
    Matrix3& operator/=(double s)          { for (int i=0;i<9;++i) data[i]/=s; return *this; }

    // =========================
    // Determinant
    // =========================

    double determinant3() const {
        return data[0]*(data[4]*data[8] - data[5]*data[7])
             - data[1]*(data[3]*data[8] - data[5]*data[6])
             + data[2]*(data[3]*data[7] - data[4]*data[6]);
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
    // Cofactor matrix
    // =========================

    Matrix3 cofactorMatrix3() const {
        return Matrix3(
             data[4]*data[8] - data[5]*data[7],
            -(data[3]*data[8] - data[5]*data[6]),
             data[3]*data[7] - data[4]*data[6],

            -(data[1]*data[8] - data[2]*data[7]),
             data[0]*data[8] - data[2]*data[6],
            -(data[0]*data[7] - data[1]*data[6]),

             data[1]*data[5] - data[2]*data[4],
            -(data[0]*data[5] - data[2]*data[3]),
             data[0]*data[4] - data[1]*data[3]
        );
    }

    Matrix3& cofactorInPlace3() {
        *this = cofactorMatrix3();
        return *this;
    }

    // =========================
    // Inverse
    // =========================

    Matrix3 inverse3() const {
        double det = determinant3();
        if (std::abs(det) < 1e-9)
            throw std::runtime_error("Matrix3 is singular.");
        Matrix3 adj = cofactorMatrix3();
        adj.transposeInPlace3();
        adj *= (1.0 / det);
        return adj;
    }

    Matrix3& invertInPlace3() {
        *this = inverse3();
        return *this;
    }

    // =========================
    // Inertia Tensor Transform
    //
    // Transforms a body-space inertia tensor to world space:
    //   I_world = R * I_body * R^T
    // where R maps body-space directions to world-space directions.
    // =========================

    static Matrix3 transformInertiaTensor3(const Matrix3& inertiaBody, const Matrix3& rotation) {
        return rotation * inertiaBody * rotation.transpose3();
    }

    void transformInertiaTensorInPlace3(const Matrix3& rotation) {
        *this = rotation * (*this) * rotation.transpose3();
    }

    // =========================
    // Symmetric Inertia Tensor Transform  (optimised)
    //
    // Inertia tensors are always symmetric: I[i][j] == I[j][i].
    // This overload exploits that fact in both the input and the output,
    // reading only 6 unique elements and writing only 6 unique elements,
    // saving ~17% multiplications vs the general path.
    //
    // ixx, ixy, ixz, iyy, iyz, izz are the upper-triangle elements:
    //   I = | ixx  ixy  ixz |
    //       | ixy  iyy  iyz |
    //       | ixz  iyz  izz |
    //
    // R is the same rotation matrix as in transformInertiaTensor3.
    // =========================

    static Matrix3 transformSymmetricInertiaTensor3(
        double ixx, double ixy, double ixz,
                    double iyy, double iyz,
                                double izz,
        const Matrix3& R)
    {
        // Name the rows of R for readability
        const double r00 = R.data[0], r01 = R.data[1], r02 = R.data[2];
        const double r10 = R.data[3], r11 = R.data[4], r12 = R.data[5];
        const double r20 = R.data[6], r21 = R.data[7], r22 = R.data[8];

        // Step 1: T = R * I  (27 mults — rows of R dotted with cols of I,
        // but cols of symmetric I == rows of I, so we reuse the same 6 values)
        const double t00 = r00*ixx + r01*ixy + r02*ixz;
        const double t01 = r00*ixy + r01*iyy + r02*iyz;
        const double t02 = r00*ixz + r01*iyz + r02*izz;

        const double t10 = r10*ixx + r11*ixy + r12*ixz;
        const double t11 = r10*ixy + r11*iyy + r12*iyz;
        const double t12 = r10*ixz + r11*iyz + r12*izz;

        const double t20 = r20*ixx + r21*ixy + r22*ixz;
        const double t21 = r20*ixy + r21*iyy + r22*iyz;
        const double t22 = r20*ixz + r21*iyz + r22*izz;

        // Step 2: result = T * R^T  — only compute the 6 unique elements
        // of the symmetric result matrix.
        // result[i][j] = T[i] · R[j]   (dot row i of T with row j of R)
        const double rxx = t00*r00 + t01*r01 + t02*r02;  // [0][0]
        const double rxy = t00*r10 + t01*r11 + t02*r12;  // [0][1] == [1][0]
        const double rxz = t00*r20 + t01*r21 + t02*r22;  // [0][2] == [2][0]
        const double ryy = t10*r10 + t11*r11 + t12*r12;  // [1][1]
        const double ryz = t10*r20 + t11*r21 + t12*r22;  // [1][2] == [2][1]
        const double rzz = t20*r20 + t21*r21 + t22*r22;  // [2][2]

        return Matrix3(
            rxx, rxy, rxz,
            rxy, ryy, ryz,
            rxz, ryz, rzz
        );
    }

    // Overload that reads the symmetric tensor from the upper triangle of a Matrix3.
    // Caller is responsible for ensuring the matrix is actually symmetric.
    static Matrix3 transformSymmetricInertiaTensor3(const Matrix3& I, const Matrix3& R) {
        return transformSymmetricInertiaTensor3(
            I.data[0], I.data[1], I.data[2],
                       I.data[4], I.data[5],
                                  I.data[8],
            R
        );
    }

    void transformSymmetricInertiaTensorInPlace3(const Matrix3& R) {
        *this = transformSymmetricInertiaTensor3(*this, R);
    }

    // =========================
    // Matrix-Vector multiplication — body defined in Math.hpp
    // =========================

    inline Vector3 operator*(const Vector3& v) const;
};

inline Matrix3 operator*(double s, const Matrix3& m) { return m * s; }

inline std::ostream& operator<<(std::ostream& os, const Matrix3& m) {
    return os << "Matrix3(["
              << m.data[0] << ", " << m.data[1] << ", " << m.data[2] << "], ["
              << m.data[3] << ", " << m.data[4] << ", " << m.data[5] << "], ["
              << m.data[6] << ", " << m.data[7] << ", " << m.data[8] << "])";
}

// =========================
// Axis-aligned rotation matrices
// =========================

inline Matrix3 rotationX3(double rad) {
    return Matrix3(
        1.0,  0.0,       0.0,
        0.0,  std::cos(rad), -std::sin(rad),
        0.0,  std::sin(rad),  std::cos(rad)
    );
}

inline Matrix3 rotationY3(double rad) {
    return Matrix3(
         std::cos(rad), 0.0, std::sin(rad),
         0.0,           1.0, 0.0,
        -std::sin(rad), 0.0, std::cos(rad)
    );
}

inline Matrix3 rotationZ3(double rad) {
    return Matrix3(
        std::cos(rad), -std::sin(rad), 0.0,
        std::sin(rad),  std::cos(rad), 0.0,
        0.0,            0.0,           1.0
    );
}

//=======================================
//-------------- Matrix4 ----------------
//=======================================

class Matrix4 {
public:
    std::array<double, 16> data;

    // =========================
    // Constructors
    // =========================

    Matrix4() { data.fill(0.0); }
    explicit Matrix4(double fill) { data.fill(fill); }

    Matrix4(double m00, double m01, double m02, double m03,
            double m10, double m11, double m12, double m13,
            double m20, double m21, double m22, double m23,
            double m30, double m31, double m32, double m33)
        : data{m00, m01, m02, m03,
               m10, m11, m12, m13,
               m20, m21, m22, m23,
               m30, m31, m32, m33} {}

    static Matrix4 Identity() {
        Matrix4 m;
        m.data[0] = 1.0; m.data[5] = 1.0; m.data[10] = 1.0; m.data[15] = 1.0;
        return m;
    }

    static Matrix4 fromRotation3(const Matrix3& r) {
        Matrix4 m = Identity();
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                m(row, col) = r(row, col);
        return m;
    }

    inline static Matrix4 fromTranslation3(const Vector3& t);
    inline static Matrix4 fromScale3(const Vector3& s);
    inline static Matrix4 composeTRS3(const Vector3& translation,
                                      const Matrix3& rotation,
                                      const Vector3& scale);

    // =========================
    // Element Access
    // =========================

    double&       operator()(int row, int col)       { return data[row*4 + col]; }
    const double& operator()(int row, int col) const { return data[row*4 + col]; }

    // =========================
    // Comparison
    // =========================

    bool operator==(const Matrix4& m) const { return data == m.data; }
    bool operator!=(const Matrix4& m) const { return data != m.data; }

    // =========================
    // Unary
    // =========================

    Matrix4 operator-() const {
        Matrix4 r;
        for (int i = 0; i < 16; ++i) r.data[i] = -data[i];
        return r;
    }

    // =========================
    // Arithmetic
    // =========================

    Matrix4 operator+(const Matrix4& m) const {
        Matrix4 r;
        for (int i = 0; i < 16; ++i) r.data[i] = data[i] + m.data[i];
        return r;
    }

    Matrix4 operator-(const Matrix4& m) const {
        Matrix4 r;
        for (int i = 0; i < 16; ++i) r.data[i] = data[i] - m.data[i];
        return r;
    }

    Matrix4 operator*(const Matrix4& m) const {
        Matrix4 r;
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                r.data[row*4+col] =
                    data[row*4+0] * m.data[0*4+col] +
                    data[row*4+1] * m.data[1*4+col] +
                    data[row*4+2] * m.data[2*4+col] +
                    data[row*4+3] * m.data[3*4+col];
        return r;
    }

    Matrix4 operator*(double s) const {
        Matrix4 r;
        for (int i = 0; i < 16; ++i) r.data[i] = data[i] * s;
        return r;
    }

    Matrix4 operator/(double s) const {
        Matrix4 r;
        for (int i = 0; i < 16; ++i) r.data[i] = data[i] / s;
        return r;
    }

    Matrix4& operator+=(const Matrix4& m)  { for (int i=0;i<16;++i) data[i]+=m.data[i]; return *this; }
    Matrix4& operator-=(const Matrix4& m)  { for (int i=0;i<16;++i) data[i]-=m.data[i]; return *this; }
    Matrix4& operator*=(const Matrix4& m)  { *this = *this * m; return *this; }
    Matrix4& operator*=(double s)          { for (int i=0;i<16;++i) data[i]*=s; return *this; }
    Matrix4& operator/=(double s)          { for (int i=0;i<16;++i) data[i]/=s; return *this; }

    // =========================
    // Transpose
    // =========================

    Matrix4 transpose4() const {
        Matrix4 r;
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                r.data[col*4+row] = data[row*4+col];
        return r;
    }

    Matrix4& transposeInPlace4() {
        for (int row = 0; row < 4; ++row)
            for (int col = row+1; col < 4; ++col)
                std::swap(data[row*4+col], data[col*4+row]);
        return *this;
    }

    // =========================
    // Determinant
    // =========================

    double determinant4() const {
        const auto& d = data;
        double c0 =  det3(d[5],d[6],d[7],  d[9],d[10],d[11],  d[13],d[14],d[15]);
        double c1 = -det3(d[4],d[6],d[7],  d[8],d[10],d[11],  d[12],d[14],d[15]);
        double c2 =  det3(d[4],d[5],d[7],  d[8],d[9], d[11],  d[12],d[13],d[15]);
        double c3 = -det3(d[4],d[5],d[6],  d[8],d[9], d[10],  d[12],d[13],d[14]);
        return d[0]*c0 + d[1]*c1 + d[2]*c2 + d[3]*c3;
    }

    // =========================
    // General Inverse
    // Throws if matrix is singular
    // =========================

    Matrix4 inverse4() const {
        const auto& d = data;
        double c[16];
        c[ 0] =  det3(d[5],d[6],d[7],  d[9],d[10],d[11],  d[13],d[14],d[15]);
        c[ 1] = -det3(d[4],d[6],d[7],  d[8],d[10],d[11],  d[12],d[14],d[15]);
        c[ 2] =  det3(d[4],d[5],d[7],  d[8],d[9], d[11],  d[12],d[13],d[15]);
        c[ 3] = -det3(d[4],d[5],d[6],  d[8],d[9], d[10],  d[12],d[13],d[14]);

        c[ 4] = -det3(d[1],d[2],d[3],  d[9],d[10],d[11],  d[13],d[14],d[15]);
        c[ 5] =  det3(d[0],d[2],d[3],  d[8],d[10],d[11],  d[12],d[14],d[15]);
        c[ 6] = -det3(d[0],d[1],d[3],  d[8],d[9], d[11],  d[12],d[13],d[15]);
        c[ 7] =  det3(d[0],d[1],d[2],  d[8],d[9], d[10],  d[12],d[13],d[14]);

        c[ 8] =  det3(d[1],d[2],d[3],  d[5],d[6],d[7],    d[13],d[14],d[15]);
        c[ 9] = -det3(d[0],d[2],d[3],  d[4],d[6],d[7],    d[12],d[14],d[15]);
        c[10] =  det3(d[0],d[1],d[3],  d[4],d[5],d[7],    d[12],d[13],d[15]);
        c[11] = -det3(d[0],d[1],d[2],  d[4],d[5],d[6],    d[12],d[13],d[14]);

        c[12] = -det3(d[1],d[2],d[3],  d[5],d[6],d[7],    d[9],d[10],d[11]);
        c[13] =  det3(d[0],d[2],d[3],  d[4],d[6],d[7],    d[8],d[10],d[11]);
        c[14] = -det3(d[0],d[1],d[3],  d[4],d[5],d[7],    d[8],d[9], d[11]);
        c[15] =  det3(d[0],d[1],d[2],  d[4],d[5],d[6],    d[8],d[9], d[10]);

        double det = d[0]*c[0] + d[1]*c[1] + d[2]*c[2] + d[3]*c[3];
        if (std::abs(det) < 1e-12)
            throw std::runtime_error("Matrix4 is singular.");

        double inv = 1.0 / det;

        Matrix4 result;
        for (int row = 0; row < 4; ++row)
            for (int col = 0; col < 4; ++col)
                result.data[row*4+col] = c[col*4+row] * inv;
        return result;
    }

    Matrix4& invertInPlace4() {
        *this = inverse4();
        return *this;
    }

    // =========================
    // Fast Rigid-Body Inverse
    // faster than the general inverse
    // =========================

    Matrix4 inverseFastRigid4() const {
        Matrix4 r;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                r(row, col) = (*this)(col, row);
        double tx = data[3], ty = data[7], tz = data[11];
        r(0,3) = -(r(0,0)*tx + r(0,1)*ty + r(0,2)*tz);
        r(1,3) = -(r(1,0)*tx + r(1,1)*ty + r(1,2)*tz);
        r(2,3) = -(r(2,0)*tx + r(2,1)*ty + r(2,2)*tz);
        r(3,3) = 1.0;
        return r;
    }

    // =========================
    // Transform helpers — bodies defined in Math.hpp
    // =========================

    inline Vector3 transformPoint3(const Vector3& p) const;

    inline Vector3 transformDirection3(const Vector3& d) const;

    // =========================
    // Matrix-Vector multiplication — body defined in Math.hpp
    // =========================

    inline Vector4 operator*(const Vector4& v) const;

private:
    static double det3(double a00, double a01, double a02,
                       double a10, double a11, double a12,
                       double a20, double a21, double a22)
    {
        return a00*(a11*a22 - a12*a21)
             - a01*(a10*a22 - a12*a20)
             + a02*(a10*a21 - a11*a20);
    }
};

inline Matrix4 operator*(double s, const Matrix4& m) { return m * s; }

inline std::ostream& operator<<(std::ostream& os, const Matrix4& m) {
    return os << "Matrix4(["
              << m.data[0]  << ", " << m.data[1]  << ", " << m.data[2]  << ", " << m.data[3]  << "], ["
              << m.data[4]  << ", " << m.data[5]  << ", " << m.data[6]  << ", " << m.data[7]  << "], ["
              << m.data[8]  << ", " << m.data[9]  << ", " << m.data[10] << ", " << m.data[11] << "], ["
              << m.data[12] << ", " << m.data[13] << ", " << m.data[14] << ", " << m.data[15] << "])";
}

#endif