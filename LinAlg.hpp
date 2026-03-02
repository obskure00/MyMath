#ifndef LINALG_HPP
#define LINALG_HPP

// LinAlg.hpp — advanced linear algebra utilities for physics solvers.
// All declarations live in namespace LinAlg to avoid name pollution.
// Include via Math.hpp directly.

#include "Math.hpp"

#include <array>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <ostream>

namespace LinAlg {

// ============================================================
//  SMALL MATRIX SOLVERS
// ============================================================

// ---- solve2 ------------------------------------------------

inline Vector2 solve2(const Matrix2& A, const Vector2& b) {
    double det = A.determinant2();
    if (std::fabs(det) < 1e-12)
        throw std::runtime_error("solve2: matrix is singular.");

    double invDet = 1.0 / det;
    return Vector2(
        (b.x * A.data[3] - b.y * A.data[1]) * invDet,
        (A.data[0] * b.y - A.data[2] * b.x) * invDet
    );
}

// ---- solve3 ------------------------------------------------

inline Vector3 solve3(const Matrix3& A, const Vector3& b) {
    double det = A.determinant3();
    if (std::fabs(det) < 1e-12)
        throw std::runtime_error("solve3: matrix is singular.");

    const double* d = A.data.data();
    double invDet   = 1.0 / det;

    // det(A with column 0 replaced by b)
    double dx =  b.x * (d[4]*d[8] - d[5]*d[7])
              - d[1] * (b.y*d[8]  - d[5]*b.z)
              + d[2] * (b.y*d[7]  - d[4]*b.z);

    // det(A with column 1 replaced by b)
    double dy = d[0] * (b.y*d[8]  - d[5]*b.z)
              -  b.x * (d[3]*d[8] - d[5]*d[6])
              + d[2] * (d[3]*b.z  - b.y*d[6]);

    // det(A with column 2 replaced by b)
    double dz = d[0] * (d[4]*b.z  - b.y*d[7])
              - d[1] * (d[3]*b.z  - b.y*d[6])
              +  b.x * (d[3]*d[7] - d[4]*d[6]);

    return Vector3(dx * invDet, dy * invDet, dz * invDet);
}

// ============================================================
//  Matrix6
//  6x6 row-major matrix
// ============================================================

class Matrix6 {
    public:
        std::array<double, 36> data;

        Matrix6()                  { data.fill(0.0); }
        explicit Matrix6(double f) { data.fill(f); }

        static Matrix6 Identity() {
            Matrix6 m;
            for (int i = 0; i < 6; ++i) m.data[static_cast<size_t>(i*6+i)] = 1.0;
            return m;
        }

        // -------- Element access (row-major) --------
        double&       operator()(int r, int c)       { return data[static_cast<size_t>(r*6+c)]; }
        const double& operator()(int r, int c) const { return data[static_cast<size_t>(r*6+c)]; }

        // -------- Comparison --------
        bool operator==(const Matrix6& m) const { return data == m.data; }
        bool operator!=(const Matrix6& m) const { return data != m.data; }

        // -------- Unary --------
        Matrix6 operator-() const {
            Matrix6 r;
            for (int i = 0; i < 36; ++i) r.data[static_cast<size_t>(i)] = -data[static_cast<size_t>(i)];
            return r;
        }

        // -------- Arithmetic --------
        Matrix6 operator+(const Matrix6& m) const {
            Matrix6 r;
            for (int i = 0; i < 36; ++i) r.data[static_cast<size_t>(i)] = data[static_cast<size_t>(i)] + m.data[static_cast<size_t>(i)];
            return r;
        }
        Matrix6 operator-(const Matrix6& m) const {
            Matrix6 r;
            for (int i = 0; i < 36; ++i) r.data[static_cast<size_t>(i)] = data[static_cast<size_t>(i)] - m.data[static_cast<size_t>(i)];
            return r;
        }
        Matrix6 operator*(const Matrix6& m) const {
            Matrix6 r;
            for (int row = 0; row < 6; ++row)
                for (int col = 0; col < 6; ++col) {
                    double s = 0.0;
                    for (int k = 0; k < 6; ++k)
                        s += data[static_cast<size_t>(row*6+k)] * m.data[static_cast<size_t>(k*6+col)];
                    r.data[static_cast<size_t>(row*6+col)] = s;
                }
            return r;
        }
        Matrix6 operator*(double s) const {
            Matrix6 r;
            for (int i = 0; i < 36; ++i) r.data[static_cast<size_t>(i)] = data[static_cast<size_t>(i)] * s;
            return r;
        }
        Matrix6 operator/(double s) const {
            Matrix6 r;
            for (int i = 0; i < 36; ++i) r.data[static_cast<size_t>(i)] = data[static_cast<size_t>(i)] / s;
            return r;
        }

        Matrix6& operator+=(const Matrix6& m) { for(int i=0;i<36;++i) data[static_cast<size_t>(i)]+=m.data[static_cast<size_t>(i)]; return *this; }
        Matrix6& operator-=(const Matrix6& m) { for(int i=0;i<36;++i) data[static_cast<size_t>(i)]-=m.data[static_cast<size_t>(i)]; return *this; }
        Matrix6& operator*=(const Matrix6& m) { *this = *this * m; return *this; }
        Matrix6& operator*=(double s)         { for(int i=0;i<36;++i) data[static_cast<size_t>(i)]*=s; return *this; }
        Matrix6& operator/=(double s)         { for(int i=0;i<36;++i) data[static_cast<size_t>(i)]/=s; return *this; }

        Matrix6 transpose6() const {
            Matrix6 r;
            for (int row = 0; row < 6; ++row)
                for (int col = 0; col < 6; ++col)
                    r.data[static_cast<size_t>(col*6+row)] = data[static_cast<size_t>(row*6+col)];
            return r;
        }

        std::array<double, 6> multiplyVec(const std::array<double, 6>& v) const {
            std::array<double, 6> r{};
            for (int row = 0; row < 6; ++row)
                for (int col = 0; col < 6; ++col)
                    r[static_cast<size_t>(row)] += data[static_cast<size_t>(row*6+col)] * v[static_cast<size_t>(col)];
            return r;
        }

        static Matrix6 inverseMassMatrix(double invMass, const Matrix3& IwInv) {
            Matrix6 m;
            m(0,0) = invMass;  m(1,1) = invMass;  m(2,2) = invMass;
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    m(r+3, c+3) = IwInv(r, c);
            return m;
        }
};

inline Matrix6 operator*(double s, const Matrix6& m) { return m * s; }

inline std::ostream& operator<<(std::ostream& os, const Matrix6& m) {
    os << "Matrix6(\n";
    for (int row = 0; row < 6; ++row) {
        os << "  [";
        for (int col = 0; col < 6; ++col) {
            os << m.data[static_cast<size_t>(row*6+col)];
            if (col < 5) os << ", ";
        }
        os << "]\n";
    }
    return os << ")";
}

// ============================================================
//  solve6  (convenience wrapper around LUDecomp<6>)
// ============================================================

inline std::array<double, 6> solve6(const Matrix6& A,
                                    const std::array<double, 6>& b);

// ---- solveBlock3x3 -----------------------------------------

inline void solveBlock3x3(const Matrix3& A, const Matrix3& B,
                          const Matrix3& C, const Matrix3& D,
                          const Vector3& b1, const Vector3& b2,
                          Vector3& x1, Vector3& x2) {
    Matrix3 Ainv = A.inverse3();
    Matrix3 S    = D - C * Ainv * B;
    Vector3 rhs2 = b2 - C * (Ainv * b1);
    x2 = solve3(S, rhs2);
    x1 = Ainv * (b1 - B * x2);
}

// ============================================================
//  DECOMPOSITIONS
// ============================================================

// ============================================================
//  LUDecomp<N>
// ============================================================

template<int N>
struct LUDecomp {
    static_assert(N > 0, "LUDecomp: N must be positive");

    std::array<double, N*N> LU;
    std::array<int,    N>   piv;
    bool                    singular = false;

    static LUDecomp factorize(const std::array<double, N*N>& A) {
        LUDecomp lu;
        lu.LU  = A;
        lu.singular = false;
        for (int i = 0; i < N; ++i) lu.piv[static_cast<size_t>(i)] = i;

        for (int k = 0; k < N; ++k) {
            int    pivot_row = k;
            double pivot_val = std::fabs(lu.LU[static_cast<size_t>(k*N + k)]);
            for (int i = k+1; i < N; ++i) {
                double v = std::fabs(lu.LU[static_cast<size_t>(i*N + k)]);
                if (v > pivot_val) { pivot_val = v; pivot_row = i; }
            }

            if (pivot_row != k) {
                std::swap(lu.piv[static_cast<size_t>(k)],
                          lu.piv[static_cast<size_t>(pivot_row)]);
                for (int j = 0; j < N; ++j)
                    std::swap(lu.LU[static_cast<size_t>(k*N+j)],
                              lu.LU[static_cast<size_t>(pivot_row*N+j)]);
            }

            double pivot = lu.LU[static_cast<size_t>(k*N + k)];
            if (std::fabs(pivot) < 1e-15) {
                lu.singular = true;
                return lu;
            }

            double invPivot = 1.0 / pivot;
            for (int i = k+1; i < N; ++i) {
                lu.LU[static_cast<size_t>(i*N + k)] *= invPivot;
                for (int j = k+1; j < N; ++j)
                    lu.LU[static_cast<size_t>(i*N + j)] -=
                        lu.LU[static_cast<size_t>(i*N + k)] *
                        lu.LU[static_cast<size_t>(k*N + j)];
            }
        }
        return lu;
    }

    std::array<double, N> solve(const std::array<double, N>& b) const {
        if (singular)
            throw std::runtime_error("LUDecomp::solve: singular matrix.");

        std::array<double, N> x;
        for (int i = 0; i < N; ++i)
            x[static_cast<size_t>(i)] = b[static_cast<size_t>(piv[static_cast<size_t>(i)])];

        for (int i = 1; i < N; ++i)
            for (int j = 0; j < i; ++j)
                x[static_cast<size_t>(i)] -= LU[static_cast<size_t>(i*N + j)]
                                           * x[static_cast<size_t>(j)];

        for (int i = N-1; i >= 0; --i) {
            for (int j = i+1; j < N; ++j)
                x[static_cast<size_t>(i)] -= LU[static_cast<size_t>(i*N + j)]
                                           * x[static_cast<size_t>(j)];
            x[static_cast<size_t>(i)] /= LU[static_cast<size_t>(i*N + i)];
        }
        return x;
    }

    double determinant() const {
        if (singular) return 0.0;

        double det = 1.0;
        for (int i = 0; i < N; ++i) {
            det *= LU[static_cast<size_t>(i*N + i)];
        }

        std::array<int, N> p = piv;
        for (int i = 0; i < N; ++i) {
            while (p[static_cast<size_t>(i)] != i) {
                int j = p[static_cast<size_t>(i)];
                std::swap(p[static_cast<size_t>(i)], p[static_cast<size_t>(j)]);
                det = -det;
            }
        }
        return det;
    }
};

// ---- solve6 body (after LUDecomp) --------------------------

inline std::array<double, 6> solve6(const Matrix6& A,
                                    const std::array<double, 6>& b) {
    auto lu = LUDecomp<6>::factorize(A.data);
    if (lu.singular)
        throw std::runtime_error("solve6: matrix is singular.");
    return lu.solve(b);
}

// ============================================================
//  CholeskyDecomp<N>
// ============================================================

template<int N>
struct CholeskyDecomp {
    static_assert(N > 0, "CholeskyDecomp: N must be positive");

    std::array<double, N*N> L;
    bool valid = false;

    static CholeskyDecomp factorize(const std::array<double, N*N>& A) {
        CholeskyDecomp ch;
        ch.L.fill(0.0);

        for (int j = 0; j < N; ++j) {
            double diag = A[static_cast<size_t>(j*N + j)];
            for (int k = 0; k < j; ++k)
                diag -= ch.L[static_cast<size_t>(j*N + k)]
                      * ch.L[static_cast<size_t>(j*N + k)];

            if (diag <= 0.0) {
                ch.valid = false;
                return ch;
            }
            ch.L[static_cast<size_t>(j*N + j)] = std::sqrt(diag);

            double invLjj = 1.0 / ch.L[static_cast<size_t>(j*N + j)];

            for (int i = j+1; i < N; ++i) {
                double s = A[static_cast<size_t>(i*N + j)];
                for (int k = 0; k < j; ++k)
                    s -= ch.L[static_cast<size_t>(i*N + k)]
                       * ch.L[static_cast<size_t>(j*N + k)];
                ch.L[static_cast<size_t>(i*N + j)] = s * invLjj;
            }
        }
        ch.valid = true;
        return ch;
    }

    std::array<double, N> solve(const std::array<double, N>& b) const {
        if (!valid)
            throw std::runtime_error(
                "CholeskyDecomp::solve: decomposition failed — matrix not SPD.");

        std::array<double, N> y{}, x{};

        for (int i = 0; i < N; ++i) {
            double s = b[static_cast<size_t>(i)];
            for (int j = 0; j < i; ++j)
                s -= L[static_cast<size_t>(i*N + j)] * y[static_cast<size_t>(j)];
            y[static_cast<size_t>(i)] = s / L[static_cast<size_t>(i*N + i)];
        }

        for (int i = N-1; i >= 0; --i) {
            double s = y[static_cast<size_t>(i)];
            for (int j = i+1; j < N; ++j)
                s -= L[static_cast<size_t>(j*N + i)] * x[static_cast<size_t>(j)];
            x[static_cast<size_t>(i)] = s / L[static_cast<size_t>(i*N + i)];
        }
        return x;
    }

    double determinant() const {
        if (!valid) return 0.0;
        double d = 1.0;
        for (int i = 0; i < N; ++i)
            d *= L[static_cast<size_t>(i*N + i)];
        return d * d;
    }
};

// ============================================================
//  QRDecomp<N>
// ============================================================

template<int N>
struct QRDecomp {
    static_assert(N > 0, "QRDecomp: N must be positive");

    std::array<double, N*N> Q;
    std::array<double, N*N> R;
    bool singular = false;

    static QRDecomp factorize(const std::array<double, N*N>& A) {
        QRDecomp qr;
        qr.R = A;
        qr.Q.fill(0.0);
        for (int i = 0; i < N; ++i)
            qr.Q[static_cast<size_t>(i*N + i)] = 1.0;

        for (int k = 0; k < N-1; ++k) {
            std::array<double, N> v{};
            double norm2 = 0.0;
            for (int i = k; i < N; ++i) {
                v[static_cast<size_t>(i)] = qr.R[static_cast<size_t>(i*N + k)];
                norm2 += v[static_cast<size_t>(i)] * v[static_cast<size_t>(i)];
            }
            double norm = std::sqrt(norm2);
            if (norm < 1e-15) continue;

            if (v[static_cast<size_t>(k)] >= 0.0)
                v[static_cast<size_t>(k)] += norm;
            else
                v[static_cast<size_t>(k)] -= norm;

            double vTv = 0.0;
            for (int i = k; i < N; ++i)
                vTv += v[static_cast<size_t>(i)] * v[static_cast<size_t>(i)];
            if (vTv < 1e-30) continue;

            double beta = 2.0 / vTv;

            for (int j = k; j < N; ++j) {
                double vTRj = 0.0;
                for (int i = k; i < N; ++i)
                    vTRj += v[static_cast<size_t>(i)] * qr.R[static_cast<size_t>(i*N + j)];
                for (int i = k; i < N; ++i)
                    qr.R[static_cast<size_t>(i*N + j)] -= beta * v[static_cast<size_t>(i)] * vTRj;
            }

            for (int i = 0; i < N; ++i) {
                double Qv = 0.0;
                for (int l = k; l < N; ++l)
                    Qv += qr.Q[static_cast<size_t>(i*N + l)] * v[static_cast<size_t>(l)];
                for (int l = k; l < N; ++l)
                    qr.Q[static_cast<size_t>(i*N + l)] -= beta * Qv * v[static_cast<size_t>(l)];
            }
        }

        qr.singular = false;
        for (int i = 0; i < N; ++i)
            if (std::fabs(qr.R[static_cast<size_t>(i*N + i)]) < 1e-12)
                qr.singular = true;

        return qr;
    }

    std::array<double, N> solve(const std::array<double, N>& b) const {
        if (singular)
            throw std::runtime_error("QRDecomp::solve: matrix is rank-deficient.");

        std::array<double, N> rhs{};
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                rhs[static_cast<size_t>(i)] +=
                    Q[static_cast<size_t>(j*N + i)] * b[static_cast<size_t>(j)];

        std::array<double, N> x{};
        for (int i = N-1; i >= 0; --i) {
            double s = rhs[static_cast<size_t>(i)];
            for (int j = i+1; j < N; ++j)
                s -= R[static_cast<size_t>(i*N + j)] * x[static_cast<size_t>(j)];
            x[static_cast<size_t>(i)] = s / R[static_cast<size_t>(i*N + i)];
        }
        return x;
    }
};

// ============================================================
//  jacobiEigen3  — symmetric 3x3 eigendecomposition
//
//  Finds all three eigenvalues and eigenvectors of a real symmetric
//  3x3 matrix via cyclic Jacobi sweeps.
// ============================================================

struct JacobiResult {
    Vector3 eigenvalues;
    Matrix3 eigenvectors;
    bool    converged;
};

inline JacobiResult jacobiEigen3(const Matrix3& sym, int maxIter = 100) {
    double A[3][3] = {
        { sym.data[0], sym.data[1], sym.data[2] },
        { sym.data[3], sym.data[4], sym.data[5] },
        { sym.data[6], sym.data[7], sym.data[8] }
    };
    double V[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };

    bool converged = false;

    for (int iter = 0; iter < maxIter; ++iter) {
        int p = 0, q = 1;
        double maxOff = std::fabs(A[0][1]);
        if (std::fabs(A[0][2]) > maxOff) { maxOff = std::fabs(A[0][2]); p=0; q=2; }
        if (std::fabs(A[1][2]) > maxOff) { maxOff = std::fabs(A[1][2]); p=1; q=2; }

        if (maxOff < 1e-12) { converged = true; break; }

        double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t;
        if (theta >= 0.0)
            t =  1.0 / ( theta + std::sqrt(1.0 + theta*theta));
        else
            t = -1.0 / (-theta + std::sqrt(1.0 + theta*theta));

        double c   = 1.0 / std::sqrt(1.0 + t*t);
        double s   = t * c;
        double tau = s / (1.0 + c);

        double Apq = A[p][q];
        A[p][q] = 0.0;  A[q][p] = 0.0;
        A[p][p] -= t * Apq;
        A[q][q] += t * Apq;

        for (int r = 0; r < 3; ++r) {
            if (r == p || r == q) continue;
            double Apr = A[p][r];
            double Aqr = A[q][r];
            A[p][r] = A[r][p] = Apr - s * (Aqr + tau * Apr);
            A[q][r] = A[r][q] = Aqr + s * (Apr - tau * Aqr);
        }

        for (int r = 0; r < 3; ++r) {
            double Vrp = V[r][p];
            double Vrq = V[r][q];
            V[r][p] = Vrp - s * (Vrq + tau * Vrp);
            V[r][q] = Vrq + s * (Vrp - tau * Vrq);
        }
    }

    double evals[3] = { A[0][0], A[1][1], A[2][2] };

    for (int i = 0; i < 2; ++i) {
        for (int j = i+1; j < 3; ++j) {
            if (evals[j] > evals[i]) {
                std::swap(evals[i], evals[j]);
                for (int r = 0; r < 3; ++r)
                    std::swap(V[r][i], V[r][j]);
            }
        }
    }
    JacobiResult result;
    result.converged   = converged;
    result.eigenvalues = Vector3(evals[0], evals[1], evals[2]);
    result.eigenvectors = Matrix3(
        V[0][0], V[0][1], V[0][2],
        V[1][0], V[1][1], V[1][2],
        V[2][0], V[2][1], V[2][2]
    );
    return result;
}

// ============================================================
//  10.3  CONSTRAINT SOLVER SUPPORT
// ============================================================

// ============================================================
//  ConstraintJacobian
//
//  One row of the constraint Jacobian for a two-body constraint.
//  The 1x12 row is split into four 3-vectors:
//
//    J * [vA; omegaA; vB; omegaB]
//      =  linearA  · vA
//      +  angularA · omegaA
//      +  linearB  · vB
//      +  angularB · omegaB
// ============================================================

struct ConstraintJacobian {
    Vector3 linearA;
    Vector3 angularA;
    Vector3 linearB;
    Vector3 angularB;

    ConstraintJacobian() = default;

    ConstraintJacobian(const Vector3& lA, const Vector3& aA,
                       const Vector3& lB, const Vector3& aB)
        : linearA(lA), angularA(aA), linearB(lB), angularB(aB) {}

    static ConstraintJacobian contact(const Vector3& n,
                                      const Vector3& rA, const Vector3& rB) {
        return { n, rA.cross3(n), -n, -rB.cross3(n) };
    }

    static ConstraintJacobian friction(const Vector3& t,
                                      const Vector3& rA, const Vector3& rB) {
        return { t, rA.cross3(t), -t, -rB.cross3(t) };
    }
};

// ============================================================
//  constraintVelocity
//  Computes  J * v  — the scalar relative velocity along the constraint.
//  A negative value means the constraint is being compressed/violated.
//  Used as the input "bias" when computing the impulse magnitude.
// ============================================================

inline double constraintVelocity(const ConstraintJacobian& J,
                                 const Vector3& velA,  const Vector3& omegaA,
                                 const Vector3& velB,  const Vector3& omegaB) {
    return J.linearA.dot3(velA)   + J.angularA.dot3(omegaA)
         + J.linearB.dot3(velB)   + J.angularB.dot3(omegaB);
}

// ============================================================
//  effectiveMass
// ============================================================

inline double effectiveMass(const ConstraintJacobian& J,
                            double invMassA, const Matrix3& IwA_inv,
                            double invMassB, const Matrix3& IwB_inv) {
    return  J.linearA.lengthSquared3()  * invMassA
          + J.angularA.dot3(IwA_inv * J.angularA)
          + J.linearB.lengthSquared3()  * invMassB
          + J.angularB.dot3(IwB_inv * J.angularB);
}

// ============================================================
//  applyConstraintImpulse
// ============================================================

inline void applyConstraintImpulse(double lambda,
                                   const ConstraintJacobian& J,
                                   double invMassA, const Matrix3& IwA_inv,
                                   double invMassB, const Matrix3& IwB_inv,
                                   Vector3& velA,  Vector3& omegaA,
                                   Vector3& velB,  Vector3& omegaB) {
    velA   += J.linearA  * (lambda * invMassA);
    omegaA += IwA_inv * J.angularA * lambda;
    velB   += J.linearB  * (lambda * invMassB);
    omegaB += IwB_inv * J.angularB * lambda;
}

// ============================================================
//  assembleDelassus<M>
// ============================================================

template<int M>
std::array<double, M*M> assembleDelassus(
    const std::array<ConstraintJacobian, M>& J,
    double invMassA, const Matrix3& IwA_inv,
    double invMassB, const Matrix3& IwB_inv) {
    static_assert(M > 0, "assembleDelassus: M must be positive");

    std::array<double, M*M> K{};

    for (int i = 0; i < M; ++i) {
        for (int j = i; j < M; ++j) {
            double kij =
                  J[static_cast<size_t>(i)].linearA.dot3(J[static_cast<size_t>(j)].linearA)  * invMassA
                + J[static_cast<size_t>(i)].angularA.dot3(IwA_inv * J[static_cast<size_t>(j)].angularA)
                + J[static_cast<size_t>(i)].linearB.dot3(J[static_cast<size_t>(j)].linearB)  * invMassB
                + J[static_cast<size_t>(i)].angularB.dot3(IwB_inv * J[static_cast<size_t>(j)].angularB);

            K[static_cast<size_t>(i*M + j)] = kij;
            K[static_cast<size_t>(j*M + i)] = kij;
        }
    }
    return K;
}

// ============================================================
//  effectiveMassContribution  (single-body helper)
// ============================================================

inline double effectiveMassContribution(const ConstraintJacobian& Ji,
                                        const ConstraintJacobian& Jj,
                                        double invMass, const Matrix3& IwInv,
                                        bool isBodyA) {
    const Vector3& linI = isBodyA ? Ji.linearA  : Ji.linearB;
    const Vector3& angI = isBodyA ? Ji.angularA : Ji.angularB;
    const Vector3& linJ = isBodyA ? Jj.linearA  : Jj.linearB;
    const Vector3& angJ = isBodyA ? Jj.angularA : Jj.angularB;
    return linI.dot3(linJ) * invMass + angI.dot3(IwInv * angJ);
}

}

#endif