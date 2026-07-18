#include "descmatrat/fraction_free.hpp"

#include <gtest/gtest.h>

#include <random>

#include "descmatrat/charpoly.hpp"
#include "descmatrat/frobenius.hpp"
#include "descmatrat/linalg.hpp"
#include "descmatrat/primary.hpp"
#include "test_util.hpp"

using namespace dmr;
using dmr_test::mat;

namespace {

// Genera una matriz entera aleatoria de tamano rows x cols con un generador
// determinista (misma semilla siempre, para que las pruebas sean reproducibles).
QMat random_int_matrix(std::mt19937& rng, int rows, int cols, int lo = -5, int hi = 5) {
    std::uniform_int_distribution<int> dist(lo, hi);
    QMat A(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) A(i, j) = mpq_class(dist(rng));
    return A;
}

bool matrices_equal(const QMat& A, const QMat& B) {
    return A.rows() == B.rows() && A.cols() == B.cols() && (A == B);
}

}  // namespace

TEST(FractionFree, RrefMatchesRationalOnRandomMatrices) {
    std::mt19937 rng(12345);
    for (int trial = 0; trial < 60; ++trial) {
        int rows = 1 + trial % 5;
        int cols = 1 + (trial / 5) % 6;
        QMat A = random_int_matrix(rng, rows, cols);

        // A veces forzamos filas dependientes para cubrir el caso deficiente en rango.
        if (trial % 3 == 0 && rows >= 2) A.row(rows - 1) = A.row(0) + A.row(rows - 2);

        RrefResult a = rref(A);
        RrefResult b = ff_rref(A);
        EXPECT_EQ(a.pivots, b.pivots) << "pivotes distintos en el intento " << trial;
        EXPECT_TRUE(matrices_equal(a.R, b.R)) << "RREF distinta en el intento " << trial;
    }
}

TEST(FractionFree, RrefWithRationalEntries) {
    // Con fracciones en la entrada el resultado tambien debe coincidir.
    QMat A(2, 3);
    A << mpq_class(1, 2), mpq_class(1, 3), 1, mpq_class(1, 5), 2, mpq_class(3, 7);
    EXPECT_TRUE(matrices_equal(rref(A).R, ff_rref(A).R));
}

TEST(FractionFree, SolveMatchesRational) {
    std::mt19937 rng(999);
    for (int trial = 0; trial < 30; ++trial) {
        int n = 2 + trial % 4;
        QMat A = random_int_matrix(rng, n, n);
        if (rank(A) != n) continue;  // solo sistemas invertibles
        QMat B = random_int_matrix(rng, n, 2);
        QMat x_rat = solve(A, B);
        QMat x_ff = ff_solve(A, B);
        EXPECT_TRUE(matrices_equal(x_rat, x_ff));
        EXPECT_TRUE(matrices_equal(A * x_ff, B));  // reconstruccion exacta
    }
}

TEST(FractionFree, NullspaceMatchesRational) {
    QMat A = mat({{1, 2, 3}, {2, 4, 6}, {1, 1, 1}});
    QMat n_rat = nullspace_basis(A);
    QMat n_ff = ff_nullspace(A);
    EXPECT_TRUE(matrices_equal(n_rat, n_ff));
    EXPECT_TRUE(((A * n_ff).array() == mpq_class(0)).all());
}

TEST(FractionFree, Determinant) {
    // det(diag(2,3,4)) = 24
    EXPECT_EQ(ff_determinant(mat({{2, 0, 0}, {0, 3, 0}, {0, 0, 4}})), mpq_class(24));
    // Matriz singular: determinante cero.
    EXPECT_EQ(ff_determinant(mat({{1, 2}, {2, 4}})), mpq_class(0));

    // Comparacion con el termino independiente del polinomio caracteristico:
    // det(A) = (-1)^n * cp[n].
    std::mt19937 rng(7);
    for (int trial = 0; trial < 20; ++trial) {
        int n = 2 + trial % 4;
        QMat A = random_int_matrix(rng, n, n);
        QPoly cp = charpoly(A);
        mpq_class det_from_cp = (n % 2 == 0 ? mpq_class(1) : mpq_class(-1)) * cp[n];
        EXPECT_EQ(ff_determinant(A), det_from_cp) << "det incorrecto en intento " << trial;
    }
}

TEST(FractionFree, PrimaryBackendGivesSameAndReconstructs) {
    std::vector<QMat> cases = {
        mat({{0, -1, 0}, {1, 0, 0}, {0, 0, 2}}),
        mat({{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}}),
        mat({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}}),
    };
    for (const QMat& A : cases) {
        PrimaryDecomposition r = rational_primary_decomposition(A, Backend::Rational);
        PrimaryDecomposition f = rational_primary_decomposition(A, Backend::FractionFree);
        EXPECT_TRUE(matrices_equal(r.U, f.U));
        EXPECT_TRUE(matrices_equal(r.H, f.H));
        EXPECT_TRUE(matrices_equal(A * f.U, f.U * f.H));  // reconstruccion fiable
    }
}

TEST(FractionFree, FrobeniusBackendGivesSameAndReconstructs) {
    std::vector<QMat> cases = {
        mat({{0, -1, 0}, {1, 0, 0}, {0, 0, 2}}),
        mat({{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}}),
        mat({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}}),
        identity(3),
    };
    for (const QMat& A : cases) {
        FrobeniusForm r = frobenius_from_primary(A, Backend::Rational);
        FrobeniusForm f = frobenius_from_primary(A, Backend::FractionFree);
        EXPECT_TRUE(matrices_equal(r.P, f.P));
        EXPECT_TRUE(matrices_equal(r.F, f.F));
        EXPECT_TRUE(matrices_equal(A * f.P, f.P * f.F));  // reconstruccion fiable
    }
}
