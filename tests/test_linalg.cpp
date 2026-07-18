#include "descmatrat/linalg.hpp"

#include <gtest/gtest.h>

#include "test_util.hpp"

using namespace dmr;
using dmr_test::mat;

TEST(Linalg, RankAndPivots) {
    // Matriz de rango 2 (tercera fila = fila1 + fila2).
    QMat A = mat({{1, 2, 3}, {0, 1, 4}, {1, 3, 7}});
    EXPECT_EQ(rank(A), 2);
}

TEST(Linalg, NullspaceBasis) {
    // A tiene nucleo de dimension 1.
    QMat A = mat({{1, 2, 3}, {2, 4, 6}, {1, 1, 1}});
    QMat N = nullspace_basis(A);
    EXPECT_EQ(N.cols(), 1);
    // Cada columna del nucleo cumple A*v = 0.
    QMat prod = A * N;
    EXPECT_TRUE((prod.array() == mpq_class(0)).all());
}

TEST(Linalg, SolveInverseTimesMatrix) {
    QMat A = mat({{2, 0}, {0, 4}});
    QMat B = mat({{2, 0}, {0, 2}});
    // X = A^{-1} B = diag(1, 1/2)
    QMat X = solve(A, B);
    EXPECT_EQ(X(0, 0), mpq_class(1));
    EXPECT_EQ(X(1, 1), mpq_class(1, 2));
    EXPECT_EQ(X(0, 1), mpq_class(0));
}

TEST(Linalg, SolveThrowsOnSingular) {
    QMat A = mat({{1, 2}, {2, 4}});
    QMat B = mat({{1, 0}, {0, 1}});
    EXPECT_THROW(solve(A, B), std::runtime_error);
}

TEST(Linalg, MatPow) {
    QMat A = mat({{2, 0}, {0, 3}});
    QMat P = mat_pow(A, 3);
    EXPECT_EQ(P(0, 0), mpq_class(8));
    EXPECT_EQ(P(1, 1), mpq_class(27));
    // A^0 = identidad
    QMat I = mat_pow(A, 0);
    EXPECT_EQ(I(0, 0), mpq_class(1));
    EXPECT_EQ(I(1, 1), mpq_class(1));
}

TEST(Linalg, PolyEvalMatrix) {
    // p(x) = x - 2 evaluado en A = diag(2, 5) da diag(0, 3).
    QMat A = mat({{2, 0}, {0, 5}});
    ZPoly p = {mpz_class(1), mpz_class(-2)};
    QMat PA = poly_eval_matrix(p, A);
    EXPECT_EQ(PA(0, 0), mpq_class(0));
    EXPECT_EQ(PA(1, 1), mpq_class(3));
}
