#include "descmatrat/frobenius.hpp"

#include <gtest/gtest.h>

#include "descmatrat/charpoly.hpp"
#include "descmatrat/linalg.hpp"
#include "test_util.hpp"

using namespace dmr;
using dmr_test::mat;

namespace {

// Verifica las propiedades esenciales de la forma de Frobenius de A.
void check_frobenius(const QMat& A) {
    const int n = static_cast<int>(A.rows());
    FrobeniusForm fr = frobenius_from_primary(A);  // usa el puente y comprueba coherencia

    // 1) P es cuadrada e invertible.
    ASSERT_EQ(fr.P.rows(), n);
    ASSERT_EQ(fr.P.cols(), n);
    EXPECT_EQ(rank(fr.P), n);

    // 2) Semejanza exacta: A*P = P*F.
    EXPECT_TRUE((A * fr.P == fr.P * fr.F));

    // 3) Divisibilidad f_1 | f_2 | ... | f_t.
    for (std::size_t i = 0; i + 1 < fr.invariant_factors.size(); ++i) {
        QPoly rem = poly_divmod(fr.invariant_factors[i + 1], fr.invariant_factors[i]).second;
        EXPECT_TRUE(poly_is_zero(rem)) << "f_" << i << " no divide a f_" << (i + 1);
    }

    // 4) El producto de los factores invariantes es el polinomio caracteristico.
    QPoly prod{mpq_class(1)};
    for (const QPoly& f : fr.invariant_factors) prod = poly_mul(prod, f);
    EXPECT_TRUE(poly_equal(prod, poly_monic(charpoly(A))));

    // 5) El mayor factor invariante es el polinomio minimal: anula a A.
    const QPoly& ft = fr.invariant_factors.back();
    QMat ftA = poly_eval_matrix_q(ft, A);
    EXPECT_TRUE((ftA.array() == mpq_class(0)).all());
}

}  // namespace

TEST(Frobenius, Case3x3) {
    check_frobenius(mat({{0, -1, 0}, {1, 0, 0}, {0, 0, 2}}));
}

TEST(Frobenius, Case4x4WithJordan) {
    check_frobenius(mat({{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}}));
}

TEST(Frobenius, Case4x4Dense) {
    check_frobenius(mat({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}}));
}

TEST(Frobenius, DiagonalRepeated) {
    // diag(2,2,5): factores invariantes (x-2) y (x-2)(x-5).
    QMat A = mat({{2, 0, 0}, {0, 2, 0}, {0, 0, 5}});
    check_frobenius(A);
    FrobeniusForm fr = frobenius_form(A);
    ASSERT_EQ(fr.invariant_factors.size(), 2u);
    // f_1 = x - 2
    EXPECT_TRUE(poly_equal(fr.invariant_factors[0], dmr_test::qpoly({1, -2})));
    // f_2 = x^2 - 7x + 10
    EXPECT_TRUE(poly_equal(fr.invariant_factors[1], dmr_test::qpoly({1, -7, 10})));
}

TEST(Frobenius, IdentityForcesMultipleBlocks) {
    // La identidad 3x3 tiene tres factores invariantes iguales a (x-1).
    QMat A = identity(3);
    check_frobenius(A);
    FrobeniusForm fr = frobenius_form(A);
    EXPECT_EQ(fr.invariant_factors.size(), 3u);
    // La forma de Frobenius de la identidad es la identidad.
    EXPECT_TRUE((fr.F == identity(3)));
}

TEST(Frobenius, RecoversFromConjugation) {
    // Construimos A = S * C * S^{-1} con C companera de x^3 - 2x^2 + x - 2
    // y verificamos que se recupera esa companera como unico factor invariante.
    QPoly f = dmr_test::qpoly({1, -2, 1, -2});
    QMat C = companion(f);
    QMat S = mat({{1, 1, 0}, {0, 1, 1}, {1, 0, 1}});
    QMat A = S * C * solve(S, identity(3));
    check_frobenius(A);
    FrobeniusForm fr = frobenius_form(A);
    ASSERT_EQ(fr.invariant_factors.size(), 1u);
    EXPECT_TRUE(poly_equal(fr.invariant_factors[0], f));
}
