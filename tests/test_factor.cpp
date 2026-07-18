#include "descmatrat/factor.hpp"

#include <gtest/gtest.h>

#include "test_util.hpp"

using namespace dmr;
using dmr_test::qpoly;

TEST(Factor, IrreducibleQuadratic) {
    // x^2 + 1 es irreducible sobre Q.
    auto f = factor_rational_poly(qpoly({1, 0, 1}));
    ASSERT_EQ(f.size(), 1u);
    EXPECT_EQ(f[0].mult, 1);
    EXPECT_EQ(f[0].poly, (ZPoly{mpz_class(1), mpz_class(0), mpz_class(1)}));
}

TEST(Factor, RepeatedLinearFactor) {
    // (x-2)^2 = x^2 - 4x + 4
    auto f = factor_rational_poly(qpoly({1, -4, 4}));
    ASSERT_EQ(f.size(), 1u);
    EXPECT_EQ(f[0].mult, 2);
    EXPECT_EQ(f[0].poly, (ZPoly{mpz_class(1), mpz_class(-2)}));
}

TEST(Factor, MixedFactors) {
    // (x-2)^2 (x^2+1) = x^4 - 4x^3 + 5x^2 - 4x + 4
    auto f = factor_rational_poly(qpoly({1, -4, 5, -4, 4}));
    ASSERT_EQ(f.size(), 2u);
    // Ordenados por grado: primero (x-2), luego (x^2+1).
    EXPECT_EQ(f[0].poly, (ZPoly{mpz_class(1), mpz_class(-2)}));
    EXPECT_EQ(f[0].mult, 2);
    EXPECT_EQ(f[1].poly, (ZPoly{mpz_class(1), mpz_class(0), mpz_class(1)}));
    EXPECT_EQ(f[1].mult, 1);
}

TEST(Factor, HandlesRationalInput) {
    // (1/2)(x^2 + 1): el factor racional constante se ignora.
    QPoly p = {mpq_class(1, 2), mpq_class(0), mpq_class(1, 2)};
    auto f = factor_rational_poly(p);
    ASSERT_EQ(f.size(), 1u);
    EXPECT_EQ(f[0].poly, (ZPoly{mpz_class(1), mpz_class(0), mpz_class(1)}));
    EXPECT_EQ(f[0].mult, 1);
}
