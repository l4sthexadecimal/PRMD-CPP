#include "descmatrat/polynomial.hpp"

#include <gtest/gtest.h>

#include "test_util.hpp"

using namespace dmr;
using dmr_test::qpoly;

TEST(Polynomial, DegreeAndTrim) {
    EXPECT_EQ(poly_degree(qpoly({0, 0, 3, 1})), 1);
    EXPECT_EQ(poly_degree(qpoly({1, 0, 0})), 2);
    EXPECT_EQ(poly_degree(qpoly({0})), -1);
}

TEST(Polynomial, EvalHorner) {
    // p(x) = x^2 + 1  evaluado en 3 da 10
    EXPECT_EQ(poly_eval(qpoly({1, 0, 1}), mpq_class(3)), mpq_class(10));
    // p(x) = x - 2 evaluado en 2 da 0
    EXPECT_EQ(poly_eval(qpoly({1, -2}), mpq_class(2)), mpq_class(0));
}

TEST(Polynomial, IntegerPrimitiveClearsDenominators) {
    // (1/2) x^2 + (3/4) x + 1  ->  primitivo  2 x^2 + 3 x + 4
    QPoly p = {mpq_class(1, 2), mpq_class(3, 4), mpq_class(1)};
    ZPoly z = to_integer_primitive(p);
    ASSERT_EQ(z.size(), 3u);
    EXPECT_EQ(z[0], mpz_class(2));
    EXPECT_EQ(z[1], mpz_class(3));
    EXPECT_EQ(z[2], mpz_class(4));
}

TEST(Polynomial, IntegerPrimitiveRemovesContentAndSign) {
    // -2 x^2 - 4 x - 6  ->  primitivo con lider positivo  x^2 + 2 x + 3
    QPoly p = {mpq_class(-2), mpq_class(-4), mpq_class(-6)};
    ZPoly z = to_integer_primitive(p);
    ASSERT_EQ(z.size(), 3u);
    EXPECT_EQ(z[0], mpz_class(1));
    EXPECT_EQ(z[1], mpz_class(2));
    EXPECT_EQ(z[2], mpz_class(3));
}

TEST(Polynomial, ToString) {
    EXPECT_EQ(poly_to_string(qpoly({1, 0, 1})), "x^2 + 1");
    EXPECT_EQ(poly_to_string(qpoly({1, -2})), "x - 2");
    EXPECT_EQ(poly_to_string(qpoly({1, 0, -1, 0})), "x^3 - x");
}

TEST(Polynomial, MulAndDivmod) {
    // (x - 2)(x^2 + 1) = x^3 - 2x^2 + x - 2
    QPoly prod = poly_mul(qpoly({1, -2}), qpoly({1, 0, 1}));
    EXPECT_TRUE(poly_equal(prod, qpoly({1, -2, 1, -2})));
    // dividir de vuelta da resto cero y cociente x^2 + 1
    auto qr = poly_divmod(prod, qpoly({1, -2}));
    EXPECT_TRUE(poly_equal(qr.first, qpoly({1, 0, 1})));
    EXPECT_TRUE(poly_is_zero(qr.second));
}

TEST(Polynomial, GcdAndLcm) {
    // gcd((x-1)(x-2), (x-2)(x-3)) = x - 2
    QPoly a = poly_mul(qpoly({1, -1}), qpoly({1, -2}));
    QPoly b = poly_mul(qpoly({1, -2}), qpoly({1, -3}));
    EXPECT_TRUE(poly_equal(poly_gcd(a, b), qpoly({1, -2})));
    // lcm = (x-1)(x-2)(x-3)
    QPoly expected = poly_mul(poly_mul(qpoly({1, -1}), qpoly({1, -2})), qpoly({1, -3}));
    EXPECT_TRUE(poly_equal(poly_lcm(a, b), poly_monic(expected)));
}

TEST(Polynomial, Derivative) {
    // d/dx (x^3 - 2x^2 + x - 2) = 3x^2 - 4x + 1
    EXPECT_TRUE(poly_equal(poly_deriv(qpoly({1, -2, 1, -2})), qpoly({3, -4, 1})));
}
