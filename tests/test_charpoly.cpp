#include "descmatrat/charpoly.hpp"

#include <gtest/gtest.h>

#include "descmatrat/linalg.hpp"
#include "test_util.hpp"

using namespace dmr;
using dmr_test::mat;

TEST(Charpoly, Diagonal) {
    // A = diag(2, 3) -> (x-2)(x-3) = x^2 - 5x + 6
    QMat A = mat({{2, 0}, {0, 3}});
    QPoly cp = charpoly(A);
    ASSERT_EQ(cp.size(), 3u);
    EXPECT_EQ(cp[0], mpq_class(1));
    EXPECT_EQ(cp[1], mpq_class(-5));
    EXPECT_EQ(cp[2], mpq_class(6));
}

TEST(Charpoly, Rotation) {
    // A = [[0,-1],[1,0]] -> x^2 + 1
    QMat A = mat({{0, -1}, {1, 0}});
    QPoly cp = charpoly(A);
    ASSERT_EQ(cp.size(), 3u);
    EXPECT_EQ(cp[0], mpq_class(1));
    EXPECT_EQ(cp[1], mpq_class(0));
    EXPECT_EQ(cp[2], mpq_class(1));
}

TEST(Charpoly, CayleyHamilton) {
    // El polinomio caracteristico evaluado en la propia matriz da cero.
    QMat A = mat({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}});
    QPoly cp = charpoly(A);
    ZPoly z(cp.size());
    for (std::size_t i = 0; i < cp.size(); ++i) z[i] = cp[i].get_num();  // cp es monico entero aqui
    QMat pA = poly_eval_matrix(z, A);
    EXPECT_TRUE((pA.array() == mpq_class(0)).all());
}
