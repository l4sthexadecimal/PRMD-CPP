#include "descmatrat/primary.hpp"

#include <gtest/gtest.h>

#include <numeric>

#include "descmatrat/charpoly.hpp"
#include "descmatrat/linalg.hpp"
#include "test_util.hpp"

using namespace dmr;
using dmr_test::mat;
using dmr_test::zpoly_mul;

namespace {

// Normaliza un polinomio entero a su forma primitiva (mismo criterio que
// to_integer_primitive), para poder comparar dos polinomios salvo escalar.
ZPoly primitive_of(const ZPoly& z) {
    QPoly q(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) q[i] = mpq_class(z[i]);
    return to_integer_primitive(q);
}

// Comprueba las propiedades esenciales de una descomposicion primaria de A.
void check_decomposition(const QMat& A) {
    const int n = static_cast<int>(A.rows());
    PrimaryDecomposition d = rational_primary_decomposition(A);

    // 1) U es cuadrada e invertible.
    ASSERT_EQ(d.U.rows(), n);
    ASSERT_EQ(d.U.cols(), n);
    EXPECT_EQ(rank(d.U), n);

    // 2) Relacion de semejanza exacta: A*U = U*H.
    EXPECT_TRUE((A * d.U == d.U * d.H));

    // 3) Los tamanos de bloque suman n.
    const int total = std::accumulate(d.block_sizes.begin(), d.block_sizes.end(), 0);
    EXPECT_EQ(total, n);

    // 4) H es diagonal por bloques: fuera de los bloques todo es cero.
    int start = 0;
    for (int size : d.block_sizes) {
        for (int i = start; i < start + size; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j < start || j >= start + size) {
                    EXPECT_EQ(d.H(i, j), mpq_class(0)) << "entrada fuera de bloque no nula";
                }
            }
        }
        start += size;
    }

    // 5) Cada bloque tiene como polinomio caracteristico exactamente p_i^{a_i}.
    start = 0;
    for (std::size_t k = 0; k < d.factors.size(); ++k) {
        const int size = d.block_sizes[k];
        QMat block = d.H.block(start, start, size, size);
        QPoly cp_block = charpoly(block);

        ZPoly expected = {mpz_class(1)};
        for (int m = 0; m < d.factors[k].mult; ++m)
            expected = zpoly_mul(expected, d.factors[k].poly);

        EXPECT_EQ(to_integer_primitive(cp_block), primitive_of(expected))
            << "el bloque no corresponde al factor primario";
        start += size;
    }
}

}  // namespace

TEST(Primary, Case3x3RotationPlusScalar) {
    // charpoly = (x^2+1)(x-2): dos componentes primarias.
    QMat A = mat({{0, -1, 0}, {1, 0, 0}, {0, 0, 2}});
    check_decomposition(A);

    PrimaryDecomposition d = rational_primary_decomposition(A);
    EXPECT_EQ(d.factors.size(), 2u);
}

TEST(Primary, Case4x4WithJordan) {
    // charpoly = (x^2+1)(x-1)^2
    QMat A = mat({{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}});
    check_decomposition(A);
}

TEST(Primary, Case4x4Dense) {
    QMat A = mat({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}});
    check_decomposition(A);
}

TEST(Primary, DiagonalizableRepeatedEigenvalue) {
    // A = diag(2,2,5): (x-2)^2 (x-5). Componente de 2 tiene dimension 2.
    QMat A = mat({{2, 0, 0}, {0, 2, 0}, {0, 0, 5}});
    check_decomposition(A);
}
