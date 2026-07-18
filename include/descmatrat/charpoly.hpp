/// @file charpoly.hpp
/// @brief Polinomio caracteristico exacto de una matriz racional.
#pragma once

#include <stdexcept>

#include "descmatrat/polynomial.hpp"
#include "descmatrat/rational.hpp"

namespace dmr {

/// @brief Polinomio caracteristico monico det(x*I - A), en potencias descendentes.
///
/// Se calcula con el metodo de Faddeev-LeVerrier, que solo usa productos y
/// trazas de matrices, de modo que el resultado es exacto. Devuelve el vector
/// [1, c1, ..., cn].
inline QPoly charpoly(const QMat& A) {
    const int n = static_cast<int>(A.rows());
    if (A.cols() != n) throw std::invalid_argument("charpoly: la matriz debe ser cuadrada");

    QMat I = identity(n);
    QMat B = I;
    QPoly cp(n + 1);
    cp[0] = mpq_class(1);

    for (int k = 1; k <= n; ++k) {
        QMat AB = A * B;
        mpq_class tr(0);
        for (int i = 0; i < n; ++i) tr += AB(i, i);
        mpq_class ck = -tr / mpq_class(k);
        cp[k] = ck;
        if (k < n) B = AB + ck * I;
    }
    return cp;
}

}  // namespace dmr
