// Ayudas comunes para las pruebas.
#pragma once

#include <vector>

#include "descmatrat/polynomial.hpp"
#include "descmatrat/rational.hpp"

namespace dmr_test {

// Construye una QMat a partir de filas de enteros.
inline dmr::QMat mat(const std::vector<std::vector<long>>& rows) {
    const int n = static_cast<int>(rows.size());
    const int m = static_cast<int>(rows[0].size());
    dmr::QMat A(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) A(i, j) = mpq_class(rows[i][j]);
    return A;
}

// Construye un QPoly (descendente) a partir de enteros.
inline dmr::QPoly qpoly(const std::vector<long>& c) {
    dmr::QPoly p(c.size());
    for (std::size_t i = 0; i < c.size(); ++i) p[i] = mpq_class(c[i]);
    return p;
}

// Producto de dos polinomios enteros (potencias descendentes).
inline dmr::ZPoly zpoly_mul(const dmr::ZPoly& a, const dmr::ZPoly& b) {
    dmr::ZPoly r(a.size() + b.size() - 1, mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    return r;
}

}  // namespace dmr_test
