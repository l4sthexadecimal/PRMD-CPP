/// @file primary.hpp
/// @brief Descomposicion racional primaria de una matriz.
///
/// Dada A racional cuadrada, el espacio se parte en sus componentes primarias
///   V = ker(p_1(A)^a_1) (+) ... (+) ker(p_k(A)^a_k),
/// donde los p_i son los factores irreducibles del polinomio caracteristico y
/// a_i sus multiplicidades. Juntando bases de cada nucleo se arma la matriz de
/// cambio U, y H = U^{-1} A U queda diagonal por bloques, un bloque por
/// componente primaria.
///
/// Es la version en C++ del algoritmo 'rational_primary_decomposition' del
/// codigo MATLAB de referencia. Ademas se calculan los divisores elementales de
/// cada componente, que son el puente hacia la forma de Frobenius.
#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "descmatrat/charpoly.hpp"
#include "descmatrat/factor.hpp"
#include "descmatrat/fraction_free.hpp"
#include "descmatrat/linalg.hpp"
#include "descmatrat/polynomial.hpp"
#include "descmatrat/rational.hpp"

namespace dmr {

/// @brief Diagnostico por factor, util para inspeccionar el resultado.
struct FactorDiag {
    std::string factor;   ///< el factor irreducible como texto
    int multiplicity;     ///< su multiplicidad en el polinomio caracteristico
    int nullity;          ///< dimension del nucleo de p(A)^mult (tamano del bloque)
};

/// @brief Divisores elementales asociados a un primo p.
///
/// Cada exponente e de la lista corresponde a un bloque ciclico de orden p^e
/// dentro de la componente primaria del primo.
struct ElementaryDivisor {
    QPoly prime;                  ///< el primo p (monico)
    std::vector<int> exponents;   ///< exponentes en orden decreciente
};

/// @brief Resultado completo de la descomposicion primaria.
struct PrimaryDecomposition {
    QMat U;                                       ///< matriz de cambio de base
    QMat H;                                       ///< U^{-1} A U, diagonal por bloques
    std::vector<Factor> factors;                  ///< factores irreducibles con multiplicidad
    std::vector<FactorDiag> diagnostics;          ///< diagnostico legible por factor
    QPoly char_poly;                              ///< polinomio caracteristico de A
    std::vector<int> block_sizes;                 ///< tamano de cada bloque primario en H
    std::vector<ElementaryDivisor> elementary;    ///< divisores elementales por primo
};

/// @brief Exponentes de los divisores elementales del primo p en la matriz A.
///
/// Se obtienen de las nulidades de p(A)^j: el numero de bloques con exponente
/// mayor o igual a j es (nulidad_j menos nulidad_{j-1}) dividido por el grado de p.
inline std::vector<int> elementary_exponents(const QMat& A, const ZPoly& prime,
                                             Reducer reduce = &rref) {
    const int n = static_cast<int>(A.rows());
    const int d = static_cast<int>(prime.size()) - 1;  // grado del primo
    QMat N = poly_eval_matrix(prime, A);

    std::vector<int> blocks_ge;  // blocks_ge[j-1] = bloques con exponente >= j
    QMat Npow = identity(n);
    int prev_null = 0;
    while (true) {
        Npow = Npow * N;
        int cur_null = n - rank(Npow, reduce);
        int bge = (cur_null - prev_null) / d;
        if (bge == 0) break;
        blocks_ge.push_back(bge);
        prev_null = cur_null;
    }

    std::vector<int> exps;
    const int L = static_cast<int>(blocks_ge.size());
    for (int v = 1; v <= L; ++v) {
        int ge_v = blocks_ge[v - 1];
        int ge_next = (v < L) ? blocks_ge[v] : 0;
        int count = ge_v - ge_next;  // bloques con exponente exactamente v
        for (int t = 0; t < count; ++t) exps.push_back(v);
    }
    std::sort(exps.rbegin(), exps.rend());
    return exps;
}

/// @brief Calcula la descomposicion racional primaria de A.
/// @param A matriz racional cuadrada.
/// @param backend eleccion del algebra lineal exacta (racional o fraction-free).
/// @throw std::runtime_error si no se logra una base de n columnas independientes.
inline PrimaryDecomposition rational_primary_decomposition(const QMat& A,
                                                           Backend backend = Backend::Rational) {
    const int n = static_cast<int>(A.rows());
    if (A.cols() != n) throw std::invalid_argument("la matriz debe ser cuadrada");
    Reducer reduce = reducer_for(backend);

    PrimaryDecomposition out;
    out.char_poly = charpoly(A);
    out.factors = factor_rational_poly(out.char_poly);

    QMat U(n, 0);
    for (const Factor& f : out.factors) {
        QMat PM = poly_eval_matrix(f.poly, A);
        QMat K = mat_pow(PM, static_cast<unsigned>(f.mult));
        QMat basis = nullspace_basis(K, reduce);

        std::vector<int> keep = independent_columns(basis, reduce);
        QMat sel = select_columns(basis, keep);

        out.diagnostics.push_back(
            FactorDiag{poly_to_string(f.poly), f.mult, static_cast<int>(sel.cols())});
        out.block_sizes.push_back(static_cast<int>(sel.cols()));
        out.elementary.push_back(ElementaryDivisor{poly_monic(to_qpoly(f.poly)),
                                                   elementary_exponents(A, f.poly, reduce)});
        U = hconcat(U, sel);
    }

    // Respaldo: completar con vectores canonicos si faltan columnas.
    if (U.cols() < n) {
        QMat I = identity(n);
        for (int j = 0; j < n && U.cols() < n; ++j) {
            QMat cand = hconcat(U, QMat(I.col(j)));
            if (rank(cand, reduce) > rank(U, reduce)) U = cand;
        }
        std::vector<int> keep = independent_columns(U, reduce);
        U = select_columns(U, keep);
    }

    if (U.cols() != n)
        throw std::runtime_error(
            "no se pudieron construir n columnas independientes en la base primaria");

    out.U = U;
    out.H = solve(U, A * U, reduce);
    return out;
}

}  // namespace dmr
