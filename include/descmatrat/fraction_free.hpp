/// @file fraction_free.hpp
/// @brief Variante fraction-free (eliminacion de Bareiss) del algebra lineal.
///
/// La eliminacion racional directa puede hacer crecer mucho numeradores y
/// denominadores intermedios (coefficient swell). La eliminacion fraction-free
/// de Bareiss trabaja con enteros y divide de forma exacta por el pivote
/// anterior, de modo que los valores intermedios se mantienen acotados. Al
/// final se hace una sustitucion hacia atras para obtener la misma forma
/// escalonada reducida.
///
/// Como la forma escalonada reducida por filas es unica, ff_rref produce
/// exactamente el mismo resultado que rref; solo cambia el camino de calculo.
/// Por eso ff_rref se puede pasar como reductor a cualquier funcion de
/// linalg.hpp, y de ahi que la descomposicion primaria y la forma de Frobenius
/// admitan un backend fraction-free sin cambiar su logica.
///
/// Todo es de solo cabecera y portable (solo GMP y Eigen).
#pragma once

#include <stdexcept>
#include <vector>

#include "descmatrat/linalg.hpp"
#include "descmatrat/rational.hpp"

namespace dmr {

/// @brief Division entera exacta con verificacion.
/// @throw std::runtime_error si la division no es exacta (senal de un error de logica).
inline mpz_class exact_div(const mpz_class& num, const mpz_class& den) {
    mpz_class q, r;
    mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());
    if (r != 0) throw std::runtime_error("fraction_free: division no exacta");
    return q;
}

/// @brief Convierte una matriz racional en una rejilla de enteros.
///
/// Cada fila se multiplica por el minimo comun multiplo de los denominadores de
/// esa fila. Escalar una fila es una operacion elemental, asi que la forma
/// escalonada reducida no cambia.
inline std::vector<std::vector<mpz_class>> to_integer_grid(const QMat& M) {
    const int m = static_cast<int>(M.rows());
    const int n = static_cast<int>(M.cols());
    std::vector<std::vector<mpz_class>> Z(m, std::vector<mpz_class>(n));
    for (int i = 0; i < m; ++i) {
        mpz_class L(1);
        for (int j = 0; j < n; ++j) {
            mpz_class d = M(i, j).get_den();
            mpz_class g;
            mpz_gcd(g.get_mpz_t(), L.get_mpz_t(), d.get_mpz_t());
            L = L / g * d;
        }
        for (int j = 0; j < n; ++j) {
            mpq_class v = M(i, j) * mpq_class(L);
            Z[i][j] = v.get_num();
        }
    }
    return Z;
}

/// @brief Forma escalonada reducida por filas calculada de forma fraction-free.
///
/// Primero se limpian denominadores, luego se hace la eliminacion hacia
/// adelante de Bareiss (entera, con division exacta por el pivote anterior) y
/// por ultimo la sustitucion hacia atras para dejar la forma reducida. El
/// resultado es identico al de rref.
inline RrefResult ff_rref(QMat M) {
    const int m = static_cast<int>(M.rows());
    const int n = static_cast<int>(M.cols());
    std::vector<std::vector<mpz_class>> Z = to_integer_grid(M);

    // Eliminacion hacia adelante de Bareiss.
    mpz_class prev(1);
    int prow = 0;
    std::vector<int> pivots;
    for (int col = 0; col < n && prow < m; ++col) {
        int sel = -1;
        for (int r = prow; r < m; ++r)
            if (Z[r][col] != 0) {
                sel = r;
                break;
            }
        if (sel < 0) continue;
        if (sel != prow) std::swap(Z[prow], Z[sel]);

        for (int i = prow + 1; i < m; ++i) {
            for (int j = col + 1; j < n; ++j)
                Z[i][j] = exact_div(Z[prow][col] * Z[i][j] - Z[i][col] * Z[prow][j], prev);
            Z[i][col] = 0;
        }
        prev = Z[prow][col];
        pivots.push_back(col);
        ++prow;
    }

    // Sustitucion hacia atras sobre Q para dejar la forma reducida.
    QMat R = QMat::Zero(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j) R(i, j) = mpq_class(Z[i][j]);

    for (int idx = static_cast<int>(pivots.size()) - 1; idx >= 0; --idx) {
        const int pr = idx;
        const int pc = pivots[idx];
        mpq_class pivot = R(pr, pc);
        for (int j = 0; j < n; ++j) R(pr, j) /= pivot;
        for (int i = 0; i < m; ++i) {
            if (i == pr) continue;
            mpq_class factor = R(i, pc);
            if (factor != 0)
                for (int j = 0; j < n; ++j) R(i, j) -= factor * R(pr, j);
        }
    }
    return {R, pivots};
}

/// @brief Backends disponibles para el algebra lineal exacta.
enum class Backend {
    Rational,      ///< eliminacion racional directa (rref)
    FractionFree   ///< eliminacion de Bareiss (ff_rref)
};

/// @brief Devuelve el reductor correspondiente a un backend.
inline Reducer reducer_for(Backend b) {
    return (b == Backend::FractionFree) ? &ff_rref : &rref;
}

// --- Envoltorios comodos para usar el backend fraction-free directamente ---

/// @brief Resuelve A*X = B (A invertible) de forma fraction-free.
inline QMat ff_solve(const QMat& A, const QMat& B) { return solve(A, B, &ff_rref); }

/// @brief Base del nucleo de A calculada de forma fraction-free.
inline QMat ff_nullspace(const QMat& A) { return nullspace_basis(A, &ff_rref); }

/// @brief Rango de A calculado de forma fraction-free.
inline int ff_rank(const QMat& A) { return static_cast<int>(ff_rref(A).pivots.size()); }

/// @brief Determinante exacto de una matriz cuadrada por eliminacion de Bareiss.
///
/// Devuelve el determinante como racional (es entero si la matriz es entera).
inline mpq_class ff_determinant(const QMat& A) {
    const int n = static_cast<int>(A.rows());
    if (A.cols() != n) throw std::invalid_argument("ff_determinant: matriz no cuadrada");
    if (n == 0) return mpq_class(1);

    std::vector<std::vector<mpz_class>> Z = to_integer_grid(A);
    // det(Z) = (producto de los factores de fila L_i) * det(A).
    mpq_class row_scale(1);
    for (int i = 0; i < n; ++i) {
        mpz_class L(1);
        for (int j = 0; j < n; ++j) {
            mpz_class d = A(i, j).get_den();
            mpz_class g;
            mpz_gcd(g.get_mpz_t(), L.get_mpz_t(), d.get_mpz_t());
            L = L / g * d;
        }
        row_scale *= mpq_class(L);
    }

    mpz_class prev(1);
    int sign = 1;
    for (int col = 0; col < n; ++col) {
        int sel = -1;
        for (int r = col; r < n; ++r)
            if (Z[r][col] != 0) {
                sel = r;
                break;
            }
        if (sel < 0) return mpq_class(0);  // columna nula, determinante cero
        if (sel != col) {
            std::swap(Z[col], Z[sel]);
            sign = -sign;
        }
        for (int i = col + 1; i < n; ++i) {
            for (int j = col + 1; j < n; ++j)
                Z[i][j] = exact_div(Z[col][col] * Z[i][j] - Z[i][col] * Z[col][j], prev);
            Z[i][col] = 0;
        }
        prev = Z[col][col];
    }
    mpq_class det_scaled = mpq_class(sign) * mpq_class(prev);  // det de la rejilla entera
    return det_scaled / row_scale;
}

}  // namespace dmr
