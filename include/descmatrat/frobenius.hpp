/// @file frobenius.hpp
/// @brief Forma de Frobenius (forma racional canonica) y el puente desde la
///        descomposicion primaria.
///
/// La forma de Frobenius de A es
///   F = diag(C(f_1), ..., C(f_t)),   f_1 | f_2 | ... | f_t,
/// donde los f_i son los factores invariantes de A y C(f) es la matriz
/// companera de f. Existe una matriz de cambio P tal que P^{-1} A P = F.
///
/// Aqui hay dos caminos:
///   - frobenius_form(A): construye F y P directamente por el metodo
///     constructivo del cociente (peeling ciclico) descrito en el articulo
///     LaTeX. No necesita factorizar, asi que funciona sobre Q sin FLINT.
///   - invariant_factors_from_primary / frobenius_from_primary: usan la
///     descomposicion primaria como puente. De los divisores elementales se
///     arman los factores invariantes y con ellos la forma de Frobenius.
///
/// Ambas rutas deben producir los mismos factores invariantes; se comprueba de
/// forma exacta.
#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "descmatrat/fraction_free.hpp"
#include "descmatrat/linalg.hpp"
#include "descmatrat/polynomial.hpp"
#include "descmatrat/primary.hpp"
#include "descmatrat/rational.hpp"

namespace dmr {

/// @brief Forma de Frobenius junto con su matriz de cambio de base.
struct FrobeniusForm {
    QMat P;                                  ///< cambio de base, cumple P^{-1} A P = F
    QMat F;                                  ///< forma de Frobenius (companeras en la diagonal)
    std::vector<QPoly> invariant_factors;    ///< factores invariantes f_1 | ... | f_t (ascendente)
};

/// @brief Construye una base de un complemento de W dentro de Q^n.
///
/// Recorre los vectores canonicos y agrega cada uno que aumente el rango, hasta
/// completar n columnas junto con las columnas de W.
inline QMat complement_basis(const QMat& W, int n, Reducer reduce = &rref) {
    QMat cur = W;
    QMat comp(n, 0);
    QMat I = identity(n);
    for (int j = 0; j < n && cur.cols() < n; ++j) {
        QMat cand = hconcat(cur, QMat(I.col(j)));
        if (rank(cand, reduce) > rank(cur, reduce)) {
            cur = cand;
            comp = hconcat(comp, QMat(I.col(j)));
        }
    }
    return comp;
}

/// @brief Forma de Frobenius directa por el metodo del cociente (sin factorizar).
///
/// En cada etapa se toma el cociente V/W, se calcula su operador inducido, su
/// polinomio minimal (el mayor factor invariante restante) y un vector
/// generador maximal. La cadena ciclica del cociente se levanta a V y se
/// corrige resolviendo un sistema lineal, de modo que quede exactamente un
/// bloque companero. Al terminar se ensamblan P y F.
/// @param A matriz racional cuadrada.
/// @param backend eleccion del algebra lineal exacta (racional o fraction-free).
inline FrobeniusForm frobenius_form(const QMat& A, Backend backend = Backend::Rational) {
    const int n = static_cast<int>(A.rows());
    if (A.cols() != n) throw std::invalid_argument("frobenius_form: la matriz debe ser cuadrada");
    Reducer reduce = reducer_for(backend);

    QMat B_W(n, 0);                 // base del subespacio ya cubierto
    std::vector<QMat> blocks;       // bases de cada bloque, del mayor al menor
    std::vector<QPoly> mus;         // polinomio minimal de cada bloque

    while (B_W.cols() < n) {
        const int w = static_cast<int>(B_W.cols());
        QMat B_U = complement_basis(B_W, n, reduce);
        const int u = n - w;
        QMat full = hconcat(B_W, B_U);

        // Operador inducido en el cociente V/W respecto de B_U.
        QMat Z = solve(full, A * B_U, reduce);  // n x u
        QMat Abar = Z.bottomRows(u);            // u x u

        // Generador maximal y su orden (mayor factor invariante restante).
        MaximalVector mv = maximal_vector(Abar, reduce);
        QPoly mu = mv.minpoly;
        const int d = poly_degree(mu);

        // Cadena del cociente y su levantamiento a V.
        std::vector<QVec> y(d);
        QVec cur = mv.vector;
        for (int k = 0; k < d; ++k) {
            y[k] = B_U * cur;
            cur = Abar * cur;
        }

        // T representa A restringido a W en la base B_W.
        QMat T = (w > 0) ? solve_in_column_space(B_W, A * B_W, reduce) : QMat(0, 0);

        // Residuos rho_k = coordenadas en W de (A y_k - y_{k+1}).
        std::vector<QVec> rho(d > 0 ? d - 1 : 0);
        for (int k = 0; k + 1 < d; ++k) {
            QVec resid = A * y[k] - y[k + 1];
            rho[k] = (w > 0) ? QVec(solve_in_column_space(B_W, QMat(resid), reduce).col(0)) : QVec(0);
        }

        // sigma = coordenadas en W de (A y_{d-1} + suma a_j y_j), con mu = x^d + suma a_j x^j.
        QVec sig_vec = A * y[d - 1];
        for (int j = 0; j < d; ++j) sig_vec += mu[d - j] * y[j];
        QVec sigma =
            (w > 0) ? QVec(solve_in_column_space(B_W, QMat(sig_vec), reduce).col(0)) : QVec(0);

        // Parte particular w_k: w_0 = 0, w_{k+1} = T w_k - rho_k.
        std::vector<QVec> wk(d);
        wk[0] = QVec::Zero(w);
        for (int k = 0; k + 1 < d; ++k) wk[k + 1] = (w > 0 ? QVec(T * wk[k]) : QVec(0)) - rho[k];

        // Lado derecho del sistema de correccion y solucion de mu(T) zeta = RHS.
        QVec RHS = sigma;
        if (w > 0) RHS -= T * wk[d - 1];
        for (int j = 0; j < d; ++j) RHS -= mu[d - j] * wk[j];
        QVec zeta = (w > 0) ? solve_any(poly_eval_matrix_q(mu, T), RHS, reduce) : QVec(0);

        // Vectores corregidos v_k = y_k - B_W z_k, con z_k = T^k zeta + w_k.
        QMat Vblock(n, d);
        QVec Tk_zeta = zeta;
        for (int k = 0; k < d; ++k) {
            QVec zk = (w > 0 ? Tk_zeta : QVec(0)) + wk[k];
            QVec vk = y[k];
            if (w > 0) vk -= B_W * zk;
            Vblock.col(k) = vk;
            if (w > 0) Tk_zeta = T * Tk_zeta;
        }

        blocks.push_back(Vblock);
        mus.push_back(mu);
        B_W = hconcat(B_W, Vblock);
    }

    // Los bloques salieron del mayor al menor (f_t, f_{t-1}, ...). Se invierte
    // para dejar f_1 | f_2 | ... | f_t.
    const int t = static_cast<int>(mus.size());
    FrobeniusForm result;
    result.invariant_factors.resize(t);
    QMat P(n, 0);
    QMat F = QMat::Zero(n, n);
    int off = 0;
    for (int i = 0; i < t; ++i) {
        QPoly fi = poly_monic(mus[t - 1 - i]);
        result.invariant_factors[i] = fi;
        P = hconcat(P, blocks[t - 1 - i]);
        QMat C = companion(fi);
        F.block(off, off, C.rows(), C.cols()) = C;
        off += static_cast<int>(C.rows());
    }
    result.P = P;
    result.F = F;
    return result;
}

/// @brief Factores invariantes a partir de la descomposicion primaria (el puente).
///
/// De cada primo se toman sus divisores elementales p^{e_1} >= p^{e_2} >= ...
/// El mayor factor invariante f_t junta el mayor exponente de cada primo, el
/// siguiente el segundo mayor, y asi sucesivamente. El resultado cumple
/// f_1 | f_2 | ... | f_t.
inline std::vector<QPoly> invariant_factors_from_primary(const PrimaryDecomposition& pd) {
    std::size_t t = 0;
    for (const auto& e : pd.elementary) t = std::max(t, e.exponents.size());

    std::vector<QPoly> inv(t, QPoly{mpq_class(1)});
    // idx = 0 corresponde al mayor exponente y por tanto a f_t.
    for (std::size_t idx = 0; idx < t; ++idx) {
        QPoly f{mpq_class(1)};
        for (const auto& e : pd.elementary) {
            if (idx < e.exponents.size()) f = poly_mul(f, poly_pow(e.prime, e.exponents[idx]));
        }
        inv[t - 1 - idx] = poly_monic(f);
    }
    return inv;
}

/// @brief Compara dos listas de factores invariantes.
inline bool same_invariant_factors(const std::vector<QPoly>& a, const std::vector<QPoly>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!poly_equal(a[i], b[i])) return false;
    return true;
}

/// @brief Obtiene la forma de Frobenius usando la descomposicion primaria como puente.
///
/// Primero hace la descomposicion primaria de A, de ella deriva los factores
/// invariantes (el puente) y luego construye P y F. Se comprueba que los
/// factores invariantes del puente coinciden con los de la construccion
/// directa; si no, se lanza una excepcion.
/// @param A matriz racional cuadrada.
/// @param backend eleccion del algebra lineal exacta (racional o fraction-free).
inline FrobeniusForm frobenius_from_primary(const QMat& A, Backend backend = Backend::Rational) {
    PrimaryDecomposition pd = rational_primary_decomposition(A, backend);
    std::vector<QPoly> bridge = invariant_factors_from_primary(pd);

    FrobeniusForm fr = frobenius_form(A, backend);
    if (!same_invariant_factors(bridge, fr.invariant_factors))
        throw std::runtime_error(
            "frobenius_from_primary: el puente y la construccion directa no coinciden");
    return fr;
}

}  // namespace dmr
