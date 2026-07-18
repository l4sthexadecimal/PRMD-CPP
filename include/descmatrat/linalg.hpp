/// @file linalg.hpp
/// @brief Algebra lineal exacta sobre los racionales.
///
/// Todo se hace con eliminacion de Gauss-Jordan sobre mpq_class, asi que no hay
/// error de redondeo en ningun paso. Es la misma estrategia del codigo MATLAB
/// de referencia (rref, base del nucleo, columnas independientes, solve
/// exacto), mas las piezas que hacen falta para construir la forma de Frobenius
/// (matriz companera, polinomio minimal de un vector y de un operador, y vector
/// generador maximal).
///
/// Las operaciones que se apoyan en la forma escalonada reducida aceptan un
/// "reductor" (Reducer): una funcion que calcula esa forma. Por defecto se usa
/// la version racional rref, pero se puede pasar la version fraction-free
/// (ff_rref, en fraction_free.hpp) para reducir el crecimiento de coeficientes.
/// Como la forma escalonada reducida es unica, ambos reductores dan el mismo
/// resultado y las funciones de arriba son intercambiables.
///
/// Funciones inline: la biblioteca es de solo cabecera.
#pragma once

#include <stdexcept>
#include <vector>

#include "descmatrat/polynomial.hpp"
#include "descmatrat/rational.hpp"

namespace dmr {

/// @brief Resultado de la forma escalonada reducida por filas.
struct RrefResult {
    QMat R;                    ///< la matriz ya reducida
    std::vector<int> pivots;   ///< indices de las columnas pivote, en orden creciente
};

/// @brief Forma escalonada reducida por filas (Gauss-Jordan racional exacto).
inline RrefResult rref(QMat M) {
    const int rows = static_cast<int>(M.rows());
    const int cols = static_cast<int>(M.cols());
    std::vector<int> pivots;
    int row = 0;
    for (int col = 0; col < cols && row < rows; ++col) {
        int sel = -1;
        for (int i = row; i < rows; ++i) {
            if (M(i, col) != 0) {
                sel = i;
                break;
            }
        }
        if (sel < 0) continue;

        M.row(row).swap(M.row(sel));

        mpq_class pivot = M(row, col);
        for (int j = 0; j < cols; ++j) M(row, j) /= pivot;

        for (int i = 0; i < rows; ++i) {
            if (i == row) continue;
            mpq_class factor = M(i, col);
            if (factor != 0)
                for (int j = 0; j < cols; ++j) M(i, j) -= factor * M(row, j);
        }
        pivots.push_back(col);
        ++row;
    }
    return {M, pivots};
}

/// @brief Tipo de un reductor a forma escalonada reducida.
using Reducer = RrefResult (*)(QMat);

/// @brief Rango exacto (numero de columnas pivote).
inline int rank(const QMat& A, Reducer reduce = &rref) {
    return static_cast<int>(reduce(A).pivots.size());
}

/// @brief Indices de un conjunto maximal de columnas independientes.
inline std::vector<int> independent_columns(const QMat& A, Reducer reduce = &rref) {
    return reduce(A).pivots;
}

/// @brief Copia las columnas indicadas en una nueva matriz.
inline QMat select_columns(const QMat& A, const std::vector<int>& cols) {
    QMat out(A.rows(), static_cast<int>(cols.size()));
    for (std::size_t k = 0; k < cols.size(); ++k) out.col(static_cast<int>(k)) = A.col(cols[k]);
    return out;
}

/// @brief Concatena dos matrices lado a lado (misma cantidad de filas).
inline QMat hconcat(const QMat& A, const QMat& B) {
    if (A.cols() == 0) return B;
    if (B.cols() == 0) return A;
    if (A.rows() != B.rows()) throw std::invalid_argument("hconcat: filas distintas");
    QMat out(A.rows(), A.cols() + B.cols());
    out.leftCols(A.cols()) = A;
    out.rightCols(B.cols()) = B;
    return out;
}

/// @brief Base del nucleo (espacio nulo) de A.
/// @return matriz cuyas columnas v cumplen A*v = 0; sin columnas si el nucleo es trivial.
inline QMat nullspace_basis(const QMat& A, Reducer reduce = &rref) {
    RrefResult rr = reduce(A);
    const int cols = static_cast<int>(A.cols());

    std::vector<bool> is_pivot(cols, false);
    for (int p : rr.pivots) is_pivot[p] = true;

    std::vector<int> free_cols;
    for (int c = 0; c < cols; ++c)
        if (!is_pivot[c]) free_cols.push_back(c);

    QMat B = QMat::Zero(cols, static_cast<int>(free_cols.size()));
    for (std::size_t k = 0; k < free_cols.size(); ++k) {
        const int f = free_cols[k];
        B(f, static_cast<int>(k)) = mpq_class(1);
        for (std::size_t i = 0; i < rr.pivots.size(); ++i)
            B(rr.pivots[i], static_cast<int>(k)) = -rr.R(static_cast<int>(i), f);
    }
    return B;
}

/// @brief Resuelve A*X = B de forma exacta con A cuadrada e invertible.
/// @throw std::runtime_error si A no es invertible.
inline QMat solve(const QMat& A, const QMat& B, Reducer reduce = &rref) {
    const int n = static_cast<int>(A.rows());
    if (A.cols() != n) throw std::invalid_argument("solve: A debe ser cuadrada");
    if (B.rows() != n) throw std::invalid_argument("solve: dimensiones incompatibles");
    const int m = static_cast<int>(B.cols());

    QMat aug(n, n + m);
    aug.leftCols(n) = A;
    aug.rightCols(m) = B;

    RrefResult rr = reduce(aug);
    if (static_cast<int>(rr.pivots.size()) != n)
        throw std::runtime_error("solve: la matriz no es invertible");
    for (int i = 0; i < n; ++i)
        if (rr.pivots[i] != i) throw std::runtime_error("solve: la matriz no es invertible");
    return rr.R.rightCols(m);
}

/// @brief Resuelve B*X = RHS cuando B tiene columnas independientes y RHS esta en su espacio columna.
/// @throw std::runtime_error si B no tiene rango columna completo o RHS no pertenece a su imagen.
inline QMat solve_in_column_space(const QMat& B, const QMat& RHS, Reducer reduce = &rref) {
    const int n = static_cast<int>(B.rows());
    const int w = static_cast<int>(B.cols());
    const int m = static_cast<int>(RHS.cols());
    if (w == 0) {
        if (RHS.size() > 0 && !(RHS.array() == mpq_class(0)).all())
            throw std::runtime_error("solve_in_column_space: RHS no esta en el espacio (base vacia)");
        return QMat(0, m);
    }
    QMat aug(n, w + m);
    aug.leftCols(w) = B;
    aug.rightCols(m) = RHS;
    RrefResult rr = reduce(aug);
    if (static_cast<int>(rr.pivots.size()) != w)
        throw std::runtime_error("solve_in_column_space: base sin rango columna completo o RHS fuera de imagen");
    for (int i = 0; i < w; ++i)
        if (rr.pivots[i] != i)
            throw std::runtime_error("solve_in_column_space: RHS fuera del espacio columna");
    return rr.R.block(0, w, w, m);
}

/// @brief Devuelve una solucion particular de M*x = rhs suponiendo el sistema consistente.
///
/// Las variables libres se fijan en cero. Util cuando M puede ser singular pero
/// el sistema es compatible (por ejemplo en el paso de correccion de Frobenius).
/// @throw std::runtime_error si el sistema es incompatible.
inline QVec solve_any(const QMat& M, const QVec& rhs, Reducer reduce = &rref) {
    const int rows = static_cast<int>(M.rows());
    const int cols = static_cast<int>(M.cols());
    QMat aug(rows, cols + 1);
    if (cols > 0) aug.leftCols(cols) = M;
    aug.col(cols) = rhs;
    RrefResult rr = reduce(aug);
    for (int p : rr.pivots)
        if (p == cols) throw std::runtime_error("solve_any: sistema incompatible");
    QVec x = QVec::Zero(cols);
    for (std::size_t i = 0; i < rr.pivots.size(); ++i) {
        int c = rr.pivots[i];
        if (c < cols) x(c) = rr.R(static_cast<int>(i), cols);
    }
    return x;
}

/// @brief Potencia entera de una matriz cuadrada: A^e (con A^0 = identidad).
inline QMat mat_pow(const QMat& A, unsigned e) {
    QMat result = identity(static_cast<int>(A.rows()));
    QMat base = A;
    while (e > 0) {
        if (e & 1u) result = result * base;
        e >>= 1u;
        if (e > 0) base = base * base;
    }
    return result;
}

/// @brief Evalua un polinomio entero (potencias descendentes) en la matriz A, es decir p(A).
inline QMat poly_eval_matrix(const ZPoly& p, const QMat& A) {
    const int n = static_cast<int>(A.rows());
    QMat I = identity(n);
    QMat PM = QMat::Zero(n, n);
    for (const auto& c : p) PM = PM * A + mpq_class(c) * I;
    return PM;
}

/// @brief Evalua un polinomio racional (potencias descendentes) en la matriz A.
inline QMat poly_eval_matrix_q(const QPoly& p, const QMat& A) {
    const int n = static_cast<int>(A.rows());
    QMat I = identity(n);
    QMat PM = QMat::Zero(n, n);
    for (const auto& c : p) PM = PM * A + c * I;
    return PM;
}

/// @brief Matriz companera de un polinomio monico f.
///
/// En la base (v, A v, ..., A^{d-1} v) de un subespacio ciclico, la relacion
/// A*v_k = v_{k+1} y A*v_{d-1} = -suma(a_j v_j) produce esta matriz, donde los
/// a_j son los coeficientes de f = x^d + a_{d-1} x^{d-1} + ... + a_0.
inline QMat companion(const QPoly& f_in) {
    QPoly f = poly_monic(f_in);
    const int d = poly_degree(f);
    if (d < 1) throw std::invalid_argument("companion: grado invalido");
    QMat C = QMat::Zero(d, d);
    for (int i = 1; i < d; ++i) C(i, i - 1) = mpq_class(1);  // subdiagonal de unos
    for (int i = 0; i < d; ++i) C(i, d - 1) = -f[d - i];     // ultima columna: -a_i
    return C;
}

/// @brief Indica si el vector v pertenece al espacio columna de B (con B de rango columna completo).
inline bool in_column_space(const QMat& B, const QVec& v, Reducer reduce = &rref) {
    QMat aug(B.rows(), B.cols() + 1);
    if (B.cols() > 0) aug.leftCols(B.cols()) = B;
    aug.col(static_cast<int>(B.cols())) = v;
    return rank(aug, reduce) == static_cast<int>(B.cols());
}

/// @brief Polinomio minimal de un vector bajo un operador.
///
/// Es el polinomio monico de menor grado m tal que m(M) v = 0. Se calcula
/// construyendo la sucesion de Krylov v, M v, M^2 v, ... hasta la primera
/// dependencia lineal.
inline QPoly vector_min_poly(const QMat& M, const QVec& v, Reducer reduce = &rref) {
    const int n = static_cast<int>(M.rows());
    if ((v.array() == mpq_class(0)).all()) return QPoly{mpq_class(1)};

    std::vector<QVec> ks;
    ks.push_back(v);
    while (true) {
        const int m = static_cast<int>(ks.size());
        if (m - 1 > 0) {
            QMat Kprev(n, m - 1);
            for (int i = 0; i < m - 1; ++i) Kprev.col(i) = ks[i];
            const QVec& last = ks.back();
            if (in_column_space(Kprev, last, reduce)) {
                QMat a = solve_in_column_space(Kprev, QMat(last), reduce);
                const int d = m - 1;
                QPoly mp(d + 1);
                mp[0] = mpq_class(1);
                for (int j = 0; j < d; ++j) mp[d - j] = -a(j, 0);
                return mp;
            }
        }
        ks.push_back(M * ks.back());
        if (static_cast<int>(ks.size()) > n + 1)
            throw std::runtime_error("vector_min_poly: guardia de tamano");
    }
}

/// @brief Resultado de la busqueda de un vector generador maximal.
struct MaximalVector {
    QVec vector;   ///< un vector cuyo orden es el polinomio minimal del operador
    QPoly minpoly; ///< el polinomio minimal del operador
};

/// @brief Construye un vector cuyo polinomio minimal es el minimal del operador.
///
/// Recorre los vectores canonicos y los va combinando: si e_i aporta factores
/// nuevos al orden acumulado, se busca un escalar entero t tal que el vector
/// corriente mas t*e_i alcance el minimo comun multiplo de los ordenes. Al
/// terminar, el orden acumulado es el polinomio minimal del operador y el
/// vector obtenido lo realiza. Todo se verifica de forma exacta antes de
/// aceptar cada paso.
inline MaximalVector maximal_vector(const QMat& M, Reducer reduce = &rref) {
    const int n = static_cast<int>(M.rows());
    if (n == 0) return {QVec(0), QPoly{mpq_class(1)}};

    QVec c = unit_vector(n, 0);
    QPoly g = vector_min_poly(M, c, reduce);
    for (int i = 1; i < n; ++i) {
        QVec ei = unit_vector(n, i);
        QPoly mi = vector_min_poly(M, ei, reduce);
        QPoly L = poly_lcm(g, mi);
        if (poly_degree(L) <= poly_degree(g)) continue;

        bool found = false;
        if (poly_equal(mi, L)) {
            c = ei;
            g = L;
            found = true;
        }
        for (int step = 1; step <= 200 && !found; ++step) {
            for (int sign : {1, -1}) {
                QVec cand = c + mpq_class(sign * step) * ei;
                if (poly_equal(vector_min_poly(M, cand, reduce), L)) {
                    c = cand;
                    g = L;
                    found = true;
                    break;
                }
            }
        }
        if (!found) throw std::runtime_error("maximal_vector: no se hallo combinacion adecuada");
    }
    return {c, g};
}

}  // namespace dmr
