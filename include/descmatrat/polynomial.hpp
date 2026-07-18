/// @file polynomial.hpp
/// @brief Polinomios en una variable con coeficientes exactos.
///
/// Los coeficientes se guardan en potencias descendentes, igual que en el
/// codigo MATLAB de referencia: el primer elemento es el de mayor grado y el
/// ultimo es el termino independiente. Hay dos alias:
///   - QPoly: coeficientes racionales (mpq_class).
///   - ZPoly: coeficientes enteros (mpz_class), usados en forma primitiva.
///
/// Todas las funciones son inline para que la biblioteca sea de solo cabecera y
/// se pueda usar igual en Windows y en Linux.
#pragma once

#include <gmpxx.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dmr {

/// Polinomio con coeficientes racionales, en potencias descendentes.
using QPoly = std::vector<mpq_class>;
/// Polinomio con coeficientes enteros, en potencias descendentes.
using ZPoly = std::vector<mpz_class>;

/// @brief Quita los ceros a la izquierda (coeficientes de mayor grado nulos).
/// @param p polinomio de entrada.
/// @return el polinomio sin ceros lider; el polinomio nulo queda como {0}.
inline QPoly poly_trim(const QPoly& p) {
    std::size_t i = 0;
    while (i < p.size() && p[i] == 0) ++i;
    if (i == p.size()) return QPoly{mpq_class(0)};
    return QPoly(p.begin() + i, p.end());
}

/// @brief Grado del polinomio ya recortado.
/// @return el grado, o -1 para el polinomio nulo.
inline int poly_degree(const QPoly& p) {
    QPoly t = poly_trim(p);
    if (t.size() == 1 && t[0] == 0) return -1;
    return static_cast<int>(t.size()) - 1;
}

/// @brief Indica si el polinomio es el polinomio nulo.
inline bool poly_is_zero(const QPoly& p) { return poly_degree(p) < 0; }

/// @brief Evalua el polinomio en un punto usando el esquema de Horner.
inline mpq_class poly_eval(const QPoly& p, const mpq_class& x) {
    mpq_class acc(0);
    for (const auto& c : p) acc = acc * x + c;
    return acc;
}

/// @brief Compara dos polinomios salvo ceros lider.
inline bool poly_equal(const QPoly& a, const QPoly& b) {
    QPoly ta = poly_trim(a), tb = poly_trim(b);
    if (ta.size() != tb.size()) return false;
    for (std::size_t i = 0; i < ta.size(); ++i)
        if (ta[i] != tb[i]) return false;
    return true;
}

/// @brief Divide el polinomio por su coeficiente lider para volverlo monico.
inline QPoly poly_monic(const QPoly& p) {
    QPoly t = poly_trim(p);
    if (poly_is_zero(t)) return t;
    mpq_class lead = t[0];
    for (auto& c : t) c /= lead;
    return t;
}

/// @brief Producto de dos polinomios.
inline QPoly poly_mul(const QPoly& a, const QPoly& b) {
    QPoly ta = poly_trim(a), tb = poly_trim(b);
    if (poly_is_zero(ta) || poly_is_zero(tb)) return QPoly{mpq_class(0)};
    QPoly r(ta.size() + tb.size() - 1, mpq_class(0));
    for (std::size_t i = 0; i < ta.size(); ++i)
        for (std::size_t j = 0; j < tb.size(); ++j) r[i + j] += ta[i] * tb[j];
    return r;
}

/// @brief Division con resto: encuentra q y r tales que a = q*b + r con grado(r) < grado(b).
/// @return el par (cociente, resto).
inline std::pair<QPoly, QPoly> poly_divmod(const QPoly& A, const QPoly& B) {
    QPoly a = poly_trim(A), b = poly_trim(B);
    // Trabajamos en potencias ascendentes, es mas comodo para dividir.
    // asc[i] es el coeficiente de x^i.
    const int db = poly_degree(b);
    if (db < 0) throw std::invalid_argument("poly_divmod: division por cero");
    const int da = poly_degree(a);
    if (da < 0 || da < db) return {QPoly{mpq_class(0)}, a};

    std::vector<mpq_class> asc_a(da + 1), asc_b(db + 1);
    for (int i = 0; i <= da; ++i) asc_a[i] = a[da - i];
    for (int i = 0; i <= db; ++i) asc_b[i] = b[db - i];

    std::vector<mpq_class> q(da - db + 1, mpq_class(0));
    std::vector<mpq_class> r = asc_a;
    for (int i = da; i >= db; --i) {
        mpq_class coef = r[i] / asc_b[db];
        q[i - db] = coef;
        for (int j = 0; j <= db; ++j) r[i - db + j] -= coef * asc_b[j];
    }

    // Volvemos a potencias descendentes.
    QPoly qd(q.size());
    for (std::size_t i = 0; i < q.size(); ++i) qd[i] = q[q.size() - 1 - i];
    QPoly rd(db);  // resto de grado < db
    for (int i = 0; i < db; ++i) rd[i] = r[db - 1 - i];
    return {poly_trim(qd), poly_trim(rd)};
}

/// @brief Maximo comun divisor monico de dos polinomios (algoritmo de Euclides).
inline QPoly poly_gcd(const QPoly& A, const QPoly& B) {
    QPoly a = poly_trim(A), b = poly_trim(B);
    while (!poly_is_zero(b)) {
        QPoly r = poly_divmod(a, b).second;
        a = b;
        b = r;
    }
    if (poly_is_zero(a)) return a;
    return poly_monic(a);
}

/// @brief Minimo comun multiplo monico de dos polinomios.
inline QPoly poly_lcm(const QPoly& a, const QPoly& b) {
    if (poly_is_zero(a) || poly_is_zero(b)) return QPoly{mpq_class(0)};
    QPoly g = poly_gcd(a, b);
    QPoly prod = poly_mul(a, b);
    return poly_monic(poly_divmod(prod, g).first);
}

/// @brief Potencia entera de un polinomio: p^e (con p^0 = 1).
inline QPoly poly_pow(const QPoly& p, int e) {
    QPoly result{mpq_class(1)};
    for (int i = 0; i < e; ++i) result = poly_mul(result, p);
    return result;
}

/// @brief Suma de dos polinomios.
inline QPoly poly_add(const QPoly& a, const QPoly& b) {
    const int na = static_cast<int>(a.size()), nb = static_cast<int>(b.size());
    const int n = std::max(na, nb);
    QPoly r(n, mpq_class(0));
    for (int i = 0; i < na; ++i) r[n - na + i] += a[i];
    for (int i = 0; i < nb; ++i) r[n - nb + i] += b[i];
    return poly_trim(r);
}

/// @brief Resta de dos polinomios (a - b).
inline QPoly poly_sub(const QPoly& a, const QPoly& b) {
    QPoly nb = b;
    for (auto& c : nb) c = -c;
    return poly_add(a, nb);
}

/// @brief Derivada del polinomio.
inline QPoly poly_deriv(const QPoly& p) {
    const int deg = poly_degree(p);
    if (deg < 1) return QPoly{mpq_class(0)};
    QPoly t = poly_trim(p);
    QPoly d(deg);
    for (int i = 0; i < deg; ++i) d[i] = t[i] * mpq_class(deg - i);
    return poly_trim(d);
}

/// @brief Pasa un polinomio racional a su forma entera primitiva.
///
/// El procedimiento: multiplica por el minimo comun multiplo de los
/// denominadores, divide por el maximo comun divisor de los numeradores y deja
/// el coeficiente lider positivo. El resultado representa el mismo polinomio
/// salvo un factor racional no nulo.
inline ZPoly to_integer_primitive(const QPoly& p_in) {
    QPoly p = poly_trim(p_in);

    mpz_class lcm_den(1);
    for (const auto& c : p) {
        mpz_class den = c.get_den();
        mpz_class g;
        mpz_gcd(g.get_mpz_t(), lcm_den.get_mpz_t(), den.get_mpz_t());
        lcm_den = lcm_den / g * den;
    }

    ZPoly z(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) {
        mpq_class scaled = p[i] * mpq_class(lcm_den);
        z[i] = scaled.get_num();
    }

    mpz_class content(0);
    for (const auto& c : z) {
        mpz_class a = abs(c);
        mpz_gcd(content.get_mpz_t(), content.get_mpz_t(), a.get_mpz_t());
    }
    if (content == 0) content = 1;
    for (auto& c : z) c /= content;

    if (!z.empty() && z[0] < 0)
        for (auto& c : z) c = -c;
    return z;
}

/// @brief Convierte un polinomio entero a racional.
inline QPoly to_qpoly(const ZPoly& z) {
    QPoly q(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) q[i] = mpq_class(z[i]);
    return q;
}

/// @brief Representacion legible del polinomio, por ejemplo "x^2 + 1" o "x - 2".
///
/// Sirve para coeficientes enteros o racionales.
template <class T>
std::string poly_to_string(const std::vector<T>& p, const std::string& var = "x") {
    const int n = static_cast<int>(p.size()) - 1;
    std::string out;
    bool first = true;
    for (int i = 0; i < static_cast<int>(p.size()); ++i) {
        T c = p[i];
        const int deg = n - i;
        if (c == 0) continue;
        const bool neg = c < 0;
        T mag = neg ? -c : c;
        if (first) {
            if (neg) out += "-";
        } else {
            out += neg ? " - " : " + ";
        }
        first = false;
        const bool unit_coeff = (mag == 1);
        if (!(unit_coeff && deg > 0)) {
            std::ostringstream os;
            os << mag;
            out += os.str();
            if (deg > 0) out += "*";
        }
        if (deg == 1) {
            out += var;
        } else if (deg > 0) {
            out += var;
            out += "^";
            out += std::to_string(deg);
        }
    }
    if (first) out = "0";
    return out;
}

}  // namespace dmr
