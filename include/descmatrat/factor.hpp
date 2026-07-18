/// @file factor.hpp
/// @brief Factorizacion de polinomios sobre los racionales.
///
/// Hay dos rutas:
///   - Con FLINT (si se define la macro DESCMATRAT_USE_FLINT): factorizacion
///     exacta y completa para cualquier grado, usando fmpz_poly_factor. Es el
///     mismo motor que usa SymEngine por debajo.
///   - Sin FLINT (respaldo portable): factorizacion libre de cuadrados mas
///     extraccion de raices racionales. Es correcta pero puede dejar sin
///     separar dos factores irreducibles del mismo grado sin raiz racional
///     (por ejemplo (x^2+1)(x^2+2)). Sirve para que la biblioteca compile y
///     funcione en Windows sin FLINT en los casos habituales.
///
/// Para la forma de Frobenius no hace falta factorizar; la factorizacion solo
/// interviene en la descomposicion primaria.
#pragma once

#ifdef DESCMATRAT_USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#endif

#include <algorithm>
#include <vector>

#include "descmatrat/polynomial.hpp"

namespace dmr {

/// @brief Un factor irreducible junto con su multiplicidad.
struct Factor {
    ZPoly poly;  ///< factor en forma entera primitiva (potencias descendentes)
    int mult;    ///< cuantas veces aparece
};

/// @brief Orden determinista de factores: por grado y luego por coeficientes.
inline bool factor_less(const Factor& a, const Factor& b) {
    if (a.poly.size() != b.poly.size()) return a.poly.size() < b.poly.size();
    for (std::size_t i = 0; i < a.poly.size(); ++i)
        if (a.poly[i] != b.poly[i]) return a.poly[i] < b.poly[i];
    return false;
}

#ifdef DESCMATRAT_USE_FLINT

namespace detail {

/// Convierte un ZPoly descendente a fmpz_poly de FLINT (que guarda ascendente).
inline void to_fmpz_poly(const ZPoly& z, fmpz_poly_t out) {
    const int deg = static_cast<int>(z.size()) - 1;
    fmpz_t tmp;
    fmpz_init(tmp);
    for (int k = 0; k <= deg; ++k) {
        fmpz_set_mpz(tmp, z[deg - k].get_mpz_t());
        fmpz_poly_set_coeff_fmpz(out, k, tmp);
    }
    fmpz_clear(tmp);
}

/// Convierte un fmpz_poly de FLINT a un ZPoly descendente.
inline ZPoly from_fmpz_poly(const fmpz_poly_struct* fp) {
    const long deg = fmpz_poly_degree(fp);
    ZPoly z(static_cast<std::size_t>(deg) + 1);
    fmpz_t c;
    fmpz_init(c);
    mpz_t m;
    mpz_init(m);
    for (long k = 0; k <= deg; ++k) {
        fmpz_poly_get_coeff_fmpz(c, fp, k);
        fmpz_get_mpz(m, c);
        z[static_cast<std::size_t>(deg - k)] = mpz_class(m);
    }
    mpz_clear(m);
    fmpz_clear(c);
    return z;
}

}  // namespace detail

/// @brief Factoriza un polinomio racional en irreducibles sobre Q usando FLINT.
inline std::vector<Factor> factor_rational_poly(const QPoly& p) {
    ZPoly z = to_integer_primitive(p);

    fmpz_poly_t poly;
    fmpz_poly_init(poly);
    detail::to_fmpz_poly(z, poly);

    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, poly);

    std::vector<Factor> result;
    for (long i = 0; i < fac->num; ++i) {
        const fmpz_poly_struct* fp = fac->p + i;
        if (fmpz_poly_degree(fp) < 1) continue;  // se ignoran factores constantes
        result.push_back(Factor{detail::from_fmpz_poly(fp), static_cast<int>(fac->exp[i])});
    }

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(poly);

    std::sort(result.begin(), result.end(), factor_less);
    return result;
}

#else  // respaldo portable sin FLINT

namespace detail {

/// Divisores positivos de un entero no negativo.
inline std::vector<mpz_class> divisors(mpz_class n) {
    std::vector<mpz_class> out;
    if (n <= 0) return out;
    for (mpz_class k = 1; k * k <= n; ++k) {
        if (n % k == 0) {
            out.push_back(k);
            if (k * k != n) out.push_back(n / k);
        }
    }
    return out;
}

/// Factorizacion libre de cuadrados (Yun) de un polinomio monico sobre Q.
/// Devuelve pares (parte libre de cuadrados monica, multiplicidad).
inline std::vector<std::pair<QPoly, int>> square_free(const QPoly& P_in) {
    std::vector<std::pair<QPoly, int>> out;
    QPoly P = poly_monic(P_in);
    if (poly_degree(P) < 1) return out;

    QPoly dP = poly_deriv(P);
    QPoly a0 = poly_gcd(P, dP);
    QPoly b = poly_divmod(P, a0).first;
    QPoly c = poly_divmod(dP, a0).first;
    QPoly d = poly_sub(c, poly_deriv(b));
    int i = 1;
    while (poly_degree(b) > 0) {
        QPoly a = poly_gcd(b, d);
        if (poly_degree(a) > 0) out.push_back({poly_monic(a), i});
        QPoly b1 = poly_divmod(b, a).first;
        QPoly c1 = poly_divmod(d, a).first;
        d = poly_sub(c1, poly_deriv(b1));
        b = b1;
        ++i;
    }
    return out;
}

/// Extrae los factores lineales racionales de un polinomio primitivo entero.
/// Deja en 'rest' el resto sin raices racionales.
inline std::vector<ZPoly> linear_factors(ZPoly P, ZPoly& rest) {
    std::vector<ZPoly> factors;
    QPoly cur = poly_monic(to_qpoly(P));
    while (poly_degree(cur) >= 1) {
        ZPoly z = to_integer_primitive(cur);
        mpz_class a0 = z.back();
        mpz_class an = z.front();
        bool found = false;
        // Raiz cero.
        if (a0 == 0) {
            factors.push_back(ZPoly{mpz_class(1), mpz_class(0)});  // x
            cur = poly_divmod(cur, QPoly{mpq_class(1), mpq_class(0)}).first;
            continue;
        }
        std::vector<mpz_class> num = divisors(abs(a0));
        std::vector<mpz_class> den = divisors(abs(an));
        for (const auto& pn : num) {
            for (const auto& qd : den) {
                for (int sign : {1, -1}) {
                    mpq_class root(sign * pn, qd);
                    root.canonicalize();
                    if (poly_eval(cur, root) == 0) {
                        // factor (x - root) llevado a forma primitiva entera
                        QPoly lin{mpq_class(1), -root};
                        factors.push_back(to_integer_primitive(lin));
                        cur = poly_divmod(cur, lin).first;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (!found) break;
    }
    rest = (poly_degree(cur) >= 1) ? to_integer_primitive(cur) : ZPoly{};
    return factors;
}

}  // namespace detail

/// @brief Factoriza un polinomio racional (respaldo sin FLINT).
inline std::vector<Factor> factor_rational_poly(const QPoly& p) {
    std::vector<Factor> result;
    auto add = [&](const ZPoly& poly, int mult) {
        for (auto& f : result) {
            if (f.poly == poly) {
                f.mult += mult;
                return;
            }
        }
        result.push_back(Factor{poly, mult});
    };

    auto sf = detail::square_free(to_qpoly(to_integer_primitive(p)));
    for (const auto& part : sf) {
        const QPoly& s = part.first;
        const int mult = part.second;
        ZPoly rest;
        std::vector<ZPoly> lins = detail::linear_factors(to_integer_primitive(s), rest);
        for (const auto& lf : lins) add(lf, mult);
        if (!rest.empty() && (rest.size() > 1)) add(rest, mult);
    }

    std::sort(result.begin(), result.end(), factor_less);
    return result;
}

#endif  // DESCMATRAT_USE_FLINT

}  // namespace dmr
