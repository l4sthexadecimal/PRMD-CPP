/// @file descmatrat.hpp
/// @brief Cabecera paraguas: incluye toda la biblioteca DescMatRat.
///
/// Con incluir este archivo se dispone de todas las funciones. La biblioteca es
/// de solo cabecera y portable (Windows y Linux); solo depende de GMP y Eigen.
/// FLINT es opcional y solo mejora la factorizacion de polinomios (definir la
/// macro DESCMATRAT_USE_FLINT para activarlo).
///
/// @mainpage DescMatRat
///
/// Descomposicion racional de matrices con invariantes por bloques, en
/// aritmetica exacta sobre los racionales.
///
/// Funciones principales:
///   - dmr::rational_primary_decomposition : descomposicion racional primaria
///     de una matriz A. Devuelve la base U y H = U^{-1} A U diagonal por bloques.
///   - dmr::frobenius_form : forma de Frobenius (forma racional canonica) con su
///     matriz de cambio P. Construccion directa, no necesita factorizar.
///   - dmr::frobenius_from_primary : misma forma de Frobenius, pero obtenida
///     usando la descomposicion primaria como puente hacia los factores
///     invariantes.
///
/// Las tres aceptan un parametro dmr::Backend para elegir el algebra lineal:
/// racional directa o fraction-free (eliminacion de Bareiss), esta ultima util
/// para acotar el crecimiento de coeficientes. Ambas dan el mismo resultado
/// exacto.
///
/// Todo el calculo es exacto: el escalar es mpq_class (GMP) y no interviene en
/// ningun momento la aritmetica de punto flotante.
///
/// @section threadsafety Seguridad ante hilos
///
/// Las funciones son reentrantes y no comparten estado mutable: cada llamada
/// trabaja solo sobre sus argumentos y sus variables locales. Por eso es seguro
/// llamarlas desde varios hilos a la vez mientras cada hilo use sus propios
/// objetos (por ejemplo, resolver un lote de matrices independientes con
/// OpenMP o std::thread). GMP opera de forma segura sobre objetos distintos y
/// la biblioteca no fija asignadores de memoria globales.
///
/// Si se activa FLINT (macro DESCMATRAT_USE_FLINT), conviene que cada hilo
/// llame a flint_cleanup() al terminar para liberar su cache local; esto es
/// higiene de memoria, no afecta la correctitud. No se paraleliza dentro de una
/// sola llamada: la paralelizacion queda en manos de quien usa la biblioteca.
#pragma once

#include "descmatrat/charpoly.hpp"
#include "descmatrat/factor.hpp"
#include "descmatrat/fraction_free.hpp"
#include "descmatrat/frobenius.hpp"
#include "descmatrat/linalg.hpp"
#include "descmatrat/polynomial.hpp"
#include "descmatrat/primary.hpp"
#include "descmatrat/rational.hpp"
