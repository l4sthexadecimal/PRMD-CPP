/// @file rational.hpp
/// @brief Tipos base para trabajar con matrices de racionales exactos.
///
/// El escalar exacto es mpq_class de GMP y el contenedor de matrices es Eigen.
/// Para que Eigen sepa manejar mpq_class hay que declararle sus rasgos
/// numericos (NumTraits); con eso el producto y la suma de matrices funcionan
/// de forma exacta, sin ningun redondeo.
///
/// El codigo es de solo cabecera y no usa nada especifico del sistema
/// operativo, asi que compila igual en Windows y en Linux.
#pragma once

#include <gmpxx.h>

#include <Eigen/Core>

namespace Eigen {

/// @brief Rasgos numericos minimos para usar mpq_class como escalar de Eigen.
///
/// No se define un epsilon real porque la aritmetica es exacta.
template <>
struct NumTraits<mpq_class> : GenericNumTraits<mpq_class> {
    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    typedef mpq_class Nested;
    typedef mpq_class Literal;

    static inline Real epsilon() { return mpq_class(0); }
    static inline Real dummy_precision() { return mpq_class(0); }
    static inline int digits10() { return 0; }

    enum {
        IsInteger = 0,
        IsSigned = 1,
        IsComplex = 0,
        RequireInitialization = 1,
        ReadCost = 6,
        AddCost = 150,
        MulCost = 100
    };
};

}  // namespace Eigen

/// @brief Espacio de nombres de la biblioteca DescMatRat.
namespace dmr {

/// Matriz densa de racionales exactos.
using QMat = Eigen::Matrix<mpq_class, Eigen::Dynamic, Eigen::Dynamic>;
/// Vector columna de racionales exactos.
using QVec = Eigen::Matrix<mpq_class, Eigen::Dynamic, 1>;

/// @brief Matriz identidad de tamano n.
inline QMat identity(int n) {
    QMat I = QMat::Zero(n, n);
    for (int i = 0; i < n; ++i) I(i, i) = mpq_class(1);
    return I;
}

/// @brief Vector canonico e_i de longitud n (todo cero salvo un 1 en la posicion i).
inline QVec unit_vector(int n, int i) {
    QVec e = QVec::Zero(n);
    e(i) = mpq_class(1);
    return e;
}

}  // namespace dmr
