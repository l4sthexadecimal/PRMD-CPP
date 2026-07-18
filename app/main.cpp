/// @file main.cpp
/// @brief Demo de la descomposicion racional primaria y de la forma de Frobenius.
///
/// Reproduce los casos de los scripts MATLAB (demo_rational_decomposition.m y
/// mattest.m) y verifica de forma exacta que:
///   - A*U = U*H para la descomposicion primaria,
///   - A*P = P*F para la forma de Frobenius,
///   - los factores invariantes del puente coinciden con los directos.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "descmatrat/descmatrat.hpp"

using namespace dmr;

namespace {

/// Imprime una matriz racional alineando por columnas.
void print_matrix(const std::string& name, const QMat& M) {
    std::vector<std::size_t> width(M.cols(), 0);
    std::vector<std::vector<std::string>> cells(M.rows(), std::vector<std::string>(M.cols()));
    for (int i = 0; i < M.rows(); ++i)
        for (int j = 0; j < M.cols(); ++j) {
            cells[i][j] = M(i, j).get_str();
            width[j] = std::max(width[j], cells[i][j].size());
        }
    std::cout << name << " =\n";
    for (int i = 0; i < M.rows(); ++i) {
        std::cout << "  [ ";
        for (int j = 0; j < M.cols(); ++j) {
            std::cout << std::setw(static_cast<int>(width[j])) << cells[i][j];
            if (j + 1 < M.cols()) std::cout << "  ";
        }
        std::cout << " ]\n";
    }
}

/// Construye una QMat a partir de una lista de filas.
QMat make_matrix(const std::vector<std::vector<long>>& rows) {
    const int n = static_cast<int>(rows.size());
    const int m = static_cast<int>(rows[0].size());
    QMat A(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) A(i, j) = mpq_class(rows[i][j]);
    return A;
}

void expect(bool ok, const std::string& what) {
    std::cout << what << " : " << (ok ? "OK" : "FALLO") << "\n";
    if (!ok) {
        std::cerr << "Verificacion fallida: " << what << "\n";
        std::exit(1);
    }
}

void run_case(const std::string& title, const QMat& A) {
    std::cout << "==== " << title << " ====\n";
    print_matrix("A", A);

    // Descomposicion primaria.
    PrimaryDecomposition d = rational_primary_decomposition(A);
    std::cout << "\npolinomio caracteristico: " << poly_to_string(d.char_poly) << "\n";
    std::cout << "factores irreducibles:\n";
    for (const FactorDiag& fd : d.diagnostics)
        std::cout << "  (" << fd.factor << ")^" << fd.multiplicity
                  << "   dim del bloque = " << fd.nullity << "\n";
    std::cout << "\n";
    print_matrix("U", d.U);
    std::cout << "\n";
    print_matrix("H (diagonal por bloques)", d.H);
    expect(A * d.U == d.U * d.H, "\nverificacion A*U == U*H");

    // Forma de Frobenius usando la descomposicion primaria como puente.
    FrobeniusForm fr = frobenius_from_primary(A);
    std::cout << "\nfactores invariantes:\n";
    for (const QPoly& f : fr.invariant_factors) std::cout << "  " << poly_to_string(f) << "\n";
    std::cout << "\n";
    print_matrix("P", fr.P);
    std::cout << "\n";
    print_matrix("F (forma de Frobenius)", fr.F);
    expect(A * fr.P == fr.P * fr.F, "\nverificacion A*P == P*F");

    std::cout << "\n";
}

}  // namespace

// Comprueba que el backend fraction-free da el mismo resultado exacto que el
// racional y que tambien reconstruye A*P = P*F.
void run_backend_check(const QMat& A) {
    std::cout << "==== Backend fraction-free (variante de Bareiss) ====\n";
    FrobeniusForm rat = frobenius_from_primary(A, Backend::Rational);
    FrobeniusForm ff = frobenius_from_primary(A, Backend::FractionFree);
    expect(rat.P == ff.P && rat.F == ff.F, "coincide con el backend racional (P y F)");
    expect(A * ff.P == ff.P * ff.F, "reconstruccion A*P == P*F (fraction-free)");
    std::cout << "\n";
}

int main() {
    std::cout << "Demo de descomposicion racional primaria y forma de Frobenius\n\n";

    run_case("Caso 3x3", make_matrix({{0, -1, 0}, {1, 0, 0}, {0, 0, 2}}));
    run_case("Caso 4x4",
             make_matrix({{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}}));
    QMat dense = make_matrix({{-1, -3, -3, -3}, {3, 5, 3, 3}, {-1, -1, 1, -1}, {1, 1, 1, 3}});
    run_case("Caso 4x4 denso", dense);

    run_backend_check(dense);

    std::cout << "Todos los casos verificados.\n";
    return 0;
}
