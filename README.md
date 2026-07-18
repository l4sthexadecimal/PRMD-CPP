# DescMatRat

Descomposicion racional de matrices con invariantes por bloques, en **aritmetica
exacta** sobre los numeros racionales.

<!-- Reemplaza OWNER/REPO por tu usuario y repositorio para que el badge funcione. -->
![CI](https://github.com/OWNER/REPO/actions/workflows/ci.yml/badge.svg)

DescMatRat toma una matriz racional cuadrada `A` y calcula, sin ningun
redondeo:

- su **descomposicion racional primaria**: una base `U` y `H = U^{-1} A U`
  diagonal por bloques, un bloque por componente primaria;
- su **forma de Frobenius** (forma racional canonica) `F` con la matriz de
  cambio `P`, tal que `P^{-1} A P = F`;
- el **puente** entre ambas: de la descomposicion primaria se obtienen los
  factores invariantes y con ellos la forma de Frobenius.

Es una biblioteca **de solo cabecera** y **portable** (Windows y Linux). Solo
depende de GMP y Eigen; FLINT es opcional.

## Por que

La forma de Frobenius y la descomposicion primaria son objetos clasicos del
algebra lineal (teoria de modulos sobre un dominio de ideales principales). En
una implementacion exacta interesan tres cosas que aqui se cuidan:

- reconstruir de forma explicita la base de cambio (`U`, `P`), no solo los
  invariantes polinomicos;
- mantener todo exacto, para poder verificar `A U = U H` y `A P = P F` sin
  tolerancias;
- controlar el crecimiento de coeficientes, para lo cual se ofrece una variante
  fraction-free (eliminacion de Bareiss).

## Un vistazo rapido

```cpp
#include "descmatrat/descmatrat.hpp"
using namespace dmr;

int main() {
    QMat A(3, 3);
    A << 0, -1, 0,
         1,  0, 0,
         0,  0, 2;

    // Descomposicion racional primaria: A*U = U*H, H diagonal por bloques.
    PrimaryDecomposition d = rational_primary_decomposition(A);

    // Forma de Frobenius usando la primaria como puente: A*P = P*F.
    FrobeniusForm fr = frobenius_from_primary(A);
    for (const QPoly& f : fr.invariant_factors)
        std::cout << poly_to_string(f) << "\n";   // factores invariantes
}
```

Para el ejemplo anterior el polinomio caracteristico es `x^3 - 2x^2 + x - 2` y,
como la matriz es ciclica, el unico factor invariante es ese mismo polinomio.

## Funciones principales

| Funcion | Que devuelve |
| --- | --- |
| `rational_primary_decomposition(A, backend)` | `U`, `H = U^{-1} A U`, factores, divisores elementales |
| `frobenius_form(A, backend)` | `P`, `F`, factores invariantes (construccion directa) |
| `frobenius_from_primary(A, backend)` | lo mismo, usando la primaria como puente |
| `invariant_factors_from_primary(pd)` | factores invariantes a partir de la primaria |

Cada una acepta un `Backend` opcional:

- `Backend::Rational` (por defecto): eliminacion racional directa.
- `Backend::FractionFree`: eliminacion de Bareiss. Trabaja con enteros y divide
  de forma exacta por el pivote anterior, lo que acota el crecimiento de
  numeradores y denominadores intermedios. Da exactamente el mismo resultado.

Tambien se exponen las primitivas fraction-free por separado: `ff_rref`,
`ff_solve`, `ff_nullspace`, `ff_rank` y `ff_determinant`.

## Dependencias

| Biblioteca | Obligatoria | Para que |
| --- | --- | --- |
| GMP (`mpq_class`) | si | escalar racional exacto |
| Eigen | si | contenedor de matrices |
| FLINT | no | factorizacion exacta y completa sobre Q |

Sin FLINT se usa un factorizador de respaldo incluido (libre de cuadrados mas
raices racionales), suficiente para los casos habituales. La forma de Frobenius
no necesita factorizar en ningun caso. La macro `DESCMATRAT_USE_FLINT` activa
FLINT y CMake la define sola si lo encuentra.

## Compilar y probar

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/demo
```

En Debian o Ubuntu las dependencias se instalan con:

```bash
sudo apt-get install cmake ninja-build g++ \
     libgmp-dev libeigen3-dev libflint-dev libgtest-dev
```

## Seguridad ante hilos (thread-safety)

La biblioteca es **reentrante** y **thread-safe**: las funciones no comparten
estado mutable, cada llamada solo toca sus argumentos y sus variables locales, y
no se fijan asignadores de memoria globales. Se puede llamar desde varios hilos a
la vez siempre que cada hilo use sus propios objetos. Un patron tipico es
resolver un lote de matrices independientes en paralelo:

```cpp
#include "descmatrat/descmatrat.hpp"
using namespace dmr;

std::vector<QMat> entradas = ...;
std::vector<FrobeniusForm> salidas(entradas.size());

#pragma omp parallel for
for (int k = 0; k < (int)entradas.size(); ++k)
    salidas[k] = frobenius_from_primary(entradas[k]);   // seguro: datos independientes
```

Notas:

- No se paraleliza dentro de una sola llamada; la paralelizacion queda a cargo
  de quien usa la biblioteca (OpenMP, `std::thread`, etc.). Por ahora sin GPU.
- Con FLINT activo, conviene que cada hilo llame a `flint_cleanup()` al terminar
  para liberar su cache local (higiene de memoria, no afecta la correctitud).

Hay una comparacion medida de la version serial frente a la paralela en el
proyecto de benchmark (carpeta hermana `benchmark/`).

## Documentacion

```bash
doxygen Doxyfile          # o: cmake --build build --target docs
# salida en docs/html/index.html
```

## Estructura del proyecto

```
include/descmatrat/   biblioteca de solo cabecera
  descmatrat.hpp        cabecera paraguas que incluye todo
  rational.hpp          tipos base (mpq_class, QMat) y rasgos para Eigen
  polynomial.hpp        polinomios: aritmetica, gcd, lcm, forma primitiva
  linalg.hpp            rref, nucleo, solve, companera, minimales, backend
  fraction_free.hpp     variante de Bareiss (ff_rref, ff_solve, ...)
  charpoly.hpp          polinomio caracteristico (Faddeev-LeVerrier)
  factor.hpp            factorizacion sobre Q (FLINT o respaldo)
  primary.hpp           descomposicion primaria y divisores elementales
  frobenius.hpp         forma de Frobenius y puente desde la primaria
app/main.cpp          demo con los casos de referencia
tests/                pruebas con GoogleTest
Doxyfile              configuracion de la documentacion
```

## Como se verifica

Las pruebas comprueban propiedades exactas, no aproximadas:

- `A U = U H` y `A P = P F` de forma exacta;
- divisibilidad de los factores invariantes y su producto igual al polinomio
  caracteristico;
- que el backend fraction-free da el mismo resultado que el racional (la forma
  escalonada reducida es unica) sobre decenas de matrices aleatorias.

## Licencia

MIT. Ver el archivo [LICENSE](LICENSE).
