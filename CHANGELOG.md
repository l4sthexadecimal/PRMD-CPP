# Historial de cambios

El formato sigue la idea de "Keep a Changelog" y versionado semantico.

## [0.1.0] - 2026-07-18

### Agregado
- Descomposicion racional primaria de matrices (`rational_primary_decomposition`):
  base de cambio `U` y `H = U^{-1} A U` diagonal por bloques.
- Forma de Frobenius con matriz de cambio `P` (`frobenius_form`), por el metodo
  constructivo del cociente, sin necesidad de factorizar.
- Puente de la descomposicion primaria a la forma de Frobenius
  (`frobenius_from_primary` e `invariant_factors_from_primary`).
- Variante fraction-free (eliminacion de Bareiss) como backend seleccionable y
  como funciones propias (`ff_rref`, `ff_solve`, `ff_nullspace`, `ff_rank`,
  `ff_determinant`).
- Aritmetica exacta de polinomios (gcd, lcm, division, forma primitiva) y
  factorizacion sobre Q con FLINT opcional mas respaldo portable.
- Biblioteca de solo cabecera y portable (Windows y Linux); solo GMP y Eigen
  son obligatorias.
- Documentacion Doxygen, pruebas con GoogleTest y CI de GitHub Actions.
- Garantia de seguridad ante hilos (reentrante, sin estado mutable compartido),
  documentada, para poder paralelizar por lotes con OpenMP u otras herramientas.
