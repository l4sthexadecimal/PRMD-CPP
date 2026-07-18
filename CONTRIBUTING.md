# Como contribuir

Gracias por tu interes en DescMatRat. Estas son las pautas basicas.

## Antes de enviar cambios

1. Compila y corre las pruebas:

   ```bash
   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
   cmake --build build --parallel
   ctest --test-dir build --output-on-failure
   ```

   Todo debe quedar en verde. La CI hace lo mismo en cada pull request.

2. Manten el codigo exacto. No introduzcas aritmetica de punto flotante en el
   nucleo de calculo; el escalar es `mpq_class`.

3. Si agregas una funcion nueva, agrega tambien su prueba en `tests/` y su
   documentacion Doxygen.

## Estilo

- C++20, solo cabecera. Nada especifico de un sistema operativo.
- Formatea con `clang-format` usando el `.clang-format` del repositorio.
- Comentarios claros y directos, sin adornos ni marcos de simbolos.

## Reportar problemas

Incluye un caso minimo reproducible: la matriz de entrada, la funcion llamada y
el resultado esperado frente al obtenido.
