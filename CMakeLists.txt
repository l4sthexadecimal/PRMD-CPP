cmake_minimum_required(VERSION 3.16)
project(descmatrat LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Warnings activados, como pide la guia del proyecto.
add_compile_options(-Wall -Wextra -Wpedantic)

# --- Dependencias obligatorias -------------------------------------------

# Eigen (solo cabeceras).
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# GMP y su interfaz C++ (gmpxx). No trae config de CMake, se busca a mano.
find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h)
find_library(GMPXX_LIBRARY NAMES gmpxx)
find_library(GMP_LIBRARY NAMES gmp)
if(NOT GMPXX_INCLUDE_DIR OR NOT GMPXX_LIBRARY OR NOT GMP_LIBRARY)
    message(FATAL_ERROR "No se encontro GMP/gmpxx")
endif()

# --- Biblioteca de solo cabecera -----------------------------------------

add_library(descmatrat INTERFACE)
target_include_directories(descmatrat INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${GMPXX_INCLUDE_DIR}
)
target_link_libraries(descmatrat INTERFACE
    Eigen3::Eigen
    ${GMPXX_LIBRARY}
    ${GMP_LIBRARY}
)

# --- Dependencia opcional: FLINT -----------------------------------------
# Mejora la factorizacion de polinomios sobre Q. Si no esta, se usa el
# factorizador de respaldo incluido en la biblioteca.
find_path(FLINT_INCLUDE_DIR NAMES flint/fmpz.h)
find_library(FLINT_LIBRARY NAMES flint)
if(FLINT_INCLUDE_DIR AND FLINT_LIBRARY)
    message(STATUS "FLINT encontrado: factorizacion exacta completa activada")
    target_compile_definitions(descmatrat INTERFACE DESCMATRAT_USE_FLINT)
    target_include_directories(descmatrat INTERFACE ${FLINT_INCLUDE_DIR})
    target_link_libraries(descmatrat INTERFACE ${FLINT_LIBRARY})
else()
    message(STATUS "FLINT no encontrado: se usa el factorizador de respaldo")
endif()

# --- Ejecutable demo ------------------------------------------------------

add_executable(demo app/main.cpp)
target_link_libraries(demo PRIVATE descmatrat)

# --- Pruebas --------------------------------------------------------------

enable_testing()
find_package(GTest)
if(GTest_FOUND)
    add_executable(unit_tests
        tests/test_polynomial.cpp
        tests/test_linalg.cpp
        tests/test_charpoly.cpp
        tests/test_factor.cpp
        tests/test_primary.cpp
        tests/test_frobenius.cpp
        tests/test_fraction_free.cpp
    )
    target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_link_libraries(unit_tests PRIVATE descmatrat GTest::gtest GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(unit_tests)
else()
    message(WARNING "GTest no encontrado: se omiten las pruebas unitarias")
endif()

# --- Documentacion Doxygen (opcional) ------------------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generando documentacion con Doxygen"
        VERBATIM)
endif()
