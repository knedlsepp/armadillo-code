set (ARMADILLO_TEST_SOURCES attributes.cpp
                            bounds.cpp
                            catch.hpp
                            decomp_eig_gen.cpp
                            decomp_eig_sym.cpp
                            expr_elem.cpp
                            expr_misc.cpp
                            fn_abs.cpp
                            fn_accu.cpp
                            fn_all.cpp
                            fn_any.cpp
                            fn_as_scalar.cpp
                            fn_clamp.cpp
                            fn_cond.cpp
                            fn_conj.cpp
                            fn_conv.cpp
                            fn_conv_to.cpp
                            fn_cor.cpp
                            fn_cov.cpp
                            fn_cross.cpp
                            fn_cumprod.cpp
                            fn_cumsum.cpp
                            fn_det.cpp
                            fn_diagmat.cpp
                            fn_diagvec.cpp
                            fn_diff.cpp
                            fn_dot.cpp
                            fn_eps.cpp
                            fn_expmat.cpp
                            fn_find.cpp
                            fn_find_finite.cpp
                            fn_find_nonfinite.cpp
                            fn_find_unique.cpp
                            fn_flip.cpp
                            fn_interp1.cpp
                            fn_is_finite.cpp
                            fn_sum.cpp
                            fn_symmat.cpp
                            fn_trace.cpp
                            fn_trans.cpp
                            fn_vectorise.cpp
                            gen_linspace.cpp
                            gen_ones.cpp
                            gen_randu.cpp
                            gen_zeros.cpp
                            init_auxmem.cpp
                            init_fill.cpp
                            init_misc.cpp
                            instantiation.cpp
                            mat_minus.cpp
                            mat_mul_cx.cpp
                            mat_mul_real.cpp
                            mat_neg.cpp
                            mat_plus.cpp
                            )
add_executable(armadillo_unit_tests main.cpp ${ARMADILLO_TEST_SOURCES})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# set_property(TARGET armadillo_unit_tests PROPERTY CXX_STANDARD 11)
target_link_libraries(armadillo_unit_tests armadillo)
add_test(NAME armadillo_unit_tests COMMAND armadillo_unit_tests -s --force-colour)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose  )

# Add code-coverage summary for unit tests
option(ARMADILLO_COVERAGE "Enable coverage report analysis" OFF)
# Code-coverage in case of gcc-debug build
if((${ARMADILLO_COVERAGE} MATCHES ON) )
  if ((CMAKE_BUILD_TYPE MATCHES Debug) AND ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 --coverage -fno-inline -fno-inline-small-functions -fno-default-inline")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 --coverage -fno-inline -fno-inline-small-functions -fno-default-inline")
    set(CMAKE_EXE_LINKER_FLAGS "--coverage")
    add_executable(armadillo_unit_tests_coverage main_coverage.cpp ${ARMADILLO_TEST_SOURCES})
    target_link_libraries(armadillo_unit_tests_coverage armadillo)
    add_custom_target(coverage COMMAND armadillo_unit_tests_coverage
                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                               COMMENT "Running code coverage."
                     )
  else()
    message( FATAL_ERROR "Code coverage only implemented for debug builds using GCC")
  endif()
endif()
