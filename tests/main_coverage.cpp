#define CATCH_CONFIG_RUNNER // we will define main()
#include "catch.hpp"
#include <armadillo>

int main(int argc, char **argv) {
  std::cout << "Armadillo version: " << arma::arma_version::as_string() << '\n';

  const char *args[] = {"test_main_coverage",
                        "*"}; // We want to run all tests (also hidden ones!)
  int result = Catch::Session().run(2, args);
  // We ignore the result and return 0
  // This is needed so the code coverage works correctly.
  return 0;
}
