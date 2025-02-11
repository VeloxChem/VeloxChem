#include <iostream>
#include <cstdlib>

#include <xtb.h>

int
main()
{
  auto ver = xtb_getAPIVersion();
  std::cout << "API version: " << ver << std::endl;

  return EXIT_SUCCESS;
}
