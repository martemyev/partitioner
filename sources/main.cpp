#include "parameters.h"
#include "fine_mesh.h"
#include <iostream>

int main(int argc, char **argv)
{
  Parameters param(argc, argv);
  param.establish_environment();

#if defined(DEBUG)
  std::cout << param.print() << std::endl;
#endif

  FineMesh fmesh(&param);
  fmesh.read(param.MESH_FILE);
  fmesh.make_partition();
  fmesh.write(param.MESH_FILE_OUT);

  return 0;
}

