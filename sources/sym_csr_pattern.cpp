#include "sym_csr_pattern.h"
#include "triangle.h"
#include "auxiliary_functions.h"
#include "dof_handler.h"
#include "fine_mesh.h"
#include <set>

SymCSRPattern::SymCSRPattern()
{ }



SymCSRPattern::~SymCSRPattern()
{ }



void SymCSRPattern::make_sparse_format(const DoFHandler &dof_handler)
{
  clear(); // free some memory in the case if sparse format was already made

  // the number of rows of the matrix (and its order) of connectivity between degrees of freedom
  _order = dof_handler.n_dofs();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all dofs on them
  expect(dof_handler.fmesh()->n_triangles() != 0, "");
  for (int cell = 0; cell < dof_handler.fmesh()->n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh()->triangle(cell);

    expect(triangle.n_dofs() != 0, "");
    for (int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (int dj = 0; dj < triangle.n_dofs(); ++dj)
      {
        const unsigned int dof_j = triangle.dof(dj); // the number of the second degree of freedom

        // since it's a symmetric pattern we consider only lower triangle (without diagonal too).
        // insert the values in the corresponding places
        if (dof_i > dof_j)
          connect[dof_i].insert(dof_j);
      }
    }
  }

  // initialization of the pattern
  pattern_initialization(connect);

  // free the memory
  for (int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}



void SymCSRPattern::make_sparse_format(const FineMesh &fmesh)
{
  clear(); // free some memory in the case if sparse format was already made

  // the number of rows of the matrix (and its order) of connectivity between mesh vertices
  _order = fmesh.n_vertices();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all vertices on them
  expect(fmesh.n_triangles() != 0, "");
  for (int cell = 0; cell < fmesh.n_triangles(); ++cell)
  {
    Triangle triangle = fmesh.triangle(cell);

    for (int di = 0; di < Triangle::n_vertices; ++di)
    {
      const unsigned int ver_i = triangle.vertex(di); // the number of the first vertex
      for (int dj = 0; dj < Triangle::n_vertices; ++dj)
      {
        const unsigned int ver_j = triangle.vertex(dj); // the number of the second vertex

        // since it's a symmetric pattern we consider only lower triangle (without diagonal too).
        // insert the values in the corresponding places
        if (ver_i > ver_j)
          connect[ver_i].insert(ver_j);
      }
    }
  }

  // initialization of the pattern
  pattern_initialization(connect);

  // free the memory
  for (int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}



int SymCSRPattern::find(unsigned int num_row, unsigned num_col) const
{
  // we just add one important check
  expect(num_row > num_col, "For symmetric pattern the number of a row has to be bigger than the number of a column");
  // in the rest the function is the same as its parental variant
  return CSRPattern::find(num_row, num_col);
}
