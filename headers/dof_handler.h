#ifndef DOF_HANDLER_H
#define DOF_HANDLER_H

#include "point.h"
#include "edge.h"
//#include "sym_csr_pattern.h"
#include <vector>


class FineMesh;
class Parameters;
class SymCSRPattern;


class DoFHandler
{
public:
            /**
             * Constructor
             */
  DoFHandler(FineMesh *fmesh);

            /**
             * Destructor
             */
  ~DoFHandler();

            /**
             * Distribute degrees of freedom according to fine mesh elements (connectivity of vertices),
             * order of basis functions (from param object), and type of finite element (CG, DG, etc).
             */
  void distribute_dofs(const Parameters &param);

            /**
             * Edge numeration
             */
  void numerate_edges();

            /**
             * Get the number of the degrees of freedom
             */
  unsigned int n_dofs() const;

            /**
             * Get the degree of freedom (its copy)
             * @param number - the serial number of the dof
             */
  Point dof(unsigned int number) const;

            /**
             * Get a constant reference for the whole list of dofs
             */
  const std::vector<Point>& dofs() const;

            /**
             * Const pointer to the fine mesh for access to mesh's functions
             */
  const FineMesh* fmesh() const;

            /**
             * Get the number of dof associated with mesh vertex
             * @param ver - the number of the vertex
             * @param num - the serial number of the dof
             */
  unsigned int vertices_dofs(unsigned int ver, unsigned int num) const;

            /**
             * Get the vector of dofs associated with one mesh vertex
             * @param ver - the number of the mesh vertex
             */
  const std::vector<unsigned int>& vertices_dofs(unsigned int ver) const;

  unsigned int n_dis_edges() const;
  unsigned int n_con_edges() const;
  Edge dis_edge(unsigned int num) const;
  Edge con_edge(unsigned int num) const;
  void boundary_dofs(const std::vector<unsigned int> &b_vertices,
                     std::vector<int> &bound_dofs) const;


private:
            /**
             * Finite triangular mesh which the dof handler is connected with
             */
  FineMesh *_fmesh;

            /**
             * Degrees of freedom
             */
  std::vector<Point> _dofs;

            /**
             * Edges for CG method
             */
  std::vector<Edge> _cg_edges;

            /**
             * Edges for DG method. They are called not '_dg_edges' to avoid confusion with '_cg_edges'.
             */
  std::vector<Edge> _discon_edges;

            /**
             * A map (an array of vectors) between mesh vertices and DG dofs
             */
  std::vector<unsigned int> *_vertices_dofs;

            /**
             * Copy constructor is private to prevent copying of dof handlers
             */
  DoFHandler(const DoFHandler&);

            /**
             * Copy assignment operator is also private to prevent copying
             */
  DoFHandler& operator=(const DoFHandler&);

            /**
             * Distribute degrees of freedom in case of first order basis functions
             */
  void distribute_first(const Parameters &param);

            /**
             * Distribute degrees of freedom in case of first order basis functions
             * for Discontinuous Galerkin method
             */
  void distribute_dg_first(const Parameters &param);
};



// =======================================
//
//
//
// =======================================

void associated_dg_edges(unsigned int ver_i,
                         unsigned int ver_j,
                         const std::vector<unsigned int> *vertices_dofs,
                         const SymCSRPattern &dof_pattern,
                         std::vector<unsigned int> &assoc_dg_edges);


#endif // DOF_HANDLER_H
