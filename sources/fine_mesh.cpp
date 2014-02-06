#include "fine_mesh.h"
#include "auxiliary_functions.h"
#include <fstream>
#include <map>


FineMesh::FineMesh(Parameters *param)
  : _param(param)
{ }



FineMesh::~FineMesh()
{
  clear();
}



void FineMesh::clear()
{
  _vertices.clear();
  _triangles.clear();
  _lines.clear();
  _edges.clear();
  _boundary_vertices.clear();
}



void FineMesh::read(const std::string &file)
{
  std::ifstream in(file.c_str());
  require(in, "File " + file + " cannot be opened!");

  clear(); // if some mesh has been already read before, we delete it

  std::string str;
  in >> str; // the first string of Gmsh file is "$MeshFormat"
  expect(str == "$MeshFormat",
         "The first string of the Gmsh file " + file +
         " doesn't equal to \"$MeshFormat\". The actual string is \"" + str + "\"");

  // read the information about the mesh
  double version;
  int binary, dsize;
  in >> version >> binary >> dsize;
  expect(version == 2.2,
         "The version of Gmsh's mesh is unexpected (" + d2s(version) + ").");
  expect(dsize == sizeof(double),
         "The size of Gmsh's double (" + d2s(dsize) +\
         ") doesn't equal to size of double type (" + d2s(sizeof(double)) + ")");

  getline(in, str); // read some empty string

  // there is additional 1 (the number - one) in binary format
  if (binary)
  {
    int one;
    in.read(reinterpret_cast<char*>(&one), sizeof(int));
    require(one == 1,
            "The binary one (" + d2s(one) + ") doesn't equal to 1!");
  }

  // we make a map between serial number of the vertex and its number in the file.
  // it will help us when we create mesh elements
  std::map<unsigned int, unsigned int> vertices_map;

  // read lines of mesh file.
  // if we face specific keyword, we'll treat the section.
  while (in >> str)
  {
    if (str == "$Nodes") // read the mesh vertices
    {
      unsigned int n_vertices; // the number of all mesh vertices (that are saved in the file)
      in >> n_vertices; // read that number
      _vertices.resize(n_vertices); // allocate the memory for mesh vertices
      getline(in, str); // read some empty string

      unsigned int number; // the number of the vertex
      double coord[Point::n_coord]; // Cartesian coordinates of the vertex (Gmsh produces 3D mesh regardless its real dimension)
      double max_coord[Point::n_coord]; // the limits of the computational domain: maximum
      double min_coord[Point::n_coord]; // and minimum

      // read vertices
      for (unsigned int ver = 0; ver < n_vertices; ++ver)
      {
        if (binary) // binary format
        {
          in.read(reinterpret_cast<char*>(&number), sizeof(unsigned int)); // the number of each node
          in.read(reinterpret_cast<char*>(coord), Point::n_coord * sizeof(double)); // node coordinates
        }
        else // ASCII format
        {
          in >> number;
          for (unsigned int i = 0; i < Point::n_coord; ++i)
            in >> coord[i];
          if (ver == 0) // for the first vertex
          {
            for (unsigned int i = 0; i < Point::n_coord; ++i)
              max_coord[i] = min_coord[i] = coord[i]; // initialization of the max and min coordinates
          }
          else // for the other vertices
          {
            for (unsigned int i = 0; i < Point::n_coord; ++i)
            {
              max_coord[i] = std::max(max_coord[i], coord[i]); // searching max and min coordinates
              min_coord[i] = std::min(min_coord[i], coord[i]);
            }
          }

        }
        _vertices[ver] = Point(coord); // save the vertex
        vertices_map[number] = ver; // add the number of vertex to the map
      }
      _max_coord = Point(max_coord); // this point can not be one of the mesh vertices if the domain has curvilinear boundaries
      _min_coord = Point(min_coord); // the same as above is true

      expect(n_vertices == vertices_map.size(),
             "Vertices numbers are not unique: n_vertices = " + d2s(n_vertices) +
             " vertices_map.size() = " + d2s(vertices_map.size()));

    } // read the vertices

    else if (str == "$Elements") // read the mesh elements
    {
      unsigned int n_elements; // the number of mesh elements
      in >> n_elements; // read that number
      getline(in, str); // empty string

      unsigned int number; // the number of the element [1, nElements]
      unsigned int el_type; // the type of the element (1 - line, 2 - triangle, etc)
      unsigned int n_tags; // the number of tags describing the element
      unsigned int phys_domain; // the physical domain where the element takes place
      unsigned int elem_domain; // the elementary domain where the element takes place
      unsigned int n_partitions; // the number of partitions in which the element takes place
      int partition; // the partition which the element belongs to
      std::vector<unsigned int> ghost_cells; // "ghost cells" are other partitions which this element is connected with

      // the map between the type of the element,
      // and the number of nodes that describe it
      std::map<unsigned int, unsigned int> type_nodes;
      type_nodes[1] = 2; // 2-nodes line
      type_nodes[2] = 3; // 3-nodes triangle
      type_nodes[3] = 4; // 4-nodes quadrangle
      type_nodes[4] = 4; // 4-nodes tetrahedron
      type_nodes[5] = 8; // 8-nodes hexahedron
      type_nodes[15]= 1; // 1-node point

      if (binary) // binary format
      {
        require(false, "Binary 2.0 format of Gmsh mesh files is not supported");
      } // binary format

      else // ASCII format
      {
        for (int el = 0; el < n_elements; ++el)
        {
          in >> number >> el_type >> n_tags;
          std::vector<int> data(n_tags); // allocate the memory for some data
          for (unsigned int i = 0; i < n_tags; ++i) // read this information
            in >> data[i];
          phys_domain = (n_tags > 0) ? data[0] : 0; // physical domain - the most important value
          elem_domain = (n_tags > 1) ? data[1] : 0; // elementary domain
          if (n_tags > 2)
          {
            // the number of partitions where this elements takes place
            n_partitions = data[2];
            expect(n_partitions >= 1,
                   "The number of tags is more than 2. That means that we have partitions. But the number of partitions is "
                   + d2s(n_partitions));
            // the partition which the element belongs to
            partition = data[3];
            // "ghost cells"
            if (n_partitions > 1) // if the element is on the boundary between the partitions, it is described by "ghost cells" as well
            {
              ghost_cells.resize(n_partitions - 1);
              for (int gc = 0; gc < n_partitions - 1; ++gc)
                ghost_cells[gc] = -data[4 + gc]; // 'minus' since ghost cells are described by number of partition with the negative sing
            }
          }

          data.clear(); // clear the memory

          // how many vertices (nodes) describe the element
          std::map<unsigned int, unsigned int>::const_iterator el_type_iter =
              type_nodes.find(el_type);

          require(el_type_iter != type_nodes.end(),
                  "This type of the Gmsh's element (" + d2s(el_type) +
                  ") in the mesh file \"" + file + "\" is unknown!");

          const unsigned int n_elem_nodes = el_type_iter->second; // the number of nodes
          std::vector<unsigned int> nodes(n_elem_nodes); // allocate memory for nodes
          for (unsigned int i = 0; i < n_elem_nodes; ++i)
          {
            in >> nodes[i]; // read the numbers of nodes
            // vertices can be numerated not sequentially (or at least not from 0)
            nodes[i] = vertices_map.find(nodes[i])->second;
          }

          // add new element in the list
          //MeshElement *new_element;
          switch (el_type)
          {
//          case 15: // 1-node point
//            points.push_back(MeshElement_ptr(new PhysPoint(nodes, phys_domain)));
//            break;
          case 1: // 2-nodes line
            _lines.push_back(Line(nodes, phys_domain));
            //new_element = new Line(nodes, phys_domain);
            break;
          case 2: // 3-nodes triangle
            _triangles.push_back(Triangle(nodes, _vertices, phys_domain, partition, ghost_cells));
            //new_element = new Triangle(nodes, phys_domain);
            break;
//          case 3: // 4-nodes quadrangle
//            quadrangles.push_back(MeshElement_ptr(new Quadrangle(nodes, phys_domain)));
//            //new_element = new Quadrangle(nodes, phys_domain);
//            break;
//          case 4: //4-nodes tetrahedron
//            tetrahedra.push_back(MeshElement_ptr(new Tetrahedron(nodes, phys_domain)));
//            //new_element = new Tetrahedron(nodes, phys_domain);
//            break;
//          case 5: // 8-nodes hexahedron
//            hexahedra.push_back(MeshElement_ptr(new Hexahedron(nodes, phys_domain)));
//            //new_element = new Hexahedron(nodes, phys_domain);
//            break;
          default:
            require(false,
                    "Unknown type of the Gmsh's element (" + d2s(el_type) +
                    ") in the file " + file + "!");
          }

          ghost_cells.clear();
          nodes.clear();

          //elements.push_back(new_element);
          //delete new_element;
        }

        // check some expectations
        expect(number == n_elements,
               "The number of the last read Gmsh's element (" + d2s(number) +\
               ") is not equal to the amount of all elements in the mesh (" + d2s(n_elements) + ")!");

      } // ASCII format

      // requirements after reading elements
      require(!_triangles.empty(),
             "There are no any 2D or 3D elements in the mesh!");
      //require(!_lines.empty(),
      //       "There are no boundary lines in the mesh!");

    } // read the elements
  } // read the mesh file

  in.close(); // close the file

  // generate the list of boundary vertices
  boundary_vertices_initialization();
}



void FineMesh::boundary_vertices_initialization()
{
  if (_lines.size() == 0) // there are no boundary lines in the mesh
  {
    // in this case the domain should be rectangular.
    const double tol = 1e-14;
    require(fabs(_min_coord.coord(0) - _param->X_BEG) < tol &&
            fabs(_min_coord.coord(1) - _param->Y_BEG) < tol &&
            fabs(_max_coord.coord(0) - _param->X_END) < tol &&
            fabs(_max_coord.coord(1) - _param->Y_END) < tol,
            "Computational domain defined by parameters (X_BEG, Y_BEG, etc) is not equal to the domain read from the mesh file");

    // find nodes that lie on the boundary of the computational domain
    for (int i = 0; i < _vertices.size(); ++i)
    {
      const double x = _vertices[i].coord(0);
      const double y = _vertices[i].coord(1);
      if ((fabs(x - _min_coord.coord(0)) < tol ||
           fabs(y - _min_coord.coord(1)) < tol ||
           fabs(x - _max_coord.coord(0)) < tol ||
           fabs(y - _max_coord.coord(1)) < tol) &&
          find(_boundary_vertices.begin(), _boundary_vertices.end(), i) == _boundary_vertices.end())
        _boundary_vertices.push_back(i);
    }
  }
  else // there are the boundary lines in the mesh
  {
    for (int lin = 0; lin < _lines.size(); ++lin)
    {
      for (int ver = 0; ver < Line::n_vertices; ++ver)
      {
        const int node = _lines[lin].vertex(ver);
        if (find(_boundary_vertices.begin(), _boundary_vertices.end(), node) == _boundary_vertices.end())
          _boundary_vertices.push_back(node);
      }
    }
  }
}



unsigned int FineMesh::n_vertices() const
{
  return _vertices.size();
}



Point FineMesh::vertex(unsigned int number) const
{
  expect(number >= 0 && number < _vertices.size(), "Incorrect input parameter");
  return _vertices[number];
}



const std::vector<Point>& FineMesh::vertices() const
{
  return _vertices;
}



unsigned int FineMesh::n_triangles() const
{
  return _triangles.size();
}



Triangle FineMesh::triangle(int number) const
{
  expect(number >= 0 && number < _triangles.size(),
         "Incorrect number of a triangle (" + d2s(number) +
         "). It should be in a range [0, " + d2s(_triangles.size()) + ")");
  return _triangles[number];
}



Triangle* FineMesh::triangle_orig(unsigned int number)
{
  expect(number >= 0 && number < _triangles.size(), "Incorrect input parameter");
  return &_triangles[number];
}



Point FineMesh::max_coord() const
{
  return _max_coord;
}



Point FineMesh::min_coord() const
{
  return _min_coord;
}



unsigned int FineMesh::n_lines() const
{
  return _lines.size();
}



Line FineMesh::line(unsigned int number) const
{
  expect(number >= 0 && number < _lines.size(), "Incorrect input parameter");
  return _lines[number];
}



std::vector<unsigned int> const& FineMesh::boundary_vertices() const
{
  return _boundary_vertices;
}



void FineMesh::write(const std::string &filename) const
{
  std::ofstream out(filename); // oopen the file for writing
  require(out, "The file " + filename + " cannot be opened");

  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  out << "$Nodes\n";
  out << _vertices.size() << "\n";
  for (int i = 0; i < _vertices.size(); ++i)
  {
    out << i + 1;
    for (int j = 0; j < Point::n_coord; ++j)
      out << " " << _vertices[i].coord(j);
    out << "\n";
  }
  out << "$EndNodes\n";
  out << "$Elements\n";
  out << _triangles.size() + _lines.size() << "\n";
  for (int i = 0; i < _triangles.size(); ++i)
  {
    out << i + 1 << " "
        << Triangle::gmsh_el_type << " "
        << 4 + _triangles[i].n_ghost_cells() << " "
        << _triangles[i].material_id() << " "
        << _triangles[i].material_id() << " "
        << 1 + _triangles[i].n_ghost_cells() << " "
        << _triangles[i].partition_id() + 1; // partition and ghost cells are numbered from 0, but they should be from 1 for Gmsh
    for (int j = 0; j < _triangles[i].n_ghost_cells(); ++j)
      out << " -" << _triangles[i].ghost_cell(j) + 1; // ghost cells should have negative number, but they are saved as positive numbers
    for (int j = 0; j < Triangle::n_vertices; ++j)
      out << " " << _triangles[i].vertex(j) + 1;
    out << "\n";
  }
  for (int i = 0; i < _lines.size(); ++i)
  {
    out << _triangles.size() + i + 1 << " "
        << Line::gmsh_el_type << " 2 "
        << _lines[i].material_id() << " "
        << _lines[i].material_id();
    for (int j = 0; j < Line::n_vertices; ++j)
      out << " " << _lines[i].vertex(j) + 1;
    out << "\n";
  }

  out.close();
}



void FineMesh::make_partition()
{
  // we are trying to divide the domain on equal rectangular partitions
  // therefore the size of these partitions is the same
  const double hx = (_max_coord.coord(0) - _min_coord.coord(0)) / _param->N_PARTITIONS_X;
  const double hy = (_max_coord.coord(1) - _min_coord.coord(1)) / _param->N_PARTITIONS_Y;

  // distribute all triangles by partitions
  for (int tr = 0; tr < _triangles.size(); ++tr)
  {
    double xcen = 0, ycen = 0; // center of the triangle
    double xver[Triangle::n_vertices];
    double yver[Triangle::n_vertices];
    for (int v = 0; v < Triangle::n_vertices; ++v)
    {
      xcen += _vertices[_triangles[tr].vertex(v)].coord(0);
      ycen += _vertices[_triangles[tr].vertex(v)].coord(1);
      xver[v] = _vertices[_triangles[tr].vertex(v)].coord(0);
      yver[v] = _vertices[_triangles[tr].vertex(v)].coord(1);
    }
    xcen /= Triangle::n_vertices;
    ycen /= Triangle::n_vertices;

    for (int v = 0; v < Triangle::n_vertices; ++v)
    {
      xver[v] = xcen;//0.5 * (xver[v] + xcen);
      yver[v] = ycen;//0.5 * (yver[v] + ycen);
    }

    int part_x = -1;
    for (int i = 0; i < _param->N_PARTITIONS_X && part_x == -1; ++i)
    {
      const double x0 = _min_coord.coord(0) + i * hx;
      const double x1 = _min_coord.coord(0) + (i + 1) * hx;
      int nv = 0;
      for (int v = 0; v < Triangle::n_vertices; ++v)
        if (xver[v] >= x0 && xver[v] <= x1)
          ++nv;
      if (nv == Triangle::n_vertices)
        part_x = i;
      //      if (xcen >= _min_coord.coord(0) + i * hx &&
      //          xcen <= _min_coord.coord(0) + (i + 1) * hx)
    }
    expect(part_x != -1, "x partition for triangle " + d2s(tr) + " is not found");

    int part_y = -1;
    for (int i = 0; i < _param->N_PARTITIONS_Y && part_y == -1; ++i)
    {
      const double y0 = _min_coord.coord(1) + i * hy;
      const double y1 = _min_coord.coord(1) + (i + 1) * hy;
      int nv = 0;
      for (int v = 0; v < Triangle::n_vertices; ++v)
        if (yver[v] >= y0 && yver[v] <= y1)
          ++nv;
      if (nv == Triangle::n_vertices)
        part_y = i;
      //      if (xcen >= _min_coord.coord(0) + i * hx &&
      //          xcen <= _min_coord.coord(0) + (i + 1) * hx)
    }
    expect(part_y != -1, "y partition for triangle " + d2s(tr) + " is not found");

    _triangles[tr].partition_id(part_y * _param->N_PARTITIONS_X + part_x);

  } // all triangles

  // now we need to make a connection between the triangles,
  // i.e. to fill the list of ghost cells for each triangle.
  // to do that we'll compare the vertices of all triangles
  for (int tr_first = 0; tr_first < _triangles.size(); ++tr_first)
  {
    const int partition_first = _triangles[tr_first].partition_id(); // partition of the triangle
    expect(partition_first != -1, "Partition of the triangle " + d2s(tr_first) + " is not initialized");
    std::vector<unsigned int> ghost_cells; // ghost cells (connections to other partitions) of this triangle

    for (int v = 0; v < Triangle::n_vertices; ++v) // consider all vertices of the triangle
    {
      const unsigned int vertex_interest = _triangles[tr_first].vertex(v); // vertex of interest

      // now we need to find this vertex in all other triangles
      for (int tr_second = 0; tr_second < _triangles.size(); ++tr_second)
      {
        const int partition_second = _triangles[tr_second].partition_id(); // partition of the second triangle
        expect(partition_second != -1, "Partition of the triangle " + d2s(tr_second) + " is not initialized");

        // if the second triangle is not our first triangle,
        // the partition of this new triangle is not the same as the partition of our first triangle
        // and this second triangle contains our vertex of interest
        if (tr_second != tr_first &&
            partition_second != partition_first &&
            _triangles[tr_second].contains_vertex(vertex_interest))
        {
          if (find(ghost_cells.begin(), ghost_cells.end(), partition_second) == ghost_cells.end())
            ghost_cells.push_back(partition_second); // add this partition to the list of the ghost cells
        }
      } // all triangles again
    } // all vertices

    if (!ghost_cells.empty()) // if there are any connections
    {
      _triangles[tr_first].ghost_cells(ghost_cells); // set the ghost cells of the triangle
      ghost_cells.clear();
    }

  } // all triangles
}
