#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "config.h"
#include "boost/filesystem.hpp"
#include <vector>
#include <string>


class Parameters
{
public:
  static const int FE_ORDER = 1;

            /**
             * Dimension of the problem.
             * The whole program is designed to solve the problem of specific dimension.
             */
  static const int DIM = 2;

            /**
             * The number of the domain that characterizes the inclusions
             */
  static const int INCL_DOMAIN = 11;

            /**
             * The limits of the 2D computational domain.
             * The points (X_BEG, Y_BEG) and (X_END, Y_END) are the mesh nodes only if the domain is rectangular.
             * In case of domain with curvilinear boundaries, these limits just show the maximal possible range of coordinates
             */
  double X_BEG, X_END, Y_BEG, Y_END;

            /**
             * Total number of rectangular-shape partitions
             */
  unsigned int N_PARTITIONS;

            /**
             * Number of partitions along x-axis
             */
  unsigned int N_PARTITIONS_X;

            /**
             * Number of partitions along y-axis
             */
  unsigned int N_PARTITIONS_Y;

            /**
             * The path to the directory where the meshes (.msh-files) are
             */
  std::string MESH_DIR;

            /**
             * The path to the directory where the geometry file (.geo-files) are
             */
  std::string GEO_DIR;

            /**
             * The path to the directory where all the results of the program will be kept
             */
  std::string RES_DIR;

            /**
             * The name of the file where the information about simulation will be kept.
             * This file is in the RES_DIR directory
             */
  std::string INFO_FILE;

            /**
             * The name of the file with a mesh. This is an input mesh.
             * This mesh file should be in the MESH_DIR directory
             */
  std::string MESH_FILE;

            /**
             * The name of the file with a mesh with partitions. This is an output mesh.
             * This mesh will be in the RES_DIR directory after program's termination.
             */
  std::string MESH_FILE_OUT;

            /**
             * Whether we need to print some info on the screen during the calculations
             */
  bool PRINT_INFO;

            /**
             * Constructor
             * @param argc - the number of command line arguments (+1 - the first argument is the name of the executable by default)
             * @param argv - the command line arguments themselves
             */
  Parameters(int argc = 0, char **argv = 0); // constructor

            /**
             * Destructor
             */
  ~Parameters();

            /**
             * Get the string with readable represantation of all (probably) parameters of the problem.
             * This string can be then output in file or standard output stream.
             */
  std::string print() const;

            /**
             * To make some useful things for calculations.
             * For example, create the structure of output directories, files, etc.
             * This procedure was separated from Parameter constructor to
             * avoid time-consuming operations with filesystem in default parameters objects
             * (that are used, for example, in testing procedures).
             * This function has to be called before real work!
             */
  void establish_environment();


private: //======================= PRIVATE =========================

            /**
             * Copy constructor. It's private to protect copying parameters, because it would be weird.
             * That's also useful to prevent copying objects of the Parameter class as function argument.
             * Therefore we always have to sent parameter object as a const reference or pointer
             */
  Parameters(const Parameters&);

            /**
             * Operator of copy assignment is private for the same reason as the copy constructor
             */
  const Parameters& operator =(const Parameters&);

            /**
             * Initialize default parameters
             */
  void default_parameters();

            /**
             * Read parameters from the command line
             * @param argc - the number of the command line arguments
             * @param argv - the command line arguments
             */
  void read_from_command_line(int argc, char **argv);

            /**
             * Generate the paths to some files and directories.
             * Some of the paths are based on the parameters, therefore
             * this function should be invoked after establishing all parameters
             */
  void generate_paths();

            /**
             * Check that the directories we are going to use for output data exist and empty
             */
  void check_clean_dirs() const;
};


#endif // PARAMETERS_H
