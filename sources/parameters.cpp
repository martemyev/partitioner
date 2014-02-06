#include "parameters.h"
#include "config.h"
#include "auxiliary_functions.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include <iostream>

namespace po = boost::program_options;



Parameters::Parameters(int argc, char **argv)
{
  default_parameters(); // initialize all parameters
  if (argc > 1)
    read_from_command_line(argc, argv); // change some (or all) parameters from default to custom ones
}



void Parameters::default_parameters()
{
  MESH_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/meshes/";
  GEO_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/geometries/";

  RES_DIR = ""; // should be changed and based on some parameters
  INFO_FILE = "info.txt"; // should be added to RES_DIR after generating of the latter

  MESH_FILE = "mesh.msh";  // should be added to MESH_DIR after establishing of the latter (means that MESH_DIR can be changed from parameter file of command line)
  MESH_FILE_OUT = ""; // the name of the output mesh is based on the name of the input mesh and some parameters (for example, the number of partitions)

  N_PARTITIONS = 1; // total number of partitions
  N_PARTITIONS_X = 1; // number of partitions along x-axis
  N_PARTITIONS_Y = 1; // number of partitions along y-axis

  X_BEG = Y_BEG = 0.;
  X_END = Y_END = 1.;

  PRINT_INFO = 1; // print an information to console on each time step
}



void Parameters::read_from_command_line(int argc, char **argv)
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("meshfile", po::value<std::string>(),  std::string("name of mesh file (" + MESH_FILE + ")").c_str())
    ("meshdir",  po::value<std::string>(),  std::string("path to a directory with meshes (" + MESH_DIR + ")").c_str())
    ("inf",      po::value<bool>(),         std::string("whether we need to print some info during calculations (" + d2s(PRINT_INFO) + ")").c_str())
    ("x1",       po::value<double>(),       std::string("X_END (" + d2s(X_END) + ")").c_str())
    ("y1",       po::value<double>(),       std::string("Y_END (" + d2s(Y_END) + ")").c_str())
    ("np",       po::value<unsigned int>(), std::string("Total number of partitions (" + d2s(N_PARTITIONS) + ")").c_str())
    ("npx",      po::value<unsigned int>(), std::string("Number of partitions along x (" + d2s(N_PARTITIONS_X) + ")").c_str())
    ("npy",      po::value<unsigned int>(), std::string("Number of partitions along y (" + d2s(N_PARTITIONS_Y) + ")").c_str())
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    exit(0);
  }

  if (vm.count("meshfile"))
    MESH_FILE = vm["meshfile"].as<std::string>();

  if (vm.count("meshdir"))
    MESH_DIR = vm["meshdir"].as<std::string>();

  if (vm.count("inf"))
    PRINT_INFO = vm["inf"].as<bool>();

  if (vm.count("x1"))
    X_END = vm["x1"].as<double>();
  if (vm.count("y1"))
    Y_END = vm["y1"].as<double>();
  require(X_END > X_BEG && Y_END > Y_BEG, "Wrong limits of the computational domain");

  if (vm.count("np"))
  {
    require(!vm.count("npx") && !vm.count("npy"),
            "If the total number of partitions if defined you cannot define"
            " the number of partitions along direction (npx and/or npy arguments)");
    N_PARTITIONS = vm["np"].as<unsigned int>();
    N_PARTITIONS_X = int(sqrt(N_PARTITIONS));
    N_PARTITIONS_Y = int(N_PARTITIONS / N_PARTITIONS_X);
  }
  else
  {
    if (vm.count("npx"))
      N_PARTITIONS_X = vm["npx"].as<unsigned int>();
    if (vm.count("npy"))
      N_PARTITIONS_Y = vm["npy"].as<unsigned int>();
    N_PARTITIONS = N_PARTITIONS_X * N_PARTITIONS_Y;
  }
  require(N_PARTITIONS_X * N_PARTITIONS_Y == N_PARTITIONS,
          "The numbers of partitions along axes (npx = " +
          d2s(N_PARTITIONS_X) + ", npy = " + d2s(N_PARTITIONS_Y) +
          ") don't suit the total number of the partitions (np = " + d2s(N_PARTITIONS) + ")");
}



Parameters::~Parameters()
{ }



std::string Parameters::print() const
{
  std::string str = "list of parameters:\n";
  str += "dim = " + d2s(DIM) + "\n";
  str += "mesh file name = " + MESH_FILE + "\n";
  str += "domain = [" + d2s(X_BEG) + ", " + d2s(X_END) + "] x [" + d2s(Y_BEG) + ", " + d2s(Y_END) + "]\n";
  str += "print_info = " + d2s(PRINT_INFO) + "\n";
  str += "partitions = " + d2s(N_PARTITIONS) + "\n";
  str += "partitions_x = " + d2s(N_PARTITIONS_X) + "\n";
  str += "partitions_y = " + d2s(N_PARTITIONS_Y) + "\n";
  return str;
}



void Parameters::establish_environment()
{
  generate_paths(); // generate all necessary path to files and directories
  check_clean_dirs(); // check the existance and clearance of some directories
}



void Parameters::generate_paths()
{
  using namespace boost::filesystem;

  MESH_FILE = MESH_DIR + "/" + MESH_FILE; // full path to the mesh file
  path mesh_file(MESH_FILE);

  RES_DIR = RESULTS_DIR + "/"; // + mesh_file.stem().string() + "/";
  // "pa" means "partitioner", that means that the mesh was partitioned by this program.
  // it's done to distinguish different partitioners like "ch" (chaco) and "me" (metis)
  MESH_FILE_OUT = RES_DIR + "/" + mesh_file.stem().string() + "_pa_" + d2s(N_PARTITIONS) + ".msh";

  INFO_FILE = RES_DIR + "/" + INFO_FILE;
}



void Parameters::check_clean_dirs() const
{
  using namespace boost::filesystem;

  path top_res_dir(RESULTS_DIR); // top directory with results (from config.h, CMakeLists.txt)
  if (!exists(top_res_dir)) // if this top directory for results doesn't exist, we create it
    create_directory(top_res_dir);

  require(RES_DIR != "", "Directory for results has no name");
//  path cur_res_dir(RES_DIR); // current directory with results
//  if (exists(cur_res_dir) && is_directory(cur_res_dir)) // if this directory exists, we need to clean it up
//    remove_all(cur_res_dir); // remove all contents of the directory and the directory itself
//  create_directory(cur_res_dir); // now create empty directory
}
