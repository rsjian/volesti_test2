#include <iostream>
#include <fstream>

//#include "random.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
//#include "random_walks/boundary_rdhr_walk.hpp"
#include "random_walks/random_walks.hpp"
#include "sampling/sampling.hpp"
#include "SDPAFormatManager.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron <Point> SPECTRAHEDRON;
typedef BoostRandomNumberGenerator<boost::mt19937, double, 3> RNGType;

int main() {
    std::string fileName("../../examples/optimization_spectrahedra/data/sdp_n2m3.txt");

    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    // read the spectrahedron
    // open a stream to read the input file
    std::ifstream in;
    in.open(fileName, std::ifstream::in);

    // read the file
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron, objFunction);

    std::list<Point> points;
    RNGType rng(spectrahedron.getLMI().dimension()); // this class provides random numbers
	int walkLen = 5;
    int numpoints = 1000;
    Point initialPoint(spectrahedron.getLMI().dimension());
    int nburns = 0;
    uniform_sampling_boundary<BRDHRWalk>(points, spectrahedron, rng, walkLen, numpoints, initialPoint, nburns);

    // print sampled points
    for (Point point : points)
        point.print();

	return 0;
}
