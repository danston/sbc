// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// Subdividing barycentric coordinates.

// Local includes.
#include "./coords/LocalR2.hpp"
#include "./coords/HarmonicR2.hpp"
#include "./coords/MeanValueR2.hpp"

#include "./extra/MeshR2.hpp"
#include "./extra/TriangulatorR2.hpp"
#include "./extra/QuadSubdivisionR2.hpp"
#include "./extra/TriangleSubdivisionR2.hpp"

using namespace gbc;

// Enumeration with different types of barycentric coordinates.
enum CoordinateType { MV, HM, LC };

// Create initial triangle mesh.
void createTriangleMesh(
	const double edgeLength, 
	const std::vector<VertexR2> &poly,
	std::vector<Face> &tf, std::vector<VertexR2> &tp) {

    // Refine the polygon to create regular mesh.
    std::vector<VertexR2> refined;
    const size_t n = poly.size();

    for (size_t i = 0; i < n; ++i) {
        refined.push_back(poly[i]);

        const size_t ip = (i + 1) % n;
        const size_t numS = ceil((poly[ip] - poly[i]).length() / edgeLength);

        for (size_t j = 1; j < numS; ++j) {

            VertexR2 vert = poly[i] + (double(j) / double(numS)) * (poly[ip] - poly[i]);
            refined.push_back(vert);
        }
    }

    // Create mesh.
    TriangulatorR2 tri(refined, edgeLength, true);
    tri.setPlanarGraph(true);

    tri.generate(tp, tf);
}

// Set default flags to the points of the mesh.
void setDefaultFlags(std::vector<VertexR2> &tp) {

	const size_t numP = tp.size();
	for (size_t i = 0; i < numP; ++i) {

		tp[i].val = 0;
		tp[i].out = -1;
		tp[i].alpha = -1.0;
		tp[i].type = INTERIOR;
	}
}

// Compute initial coordinates.
void computeInitialCoordinates(
	const double edgeLength,
	const std::vector<VertexR2> &poly,
	const std::vector<Face> &tf, std::vector<VertexR2> &tp,
	const CoordinateType coordType = HM) {

	switch (coordType) {
		
		case MV : { // mean value coordinates
			MeanValueR2 mvc(poly);
			mvc.bc(tp);
		} break;

		case HM : { // harmonic coordinates
			HarmonicR2 hmc(poly);
			hmc.setMesh(tp, tf);
			hmc.bc(tp);
		} break;

		case LC : { // local coordinates
			LocalR2 lcc(poly);
			lcc.setEdgeLength(edgeLength);
			lcc.setMesh(tp, tf);
			lcc.bc(tp);
		} break;

		default : { // default - harmonic coordinates
			HarmonicR2 hmc(poly);
			hmc.setMesh(tp, tf);
			hmc.bc(tp);
		} break;
	}

	setDefaultFlags(tp);
}

// Main function.
int main() {

    // Create polygon.
    std::vector<VertexR2> poly(8);

    poly[0] = VertexR2(0.0, 0.0);
    poly[1] = VertexR2(1.0, 0.0);
    poly[2] = VertexR2(2.0, 0.0);
    poly[3] = VertexR2(2.0, 1.0);
    poly[4] = VertexR2(1.0, 1.0);
    poly[5] = VertexR2(1.0, 2.0);
    poly[6] = VertexR2(0.0, 2.0);
    poly[7] = VertexR2(0.0, 1.0);

    const size_t numSubSteps = 3; // number of subdivision steps


    // Example 1: initialize mesh with a polygon and compute 
    // Loop coordinates based on the triangulation of the polygon
    // involving no interior vertices. We also perform one ternary
    // subdivision step as a preprocessing.
    bool makeTernary = true;
    TriangleSubdivisionR2 tsub;

    tsub.setFromPolygon(poly);
    tsub.preprocess(makeTernary);
    tsub.subdivide(numSubSteps);

    std::cout << "\nExample 1/3: DONE!" << std::endl;

    // If you need explicit faces of the mesh, do not
    // forget to create them. You may exclude creation of the face
    // neighbours by specifying the corresponding option below.
    // const bool makeFaceNeighbours = false;
    // tsub.createFaces(makeFaceNeighbours);

    // Another hint: the midpoint subdivision scheme for
    // triangle meshes can be accessed as
    // const bool midpoint = true;
    // tsub.subdivide(numSubSteps, midpoint);


    // Example 2: start from a coarse triangle mesh with some initial  
    // coordinates computed on top of this mesh and then subdivide it.

    // Create triangle mesh.
    const double edgeLength = 0.3;

    std::vector<VertexR2> tp;
    std::vector<Face> tf;

    createTriangleMesh(edgeLength, poly, tf, tp);

    // Compute initial coordinates: 
    // MV - mean value
    // HM - harmonic - default
    // LC - local
    computeInitialCoordinates(edgeLength, poly, tf, tp);

    // Initialize mesh, subdivide this mesh, and obtain
    // the corresponding barycentric coordinates.
    // Since initial mesh is dense enough, we use no ternary step here.
    makeTernary = false;

    tsub.initialize(tp, tf);
    tsub.preprocess(makeTernary);
    tsub.subdivide(numSubSteps);
    
    std::cout << "\nExample 2/3: DONE!\n" << std::endl;


    // Example 3: load from the file a quad mesh with no interior vertices and the initial
    // coordinates set according to the Lagrange property and subdivide it.

    QuadSubdivisionR2 qsub;

    std::cout << "WARNING: Change path to the folder with data in the file main.cpp!" << std::endl;

    qsub.loadMesh("path_to_the_folder/sbc/2d/data/cc/quad_test.obj");
    qsub.loadBarycentricCoordinates("path_to_the_folder/sbc/2d/data/cc/quad_test.data");

    qsub.subdivide(numSubSteps);

    std::cout << "Example 3/3: DONE!" << std::endl;
}
