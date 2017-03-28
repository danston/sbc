// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// README:
/*

    Implementation of the modified Loop triangle-based subdivision in R2
    for computing barycentric coordinates.

    This class depends on:
    1. Face.hpp
    2. VertexR2.hpp
    3. MeshR2.hpp
    4. TriangleCoordinatesR2.hpp

*/

#ifndef GBC_TRIANGLESUBDIVISIONR2_HPP
#define GBC_TRIANGLESUBDIVISIONR2_HPP

// STL includes.
#include <iostream>
#include <cassert>
#include <limits>
#include <cmath>

// Local includes.
#include "Face.hpp"
#include "VertexR2.hpp"
#include "MeshR2.hpp"
#include "TriangleCoordinatesR2.hpp"

namespace gbc {

	// Triangle-based subdivision in R2.
	// More details about this concrete implementation can be found in the paper:
    // D. Anisimov, C. Deng, and K. Hormann. Subdividing barycentric coordinates.
    // Computer Aided Geometric Design, 43:172-185, 2016.
    // Please also note that this implementation does not contain all the
    // features described in the paper!
	class TriangleSubdivisionR2 : public MeshR2 {

	public:
    	// Constructor.
		TriangleSubdivisionR2() { }

		// Load triangle mesh.
		inline void loadMesh(const std::string &path) {
			load(path);
		}

		// Initialize mesh with a polygon and triangulate it.
		void setFromPolygon(const std::vector<VertexR2> &poly) {

			clear();
			setTriangulation(poly);
			setMeshFlags(_f);
            initializeHalfedges(_f);
            buildHEConnectivity();
			setInitialBarycentricCoordinates();
		}

		// Subdivision.

		// Vertex adjustment near the polygon's corners.
		// This version of the function preprocess() is not optimized 
		// and may have slow performance in some cases!
        void preprocess(const bool makeTernary) {

        	assert(isTriangleMesh());

	        // Do a ternary step.
	        if (makeTernary) ternaryStep();
	        
	        // Save mesh vertices and faces before the corner modification
	        // in order to sample linearly barycentric coordinates later.
	        std::vector<VertexR2> tmpv = _v;
	        std::vector<Face> tmpf = _f;
	        
	        // Modify 1-ring neighbourhood around each concave corner.
	        modifyCorners(tmpv, tmpf);
        }

        // Do one subdivision step with the midpoint or modified Loop subdivision scheme.
        // The default scheme is Loop.
        void subdivideMesh(const bool midpoint = false) {

			assert(isTriangleMesh());

	        ///////////////////////////////////
	        // Do ODD step - add new points. //
	        ///////////////////////////////////

	        // Initiate some data.
	        const size_t numV  = numVertices();
	        const size_t numE  = numEdges();
	        const size_t numHE = numHalfedges();
	        
	        _v.resize(numV + numE);
	        std::vector<int> neighs(4);
	        
	        size_t idxV = numV;
	        
	        // Loop over all halfedges.
	        for (size_t i = 0; i < numHE; ++i) {
	            
	            // Consider only the "big brother" of each halfedge-pair.
	            const int i_ = _he[i].neigh;
	            
	            if ((int) i > i_) {    
	                // Local configuration:
	                /*
	                         w1
	                        /||\
	                      /  ||  \
	                    /   ^||    \
	                  w2     []     w4
	                    \   *||    /
	                      \  ||  /
	                        \||/
	                         w3
	                */
	                // The current halfedge and its direction are indicated by "*" and "^".
	                // Create an odd vertex at the position "[]".
	                // NOTE: the right triangle and w4 do not exist if "*" is at the boundary.
	                
	                neighs[0] = _he[i].dest;
	                neighs[1] = _he[_he[i].next].dest;
	                neighs[2] = _he[_he[i].prev].dest;
	                
	                // Mask for odd interior vertices.
	                if (i_ >= 0) {
	                    neighs[3] = _he[_he[i_].next].dest;

	                    if (!midpoint) createOddIntVertex(neighs, _v[idxV]);
	                    else createOddIntVertexMid(neighs, _v[idxV]);
	                }
	                // Mask for odd boundary vertices.
	                else createOddBndVertex(neighs, _v[idxV]);

	                idxV++; 
	            } // if a halfedge is the "big brother"
	            
	        } // for all halfedges
	        
	        // Update halfedge connectivity.
	        updateHEConnectivity(numV);
	        
	        /////////////////////////////////////
	        // Do EVEN step - move old points. //
	        /////////////////////////////////////
	        
	        for (size_t i = 0; i < numV; ++i) {
	            getRing(i, neighs);
	            
	            if (!midpoint) {

	                // Mask for even interior vertices.
	                if (_v[i].type == INTERIOR) moveEvenIntVertex(neighs, _v[i]);
	                // Mask for even boundary vertices.
	                else moveEvenBndVertex(neighs, _v[i]);
	            }
	        }
        }

    private:
    	// Triangulate the polygon.
    	void setTriangulation(const std::vector<VertexR2> &v) {

            TriangulatorR2 tri(v);

            tri.setPlanarGraph(true);
            tri.allowSteinerPoints(false);

            tri.generate(_v, _f);

            assert(_v.size() != 0 && _f.size() != 0);
    	}

	    // Update halfedge connectivity in the mesh after each subdivision step.
	    void updateHEConnectivity(const size_t numV) {

	        const size_t numHE = numHalfedges();
	        const size_t numF  = numHE / 3;
	        
	        _he.resize(4 * numHE);
	        
	        // Update an outgoing halfedge of each even vertex.
	        for (size_t i = 0; i < numV; ++i) {
	            const int k = _v[i].out;
	            
	            // An outgoing halfedge k is the local edge e: k % 3 in the triangle t: k / 3.
	            // The new halfedges of that triangle have indices starting at...
	            const int base = 3 * ((int) numF + 3 * (k / 3));
	            
	            // In general, an offset of the new halfedge e: (i1,j1) is 3 * j1 + i1.
	            // Here, we want e: (j,j+1) with j = k % 3.
	            _v[i].out = base + 3 * ((k + 1) % 3) + k % 3;
	        }
	        
	        // Fix the connectivity between halfedges that replace the "old" edges,
	        // set their destination, and set an outgoing halfedge for each odd vertex.
	        
	        size_t idxV = numV;
	        
	        // Loop over all halfedges.
	        for (size_t i = 0; i < numHE; ++i) {
	            
	            // Consider only the "big brother" of each halfedge-pair.
	            const int i_ = _he[i].neigh;
	            if ((int) i > i_) {
	                
	                // Get indices of halfedges that will replace the current halfedge.
	                // Those are e: (k,k+1) and e: (k,k+2) with k = i % 3 (as above).
	                int base = 3 * ((int) numF + 3 * (i / 3));
	                const int i0 = base + 3 * ((i + 1) % 3) + i % 3;
	                const int i1 = base + 3 * ((i + 2) % 3) + i % 3;
	                
	                // Set their destination.
	                _he[i0].dest = (int) idxV;
	                _he[i1].dest = _he[i].dest;
	                
	                // Set an outgoing halfedge for odd vertex on the current halfedge.
	                _v[idxV].out = i1;
	                
	                // Is it an interior halfedge?
	                if (i_ >= 0) {

	                    // Get indices of halfedges that will replace the neighbouring halfedge.
	                    base = 3 * ((int) numF + 3 * (i_ / 3));
	                    const int i0_ = base + 3 * ((i_ + 1) % 3) + i_ % 3;
	                    const int i1_ = base + 3 * ((i_ + 2) % 3) + i_ % 3;
	                    
	                    // Mutually connect new halfedges.
	                    _he[i0].neigh = i1_;
	                    _he[i1].neigh = i0_;
	                    
	                    _he[i0_].neigh = i1;
	                    _he[i1_].neigh = i0;
	                    
	                    // Also set destination of new neighbouring halfedges.
	                    _he[i0_].dest = (int) idxV;
	                    _he[i1_].dest = _he[i_].dest;
	                }
	                else {
	                    // If the current halfedge is at the boundary, then their children are, too.
	                    _he[i0].neigh = -1;
	                    _he[i1].neigh = -1;
	                }
	                idxV++;
	                
	            } // if a halfedge is the "big brother"
	            
	        } // for all halfedges
	        
	        // Fix connectivity and destination of new halfedges "inside" old triangles.
	        size_t Ej = 0;
	        size_t Ejj[3], Ejmj[3];
	        for (size_t i = 0; i < numF; ++i) {
	            
	            // Explicitely set indices of new halfedges e: (0,1), e: (1,2), e: (2,0).
	            const size_t base = 3 * (numF + 3 * i);
	            Ejmj[1] = base + 3;
	            Ejmj[2] = base + 7;
	            Ejmj[0] = base + 2;
	            
	            // And e: (0,0), e: (1,1), e: (2,2).
	            Ejj[0] = base;
	            Ejj[1] = base + 4;
	            Ejj[2] = base + 8;
	            for (size_t j = 0; j < 3; ++j, Ej++) {

	                // Connect new halfedges with each other.
	                _he[Ej].neigh     = (int) Ejj[j];
	                _he[Ejj[j]].neigh = (int) Ej;
	                
	                // Set their destination.
	                idxV = (size_t) _he[Ejmj[j]].dest;

	                _he[Ej].dest = (int) idxV;
	                _he[Ejj[(j + 1) % 3]].dest = (int) idxV;
	            }
	        }
	        
	        // Set next and previous halfedges for all new halfedges.
	        size_t idxE = numHE;
	        for (size_t i = numF; i < 4 * numF; ++i) {

	            const size_t base = idxE;
	            for (size_t j = 0; j < 3; ++j, idxE++) {

	                _he[idxE].prev = (int) base + (j + 2) % 3;
	                _he[idxE].next = (int) base + (j + 1) % 3;
	            }
	        }
	    }

	    // Create an interior odd vertex - midpoint scheme.
	    void createOddIntVertexMid(const std::vector<int> &neighs, VertexR2 &newV) const {
	        
	        assert(neighs.size() == 4);

	        // Create a new interior vertex.
	        newV  = _v[neighs[0]];
	        newV *= 0.5;
	        newV += 0.0 * _v[neighs[1]];
	        newV += 0.5 * _v[neighs[2]];
	        newV += 0.0 * _v[neighs[3]];
	        
	        // Update its flags.
	        newV.val   = 6;
	        newV.out   = -1;
	        newV.alpha = -1.0;
	        newV.type  = INTERIOR;
	    }

	    // Create an interior odd vertex - modified Loop scheme.
	    void createOddIntVertex(const std::vector<int> &neighs, VertexR2 &newV) const {

	        assert(neighs.size() == 4);
	        
	        double c[4] = { 0.375, 0.125, 0.375, 0.125 };
	        
	        // Zorin modification.
	        if (_v[neighs[0]].type == CONVEX || _v[neighs[0]].type == CONCAVE) {
	            
	            const int k        = _v[neighs[0]].val - 1;
	            const double alpha = _v[neighs[0]].alpha;
	            const double theta = (2.0 * M_PI - alpha) / k;
	            double gamma       = 0.5 - 0.25 * cos(theta);
	            
	            // If both vertices of an edge are convex or concave corners,
	            // we return to the original weights 3/8.
	            if(_v[neighs[2]].type == CONVEX || _v[neighs[2]].type == CONCAVE) gamma = 0.375;
	            
	            // Both weights must behave as 1/2.
	            c[0] = 0.75 - gamma;
	            c[2] = gamma;
	        }
	        
	        // Zorin modification.
	        if (_v[neighs[2]].type == CONVEX || _v[neighs[2]].type == CONCAVE) {
	            
	            const int k        = _v[neighs[2]].val - 1;
	            const double alpha = _v[neighs[2]].alpha;
	            const double theta = (2.0 * M_PI - alpha) / k;
	            double gamma       = 0.5 - 0.25 * cos(theta);
	            
	            // If both vertices of an edge are convex or concave corners,
	            // we return to the original weights 3/8.
	            if(_v[neighs[0]].type == CONVEX || _v[neighs[0]].type == CONCAVE) gamma = 0.375;
	            
	            // Both weights must behave as 1/2.
	            c[2] = 0.75 - gamma;
	            c[0] = gamma;
	        }
	        
	        // Create a new interior vertex.
	        newV  = _v[neighs[0]];
	        newV *= c[0];
	        newV += c[1] * _v[neighs[1]];
	        newV += c[2] * _v[neighs[2]];
	        newV += c[3] * _v[neighs[3]];
	        
	        // Update its flags.
	        newV.val   = 6;
	        newV.out   = -1;
	        newV.alpha = -1.0;
	        newV.type  = INTERIOR;
	    }

	    // Create a boundary odd vertex.
	    void createOddBndVertex(const std::vector<int> &neighs, VertexR2 &newV) const {
	        
	        assert(neighs.size() > 0);
	        
	        // Create a new boundary vertex.
	        newV  = _v[neighs[0]];
	        newV += _v[neighs[2]];
	        newV *= 0.5;
	        
	        // Update its flags.
	        newV.val   = 4;
	        newV.out   = -1;
	        newV.alpha = -1.0;
	        newV.type  = FLAT;
	    }

	    // Move an interior even vertex.
	    void moveEvenIntVertex(const std::vector<int> &neighs, VertexR2 &centre) const {

	        const size_t n = neighs.size();
	        
	        assert(n > 0);
	        assert(centre.alpha == -1.0);
	        assert(centre.type == INTERIOR);
	        
	        // Compute weights for arbitrary valency.
	        double c0  = 0.375 + 0.25 * cos(2.0 * M_PI / n);
	        c0 *= c0 * 1.6;
	        
	        // Compute the barycentre of a set of vertices.
	        VertexR2 barycentre = _v[neighs[0]];
	        for (size_t i = 1; i < n; ++i) barycentre += _v[neighs[i]];
	        
	        const double nInv = 1.0 / n;
	        
	        // Move an old interior vertex.
	        centre *= c0;
	        centre += (nInv * (1.0 - c0)) * barycentre;
	    }

	    // Move a boundary even vertex.
	    void moveEvenBndVertex(const std::vector<int> &neighs, VertexR2 &centre) const {

	        assert(centre.type != INTERIOR);
	        
	        if (centre.type == FLAT) { // FLAT vertex
	            
	            assert(centre.alpha == -1.0);
	            assert(neighs.size() > 1);
	            
	            const int n = centre.val;
	            
	            assert(n > 1);
	            assert((int) neighs.size() == n);
	            
	            centre *= 2.0;
	            centre += _v[neighs[0]];
	            centre += _v[neighs[n - 1]];
	            centre *= 0.25;
	        }
	        
	        // CONVEX AND CONCAVE vertices are interpolated!
	    }

    	// Ternary subdivision step.
	    void ternaryStep() {

	        const size_t numF  = numFaces();
	        const size_t numV  = numVertices();
	        const size_t numHE = numHalfedges();
	        const size_t numE  = numEdges();
	        
	        assert(numF != 0 && numV != 0 && numHE != 0 && numE != 0);
	        
	        std::vector<VertexR2> tmpv(numV + numF + 2 * numE);
	        std::vector<Face> tmpf(9 * numF);
	        
	        // Create temporary vertices.
	        
	        // Old vertices.
	        for (size_t i = 0; i < numV; ++i) {
	            
	            tmpv[i].x() = _v[i].x();
	            tmpv[i].y() = _v[i].y();
	            tmpv[i].b() = _v[i].b();
	        }
	        
	        // New face vertices.
	        for (size_t i = 0; i < numF; ++i) 
	        	createFaceBarycenter(_f[i].v[0], _f[i].v[1], _f[i].v[2], tmpv[numV + i]);
	        
	        // New edge vertices.
	        size_t count = 0;
	        for (size_t i = 0; i < numHE; ++i) {
	            const int i_ = _he[i].neigh;

	            if ((int) i > i_) {
	                
	                const int i1 = _he[_he[i].prev].dest;
	                const int i2 = _he[i].dest;
	                
	                assert(i1 != -1 && i2 != -1);
	                
	                trisectEdge(i1, i2, tmpv[numV + numF + count], tmpv[numV + numF + count + 1]);
	                count += 2;
	            }
	        }

	        // Create temporary faces.

	        VertexR2 tmpV1, tmpV2;
	        std::vector<int> neighs(6);
	        
	        // Create new faces.
	        count = 0;
	        for (size_t i = 0; i < numF; ++i) {
	            
	            const int i0 = _f[i].v[0];
	            const int i1 = _f[i].v[1];
	            const int i2 = _f[i].v[2];
	            
	            assert(i0 != -1 && i1 != -1 && i2 != -1);
	            
	            trisectEdge(i0, i1, tmpV1, tmpV2);
	            neighs[0] = findIndex(tmpv, tmpV1);
	            neighs[1] = findIndex(tmpv, tmpV2);

	            assert(neighs[0] != -1 && neighs[1] != -1);
	            
	            trisectEdge(i1, i2, tmpV1, tmpV2);
	            neighs[2] = findIndex(tmpv, tmpV1);
	            neighs[3] = findIndex(tmpv, tmpV2);

	            assert(neighs[2] != -1 && neighs[3] != -1);
	            
	            trisectEdge(i2, i0, tmpV1, tmpV2);
	            neighs[4] = findIndex(tmpv, tmpV1);
	            neighs[5] = findIndex(tmpv, tmpV2);

	            assert(neighs[4] != -1 && neighs[5] != -1);
	            
	            for (size_t j = 0; j < 6; ++j) {
	                
	                tmpf[count].v[0] = numV + i;
	                tmpf[count].v[1] = neighs[j];

	                tmpf[count++].v[2] = neighs[(j + 1) % 6];
	            }

	            tmpf[count].v[0] = _f[i].v[0]; tmpf[count].v[1] = neighs[0]; tmpf[count++].v[2] = neighs[5];
	            tmpf[count].v[0] = _f[i].v[1]; tmpf[count].v[1] = neighs[2]; tmpf[count++].v[2] = neighs[1];
	            tmpf[count].v[0] = _f[i].v[2]; tmpf[count].v[1] = neighs[4]; tmpf[count++].v[2] = neighs[3];
	        }
	        
	        // Clear old and initialize new mesh.

	        clear();
	        initialize(tmpv, tmpf);

	        // Copy back some flags.
	        for (size_t i = 0; i < numVertices(); ++i) {
	            if (i > numV - 1 && _v[i].type != INTERIOR) {

	                _v[i].type  = FLAT;
	                _v[i].alpha = -1.0;
	            }
	        }
	    }

	    // Modify corners of the mesh.
	    void modifyCorners(std::vector<VertexR2> &tmpv, std::vector<Face> &tmpf) {

	        const size_t numV = numVertices();
	        std::vector<int> neighs;
	        
	        assert(numV != 0);
	        
	        // Go through all the mesh vertices.
	        for (size_t i = 0; i < numV; ++i) {
	            if (_v[i].type == CONCAVE) { // Concave vertex
	                
	                // Find 1-ring neighbourhood of the corner vertex.
	                getRing(i, neighs);
	                const size_t nSize = neighs.size();
	                
	                for (size_t j = 1; j < nSize - 1; ++j) 
	                    if (_v[neighs[j]].type != INTERIOR)
	                        std::cerr << "\nERROR: Adjusted vertex is on the boundary!\n" << std::endl;

	                // Find the internal angle for the corner and split it into equal parts.
	                const double intAngle = 2.0 * M_PI - _v[i].alpha;
	                const double angle    = intAngle / (nSize - 1);
	                
	                // Find the closest neighbour to the corner.
	                double minLen = std::numeric_limits<double>::max();
	                for (size_t j = 0; j < nSize; ++j) {

	                    const double len = (_v[i] - _v[neighs[j]]).length();
	                    minLen = std::min(minLen, len);
	                }
	                
	                // Rotate the first neighbour of the corner vertex by the angle 'angle'
	                // to get equally-distant and distributed mesh vertices around this corner.
	                double tmpAngle = 0.0;
	                VertexR2 &firstV = _v[neighs[0]];
	                
	                VertexR2 tmpV;
	                for (size_t j = 0; j < nSize; ++j) {
	                    tmpAngle = j * angle;
	                    
	                    tmpV = firstV.translated(-_v[i]);
	                    tmpV = tmpV.rotated(tmpAngle);
	                    
	                    const double len = tmpV.length();
	                    const double scaleFactor = minLen / len;
	                    
	                    tmpV = tmpV.scaled(scaleFactor, scaleFactor);
	                    tmpV = tmpV.translated(_v[i]);

	                    _v[neighs[j]].x() = tmpV.x();
	                    _v[neighs[j]].y() = tmpV.y();
	                }
	                
	                // Scale mesh vertices around the corner to avoid interior foldovers.
	                scaleVertices(i, neighs);
	                
	                // Fix barycentric coordinates for all the neighbours of the corner vertex.
	                fixCoordinates(neighs, tmpv, tmpf);
	            }
	        }
	    }

	    // Create the barycenter point in a triangle.
	    void createFaceBarycenter(const int i0, const int i1, const int i2, VertexR2 &barycentre) const {
	        
	        // Find the triangle barycentre.
	        barycentre  = _v[i0];
	        barycentre += _v[i1];
	        barycentre += _v[i2];
	        
	        barycentre *= (1.0 / 3.0);
	        
	        // Update its flags.
	        barycentre.val   = 0;
	        barycentre.out   = -1;
	        barycentre.alpha = -1.0;
	        barycentre.type  = INTERIOR;
	    }

	    // Trisect an edge.
	    void trisectEdge(const int i1, const int i2, VertexR2 &newV1, VertexR2 &newV2) const {

	        // Trisect the edge.
	        const double third = 1.0 / 3.0;
	        
	        newV1  = _v[i1];
	        newV1 *= (1.0 - third);
	        newV1 += third * _v[i2];
	        
	        newV2  = _v[i1];
	        newV2 *= third;
	        newV2 += (1.0 - third) * _v[i2];
	        
	        // Update flags.
	        newV1.val   = 0;    	newV2.val  = 0;
	        newV1.out   = -1;   	newV2.out  = -1;
	        newV1.alpha = -1.0;     newV2.alpha = -1.0;
	        newV1.type  = INTERIOR; newV2.type = INTERIOR;
	    }

	    // Find vertex index.
	    int findIndex(const std::vector<VertexR2> &tmpv, const VertexR2 &query) const {

	        const size_t numV = tmpv.size();
	        
	        for (size_t i = 0; i < numV; ++i)
	            if (tmpv[i] == query) return i;
	        
	        // Pointer cannot be here.
	        assert(false);
	        return -1;
	    }

	    // Scale mesh vertices around the corner.
	    void scaleVertices(const size_t centreInd, const std::vector<int> &neighs) {

	        const size_t nSize = neighs.size();
	        
	        size_t count = 0;
	        const size_t maxCount = 25;
	        
	        // Until all the interior foldovers are gone, scale neighbours down by 1/2.
	        while (foldover(neighs)) {

	            for (size_t i = 0; i < nSize; ++i) {
	                _v[neighs[i]] += _v[centreInd];
	                _v[neighs[i]] *= 0.5;
	            }
	            ++count;
	            if (count > maxCount)
	                std::cerr << "ERROR: Max number of iterations is reached! Function: scaleVertices() failed!\n" << std::endl;
	        }
	    }

	    // Check if no mesh vertex around the corner creates a foldover.
	    bool foldover(const std::vector<int> &firstRing) {

	        const size_t nSizeFirst = firstRing.size();
	        assert(nSizeFirst != 0);
	        
	        // Go through all the neighbours of the corner vertex
	        // and find 1-ring neighbourhood for each of them.
	        std::vector<int> secondRing;
	        for (size_t i = 0; i < nSizeFirst; ++i) {

	            getRing(firstRing[i], secondRing);
	            const size_t nSizeSecond = secondRing.size();

	            size_t size = nSizeSecond;
	            assert(nSizeSecond != 0);
	            
	            if (i == 0 || i == nSizeFirst - 1) size--; // for boundary vertices the size is 1 element smaller
	            const VertexR2 &v2 = _v[firstRing[i]];
	            
	            assert(size != 0);
	            
	            // Find determinant for each triangle and check if it is negative/positive.
	            for (size_t j = 0; j < size; ++j) {

	                const VertexR2 &v3 = _v[secondRing[j]];
	                const VertexR2 &v1 = _v[secondRing[(j + 1) % nSizeSecond]];
	                
	                const double x1 = v1.x(); const double y1 = v1.y();
	                const double x2 = v2.x(); const double y2 = v2.y();
	                const double x3 = v3.x(); const double y3 = v3.y();
	                
	                // Compute the following determinant:
	                //
	                //       | x1 y1 1 |
	                // det = | x2 y2 1 |
	                //       | x3 y3 1 |
	                //
	                // if det > 0 or det == 0, no interior foldover is detected;
	                // if det < 0, an interior foldover happend;
	                
	                const double det = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
	                
	                if (det < 0.0) return true; // foldover found
	            }
	        }
	        return false; // no foldover is detected
	    }

	    // Fix barycentric coordinates around the corner.
	    void fixCoordinates(const std::vector<int> &neighs, const std::vector<VertexR2> &tmpv, const std::vector<Face> &tmpf) {

	        const size_t nSize = neighs.size();
	        const size_t numF  = tmpf.size();
	        const size_t numC  = _v[0].b().size();
	        
	        // Go through all the neighbours.
	        for (size_t i = 0; i < nSize; ++i) {
	            
	            // Go through all the faces of the mesh before a ternary step.
	            size_t count = 0;
	            for (size_t j = 0; j < numF; ++j) {
	                
	                // Compute triangle coordinates.
	                std::vector<double> lambda;

	                TriangleCoordinatesR2 tc(tmpv[tmpf[j].v[0]], tmpv[tmpf[j].v[1]], tmpv[tmpf[j].v[2]]);
	                tc.compute(_v[neighs[i]], lambda);
	                
	                // If we found a face to which current neighbour belongs,
	                // linearly sample barycentric cordinates from the vertices of this face.
	                if (lambda[0] >= 0.0 && lambda[1] >= 0.0 && lambda[2] >= 0.0) {
	                    for (size_t k = 0; k < numC; ++k) {

	                        _v[neighs[i]].b()[k] = 
	                        std::fabs(lambda[0]) * tmpv[tmpf[j].v[0]].b()[k] +
	                        std::fabs(lambda[1]) * tmpv[tmpf[j].v[1]].b()[k] +
	                        std::fabs(lambda[2]) * tmpv[tmpf[j].v[2]].b()[k];

	                    } break;
	                } else ++count;
	            }
	            
	            // If we did not find at least one face to which current neighbour belongs, return an error.
	            if (count == numF && numC != 0)
	                std::cerr << "ERROR: No appropriate face was found! fixCoordinates() function failed!\n" << std::endl;
	        }
	    }
	};

} // namespace gbc

#endif // GBC_TRIANGLESUBDIVISIONR2_HPP
