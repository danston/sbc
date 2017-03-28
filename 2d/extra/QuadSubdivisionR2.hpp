// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// README:
/*

    Implementation of the modified Catmull-Clark quad-based subdivision in R2
    for computing barycentric coordinates.

    This class depends on:
    1. Face.hpp
    2. VertexR2.hpp
    3. MeshR2.hpp

*/

#ifndef GBC_QUADSUBDIVISIONR2_HPP
#define GBC_QUADSUBDIVISIONR2_HPP

// STL includes.
#include <iostream>
#include <cassert>
#include <limits>
#include <cmath>

// Local includes.
#include "Face.hpp"
#include "VertexR2.hpp"
#include "MeshR2.hpp"

namespace gbc {

	// Quad-based subdivision in R2.
	// More details about this concrete implementation can be found in the paper:
    // D. Anisimov, C. Deng, and K. Hormann. Subdividing barycentric coordinates.
    // Computer Aided Geometric Design, 43:172-185, 2016.
    // Please also note that this implementation does not contain all the
    // features described in the paper!
	class QuadSubdivisionR2 : public MeshR2 {

	public:
    	// Constructor.
		QuadSubdivisionR2() { }

		// Load quad mesh.
		inline void loadMesh(const std::string &path) {
			load(path);
		}

		// Initialize mesh with a polygon and quadrangulate it.
		void setFromPolygon(const std::vector<VertexR2> &poly) {

			clear();
			setQuadrangulation(poly);
			setMeshFlags(_f);
            initializeHalfedges(_f);
            buildHEConnectivity();
			setInitialBarycentricCoordinates();
		}

		// Subdivision.

		// Vertex adjustment near the polygon's corners.
		// This version of the function preprocess() is not optimized 
		// and may have slow performance in some cases!
        void preprocess(const bool) {

        	assert(isQuadMesh());

	        // Save mesh vertices and faces before the corner modification
	        // in order to sample linearly barycentric coordinates later.
	        std::vector<VertexR2> tmpv = _v;
	        std::vector<Face> tmpf = _f;
	        
	        // Modify 1-ring neighbourhood around each concave corner.
	        modifyCorners(tmpv, tmpf);
        }

        // Do one subdivision step with the modified Catmull-Clark subdivision scheme.
        void subdivideMesh(const bool) {

			assert(isQuadMesh());

	        ///////////////////////////////////
	        // Do ODD step - add new points. //
	        ///////////////////////////////////
	        
	        // Initiate some data.
	        const size_t numV  = numVertices();
	        const size_t numE  = numEdges();
	        const size_t numHE = numHalfedges();
	        const size_t numF  = numHE / 4;
	        
	        _v.resize(numV + numE + numF);
	        std::vector<int> neighs(6);
	        
	        size_t idxV = numV;
	        
	        // Loop over all halfedges.
	        for (size_t k = 0; k < numHE; ++k) {
	            
	            // Consider only the "big brother" of each halfedge-pair.
	            const int k_ = _he[k].neigh;
	            if ((int) k > k_) {
	                
	                // Local configuration:
	                /*
	                 w2---w1---w5
	                 ||   ||   ||
	                 ||   ||   ||
	                 ||  ^||   ||
	                 ||   []   ||
	                 ||  *||   ||
	                 ||   ||   ||
	                 ||   ||   ||
	                 w3---w0---w4
	                 */
	                // The current halfedge and its direction are indicated by "*" and "^".
	                // Create an odd vertex at the position "[]".
	                // NOTE: the right quad, w4, and w5 do not exist if "*" is at the boundary.
	                
	                neighs[0] = _he[_he[k].prev].dest;
	                neighs[1] = _he[k].dest;
	                
	                neighs[2] = _he[_he[k].next].dest;
	                neighs[3] = _he[_he[_he[k].next].next].dest;
	                
	                // Mask for interior edge vertices.
	                if (k_ >= 0) {

	                    neighs[4] = _he[_he[k_].next].dest;
	                    neighs[5] = _he[_he[_he[k_].prev].prev].dest;
	                    createIntEdgeVertex(neighs, _v[idxV]);
	                }
	                // Mask for boundary edge vertices.
	                else createBndEdgeVertex(neighs, _v[idxV]);
	                idxV++;
	                
	            } // if a halfedge is the "big brother"
	            
	        } // for all halfedges
	        
	        // Add new face vertices.
	        int base = 0, c = 0, idxE = 0;
	        for (size_t k = 0; k < numHE; ++k) {
	            if (k % 4 == 0) {

	                base = 0;
	                if (_he[0].next == 1) idxE = _he[k].prev;
	                else idxE = _he[c].prev;
	            }
	            
	            neighs[base++] = _he[idxE].dest;
	            idxE = _he[idxE].next;
	            
	            if (base == 4) {
	                createFaceVertex(neighs, _v[idxV++]);
	                c++;
	            }
	        }
	        
	        // Update halfedge connectivity.
	        updateHEConnectivity(numV);
	        
	        /////////////////////////////////////
	        // Do EVEN step - move old points. //
	        /////////////////////////////////////
	        
	        for (size_t i = 0; i < numV; ++i) {
	            
	            // Mask for even interior vertices.
	            if (_v[i].type == INTERIOR) {
	                
	                size_t nSize = 0;
	                int first, curr, prev, next;
	                
	                // Find neighbours of the vertex.
	                curr = first = _v[i].out;
	                
	                assert(_v[i].val > 0);
	                neighs.resize(2 * (size_t) _v[i].val);
	                do {
	                    neighs[nSize++] = _he[curr].dest;
	                    
	                    next = _he[curr].next;
	                    
	                    neighs[nSize++] = _he[next].dest;
	                    
	                    prev = _he[curr].prev;
	                    curr = _he[prev].neigh;
	                    
	                } while((curr >= 0) && (curr != first));
	                
	                // Add one more neighbour for boundary vertices.
	                if(curr < 0) {
	                    curr = _he[prev].prev;
	                    neighs[nSize] = _he[curr].dest;
	                }
	                
	                moveEvenIntVertex(neighs, _v[i]);
	            }
	            // Mask for even boundary vertices.
	            else {
	                getRing(i, neighs);
	                moveEvenBndVertex(neighs, _v[i]);
	            }
	        }
		}

	private:
    	// Quadrangulate the polygon.
    	void setQuadrangulation(const std::vector<VertexR2>&) {
    		// Not yet implemented.
    	}

		// Update halfedge connectivity in the mesh after each subdivision step.
	    void updateHEConnectivity(const size_t numV) {

	        const size_t numHE = numHalfedges();
	        
	        _he.resize(4 * numHE);
	        
	        // Update connectivity for all new edges created on top of the old edges.
	        size_t idxV = numV, count = 0;
	        for (size_t k = 0; k < numHE; ++k) {

	            const int k_ = _he[k].neigh;
	            
	            // Consider only the "big brother" of each halfedge-pair.
	            if ((int) k > k_) {
	                if (k_ >= 0) { // interior edge
	                    
	                    int oldV = _he[k].dest;
	                    const int oldE = _he[k].neigh;
	                    
	                    _he[k].dest = (int) idxV;
	                    _he[count + numHE].dest = oldV;
	                    
	                    oldV = _he[oldE].dest;
	                    
	                    _he[oldE].dest = (int) idxV;
	                    _he[count + numHE + 1].dest = oldV;
	                    
	                    _he[count + numHE + 1].neigh = _he[oldE].neigh;
	                    _he[count + numHE].neigh = _he[k].neigh;
	                    
	                    _he[k].neigh = count + numHE + 1;
	                    _he[oldE].neigh = count + numHE;
	                    
	                    int next = _he[k].next;
	                    
	                    _he[count + numHE].next = next;
	                    _he[next].prev = count + numHE;
	                    
	                    next = _he[oldE].next;
	                    
	                    _he[count + numHE + 1].next = next;
	                    _he[next].prev = count + numHE + 1;
	                    
	                    _v[idxV].out = count + numHE;
	                    count++;

	                } else { // boundary edge
	                    
	                    const int oldV = _he[k].dest;
	                    
	                    _he[k].dest = (int) idxV;
	                    _he[k].neigh = -1;
	                    
	                    _he[count + numHE].dest = oldV;
	                    _he[count + numHE].neigh = -1;
	                    
	                    const int next = _he[k].next;
	                    
	                    _he[count + numHE].next = next;
	                    _he[next].prev = count + numHE;
	                    
	                    _v[idxV].out = count + numHE;
	                }
	                idxV++; count++;
	                
	            } // if a halfedge is the "big brother"
	            
	        } // for all halfedges
	        
	        // Update connectivity for all new edges created on top of the old faces.
	        const int next = _he[0].next;
	        
	        int base = 0, c = 0, idxE = 0;
	        for (size_t k = 0; k < numHE; ++k) {
	            if (k % 4 == 0) base = 0; base++;
	            
	            if (base == 4) {
	                if (next == 1) idxE = k - 3;
	                else idxE = c;
	                
	                {
	                    int offset_1 = 0, offset_2 = 7;
	                    
	                    const int init = int(count + numHE);
	                    _v[idxV].out = init + 1;
	                    
	                    for (size_t i = 0; i < 4; ++i) {

	                        const int ind   = init + offset_1;
	                        const int neigh = init + offset_1 + 1;
	                        
	                        _he[ind].neigh = neigh;
	                        _he[neigh].neigh = ind;
	                        
	                        _he[ind].dest = (int) idxV;
	                        _he[neigh].dest = _he[idxE].dest;
	                        
	                        _he[ind].prev = idxE;
	                        _he[ind].next = init + offset_2 % 8;
	                        _he[init + offset_2 % 8].prev = ind;
	                        
	                        const int tmp = _he[idxE].next;
	                        _he[idxE].next = ind;
	                        
	                        _he[neigh].next = _he[tmp].prev;
	                        _he[_he[tmp].prev].prev = neigh;
	                        
	                        idxE = tmp;
	                        offset_1 += 2;
	                        offset_2 += 2;
	                    }
	                }
	                count += 8;
	                idxV++;
	                c++;
	            }
	        }
	    }

	    // Create a face vertex.
	    void createFaceVertex(const std::vector<int> &neighs, VertexR2 &newV) const {
	        
	        assert(neighs.size() == 6);
	        
	        newV = _v[neighs[0]];
	        for (size_t i = 1; i < 4; ++i)
	            newV += _v[neighs[i]];
	        
	        newV *= 0.25;
	        
	        // Update its flags.
	        newV.val   = 4;
	        newV.out   = -1;
	        newV.alpha = -1.0;
	        newV.type  = INTERIOR;
	    }

	    // Create an interior edge vertex.
	    void createIntEdgeVertex(const std::vector<int> &neighs, VertexR2 &newV) const {

	        assert(neighs.size() == 6);
	        
	        double c[6] = { 0.375, 0.375, 0.0625, 0.0625, 0.0625, 0.0625 };
	        
	        // Zorin modification.
	        if (_v[neighs[0]].type == CONVEX || _v[neighs[0]].type == CONCAVE) {
	            
	            const int k        = _v[neighs[0]].val - 1;
	            const double alpha = _v[neighs[0]].alpha;
	            const double theta = (2.0 * M_PI - alpha) / k;
	            double gamma       = 0.375 - 0.25 * cos(theta);
	            
	            // If both vertices of an edge are convex or concave corners,
	            // we return to the original weights 3/8.
	            if (_v[neighs[1]].type == CONVEX || _v[neighs[1]].type == CONCAVE) gamma = 0.375;
	            
	            // Both weights must behave as 1/2.
	            c[0] = 0.75 - gamma;
	            c[1] = gamma;
	        }
	        
	        // Zorin modification.
	        if(_v[neighs[1]].type == CONVEX || _v[neighs[1]].type == CONCAVE) {
	            
	            const int k        = _v[neighs[1]].val - 1;
	            const double alpha = _v[neighs[1]].alpha;
	            const double theta = (2.0 * M_PI - alpha) / k;
	            double gamma       = 0.375 - 0.25 * cos(theta);
	            
	            // If both vertices of an edge are convex or concave corners,
	            // we return to the original weights 3/8.
	            if (_v[neighs[0]].type == CONVEX || _v[neighs[0]].type == CONCAVE) gamma = 0.375;
	            
	            // Both weights must behave as 1/2.
	            c[1] = 0.75 - gamma;
	            c[0] = gamma;
	        }
	        
	        newV  = _v[neighs[0]];
	        newV *= c[0];
	        newV += c[1] * _v[neighs[1]];
	        
	        for (size_t i = 2; i < 6; ++i) newV += c[i] * _v[neighs[i]];
	        
	        // Update its flags.
	        newV.val   = 4;
	        newV.out   = -1;
	        newV.alpha = -1.0;
	        newV.type  = INTERIOR;
	    }

	    // Create a boundary edge vertex.
	    void createBndEdgeVertex(const std::vector<int> &neighs, VertexR2 &newV) const {

	        assert(neighs.size() == 6);
	        
	        newV  = _v[neighs[0]];
	        newV += _v[neighs[1]];
	        newV *= 0.5;
	        
	        // Update its flags.
	        newV.val   = 3;
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
	        
	        // Update weights.
	        const int k = centre.val;
	        const int ks = k * k;
	        
	        const double alpha = - 1.0 / ks;
	        const double beta  =   4.0 / ks;
	        const double gamma = 1.0 - k * alpha - k * beta;
	        
	        centre *= gamma;
	        for (size_t i = 0; i < n; ++i) {
	            
	            if (i % 2 == 0) {
	                centre += beta * _v[neighs[i]];
	            }
	            else {
	                centre += alpha * _v[neighs[i]];
	            }
	        }
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

		// Modify corners of the mesh.
	    void modifyCorners(std::vector<VertexR2> &tmpv, std::vector<Face> &tmpf) {

	        const size_t numV = numVertices();
	        std::vector<int> neighs;
	        
	        assert(numV != 0);
	        
	        // Go through all the mesh vertices.
	        for (size_t i = 0; i < numV; ++i) {
	            if (_v[i].type == CONCAVE) { // Concave vertex
	                
	                // Find 1-ring neighbourhood of the corner vertex.
	                getFullRing(i, neighs);
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
	                    
	                    if (j % 2 == 0) {
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
	                }

	                for (size_t j = 0; j < nSize; ++j) {
	                    if (j % 2 != 0) {

	                        _v[neighs[j]].x() = _v[neighs[j - 1]].x() + _v[neighs[j + 1]].x() - _v[i].x();
	                        _v[neighs[j]].y() = _v[neighs[j - 1]].y() + _v[neighs[j + 1]].y() - _v[i].y();
	                    }
	                }
	                
	                // Scale mesh vertices around the corner to avoid interior fold-overs.
	                // scaleVertices(i, neighs); // to be added later
	                
	                // Fix barycentric coordinates for all the neighbours of the corner vertex.
	                fixCoordinates(neighs, tmpv, tmpf);
	            }
	        }
	    }

	    // Get full ring around a vertex.
	    void getFullRing(const size_t vInd, std::vector<int> &neighs) const {

	        int first, prev, curr, next;
	        
	        // Find neighbours of the vertex.
	        curr = first = _v[vInd].out;
	        
	        assert(_v[vInd].val > 0);

	        int nSize = _v[vInd].val * 2;
	        if (_v[vInd].type != INTERIOR) nSize--;

	        assert(nSize >= 0);

	        neighs.resize(nSize);

	        size_t count = 0;
	        do {
	            neighs[count++] = _he[curr].dest;
	            next = _he[curr].next;
	            neighs[count++] = _he[next].dest; 
	            
	            prev = _he[curr].prev;
	            curr = _he[prev].neigh;

	        } while((curr >= 0) && (curr != first));
	        
	        // Add one more neighbour for boundary vertices.
	        if (curr < 0) {

	            curr = _he[prev].prev;
	            neighs[count] = _he[curr].dest;
	        }
	    }

	    // Fix barycentric coordinates around the corner.
	    void fixCoordinates(const std::vector<int> &neighs, const std::vector<VertexR2> &tmpv, const std::vector<Face> &tmpf) {

	        const size_t nSize = neighs.size();
	        const size_t numF  = tmpf.size();
	        const size_t numC  = _v[0].b().size();
	        
	        std::vector<VertexR2> quad(4);

	        // Go through all the neighbours.
	        for (size_t i = 0; i < nSize; ++i) {
	            
	            // Go through all the faces of the mesh.
	            size_t count = 0;
	            for (size_t j = 0; j < numF; ++j) {
	                
	                // Compute mean value coordinates.
	                std::vector<double> lambda;

	                quad[0] = tmpv[tmpf[j].v[0]];
	                quad[1] = tmpv[tmpf[j].v[1]];
	                quad[2] = tmpv[tmpf[j].v[2]];
	                quad[3] = tmpv[tmpf[j].v[3]];

	                MeanValueR2 mvc(quad);
	                mvc.compute(_v[neighs[i]], lambda);
	                
	                // If we found a face to which current neighbour belongs,
	                // linearly sample barycentric cordinates from the vertices of this face.
	                if (lambda[0] >= 0.0 && lambda[1] >= 0.0 && lambda[2] >= 0.0 && lambda[3] >= 0.0) {
	                    for (size_t k = 0; k < numC; ++k) {

	                        _v[neighs[i]].b()[k] = std::fabs(lambda[0]) * tmpv[tmpf[j].v[0]].b()[k] +
	                                               std::fabs(lambda[1]) * tmpv[tmpf[j].v[1]].b()[k] +
	                                               std::fabs(lambda[2]) * tmpv[tmpf[j].v[2]].b()[k] +
	                                               std::fabs(lambda[3]) * tmpv[tmpf[j].v[3]].b()[k];
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

#endif // GBC_QUADSUBDIVISIONR2_HPP
