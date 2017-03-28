# Subdividing barycentric coordinates Ver. 1.0.0

### Description

This set of C++ classes provides you with an implementation of barycentric coordinates based on Loop and Catmull-Clark subdivision schemes as discussed in the paper: D. Anisimov, C. Deng, and K. Hormann. Subdividing barycentric coordinates. Computer Aided Geometric Design, 43:172-185, 2016. The paper itself can be found [here](http://www.inf.usi.ch/hormann/papers/Anisimov.2016.SBC.pdf).

##### NOTE: This code has been tested only on Mac OS!

### Run the code

In order to run the code, do the following:

1. Install [macports](https://www.macports.org/install.php)
2. Open terminal and type the following:

```bash
  sudo port install cmake
```
```bash
  cd path_to_the_folder/sbc/2d/
```
```bash
  mkdir bin
```
```bash
  cd bin
```
```bash
  cmake -DCMAKE_BUILD_TYPE=Debug ..
```
```bash
  make
```
```bash
  ./sbc
```

For the release version use instead: 

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Example

```C++
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

// Initialize mesh with a polygon and compute 
// Loop coordinates based on the triangulation of the polygon
// involving no interior vertices. We also perform one ternary
// subdivision step as a preprocessing.
bool makeTernary = true;
TriangleSubdivisionR2 tsub;

tsub.setFromPolygon(poly);
tsub.preprocess(makeTernary);
tsub.subdivide(numSubSteps);
```

##### NOTE: For more examples see main.cpp!

### Bugs

If you find any bugs, please report them to me, and I will try to fix them as soon as possible! Please also note that this code does not contain all the features described in the paper but rather gives a default set of functions to subdivide some given barycentric coordinates.
