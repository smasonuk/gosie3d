package gosie3d

import (
	"image/color"
	"math"
	"math/rand"
	"sort"

	"github.com/aquilax/go-perlin"
)

func NewCube() *Model {
	obj := NewObject_3d()
	s := 40.0 // size

	// The 8 vertices of the cube
	points := [][3]float64{
		{-s, -s, -s}, // 0
		{s, -s, -s},  // 1
		{s, s, -s},   // 2
		{-s, s, -s},  // 3
		{-s, -s, s},  // 4
		{s, -s, s},   // 5
		{s, s, s},    // 6
		{-s, s, s},   // 7
	}

	colors := []color.RGBA{
		{255, 0, 0, 255},   // Red
		{0, 255, 0, 255},   // Green
		{0, 0, 255, 255},   // Blue
		{255, 255, 0, 255}, // Yellow
		{0, 255, 255, 255}, // Cyan
		{255, 0, 255, 255}, // Magenta
	}

	// Define quads with a consistent Counter-Clockwise (CCW) winding order
	// This ensures the calculated normals point outwards from the cube.
	quads := [][]int{
		{0, 3, 2, 1}, // Front face  (Normal: 0, 0, -1)
		{1, 2, 6, 5}, // Right face  (Normal: 1, 0, 0)
		{5, 6, 7, 4}, // Back face   (Normal: 0, 0, 1)
		{4, 7, 3, 0}, // Left face   (Normal: -1, 0, 0)
		{3, 7, 6, 2}, // Top face    (Normal: 0, 1, 0)
		{4, 0, 1, 5}, // Bottom face (Normal: 0, -1, 0)
	}

	for i, q := range quads {
		face := NewFace(nil, colors[i], nil)
		// Vertices are added in the specified order to ensure correct normal
		face.AddPoint(points[q[0]][0], points[q[0]][1], points[q[0]][2])
		face.AddPoint(points[q[1]][0], points[q[1]][1], points[q[1]][2])
		face.AddPoint(points[q[2]][0], points[q[2]][1], points[q[2]][2])
		face.AddPoint(points[q[3]][0], points[q[3]][1], points[q[3]][2])

		// Use FACE_NORMAL because our winding order is now correct.
		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	obj.Finished(false, false)
	return obj
}

type Point2D struct {
	X, Z float64
}

// Extrude creates a 3D object by extruding a 2D shape defined by a set of
// vertices. It correctly handles vertices supplied in any order.
func Extrude(xp []float64, zp []float64, height float64, clr color.RGBA) *Model {
	// --- 1. Sanitize and Sort the 2D Points ---

	// Extract unique points from the input slices to avoid duplicates.
	pointSet := make(map[Point2D]bool)
	for i := 0; i < len(xp); i++ {
		pointSet[Point2D{X: xp[i], Z: zp[i]}] = true
	}

	// Create a slice from the unique points.
	var points []Point2D
	for p := range pointSet {
		points = append(points, p)
	}

	// Cannot form a polygon with less than 3 points.
	if len(points) < 3 {
		return NewObject_3d() // Return an empty object
	}

	// Calculate the centroid (average point) to sort around.
	var centerX, centerZ float64
	for _, p := range points {
		centerX += p.X
		centerZ += p.Z
	}
	centerX /= float64(len(points))
	centerZ /= float64(len(points))

	// Sort points by the angle they make with the centroid.
	// This arranges them into a continuous counter-clockwise (CCW) loop.
	sort.Slice(points, func(i, j int) bool {
		angle1 := math.Atan2(points[i].Z-centerZ, points[i].X-centerX)
		angle2 := math.Atan2(points[j].Z-centerZ, points[j].X-centerX)
		return angle1 < angle2
	})

	// --- 2. Build the 3D Object ---

	obj := NewObject_3d()
	topFace := NewFace(nil, clr, nil)
	baseFace := NewFace(nil, clr, nil)

	// Iterate through the sorted points to build all faces.
	for i := 0; i < len(points); i++ {
		// Get the start and end points for the current side panel.
		p1 := points[i]
		p2 := points[(i+1)%len(points)] // Wraps around for the last segment

		// Create the side face in CCW order for an outward normal.
		sideFace := NewFace(nil, clr, nil)
		sideFace.AddPoint(p1.X, 0, p1.Z)      // Bottom-start
		sideFace.AddPoint(p2.X, 0, p2.Z)      // Bottom-end
		sideFace.AddPoint(p2.X, height, p2.Z) // Top-end
		sideFace.AddPoint(p1.X, height, p1.Z) // Top-start
		sideFace.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(sideFace)

		// Add the current point to the top and base face polygons.
		topFace.AddPoint(p1.X, height, p1.Z)
		baseFace.AddPoint(p1.X, 0, p1.Z)
	}

	// Finalize the top and base cap faces.
	topFace.Finished(FACE_REVERSE) // CCW vertex order creates an upward (+Y) normal.
	// baseFace.Finished(FACE_REVERSE) // Must be reversed for a downward (-Y) normal.
	obj.theFaces.AddFace(topFace)
	// obj.theFaces.AddFace(baseFace)

	obj.Finished(true, true)
	return obj
}

// NewRectangle creates a new Object_3d in the shape of a cuboid (a 3D rectangle).
// It's centered at the origin and has the specified dimensions and color.
func NewRectangle(width, height, length float64, clr color.RGBA) *Model {
	// create a new, empty object to populate.
	obj := NewObject_3d()

	// calculate half-dimensions for centering the rectangle at the origin.
	w2 := width / 2.0
	h2 := height / 2.0
	l2 := length / 2.0

	// Define the 8 vertices of the cuboid based on the dimensions.
	points := [][3]float64{
		{-w2, -h2, -l2}, // 0: front-bottom-left
		{w2, -h2, -l2},  // 1: front-bottom-right
		{w2, h2, -l2},   // 2: front-top-right
		{-w2, h2, -l2},  // 3: front-top-left
		{-w2, -h2, l2},  // 4: back-bottom-left
		{w2, -h2, l2},   // 5: back-bottom-right
		{w2, h2, l2},    // 6: back-top-right
		{-w2, h2, l2},   // 7: back-top-left
	}

	quads := [][]int{
		{0, 3, 2, 1}, // Front face  (Normal: 0, 0, -1)
		{1, 2, 6, 5}, // Right face  (Normal: 1, 0, 0)
		{5, 6, 7, 4}, // Back face   (Normal: 0, 0, 1)
		{4, 7, 3, 0}, // Left face   (Normal: -1, 0, 0)
		{3, 7, 6, 2}, // Top face    (Normal: 0, 1, 0)
		{4, 0, 1, 5}, // Bottom face (Normal: 0, -1, 0)
	}

	// Create each face and add it to the object.
	for _, q := range quads {
		// Use the single color passed into the function for every face.
		face := NewFace(nil, clr, nil)

		// Add the four vertices that make up this face.
		face.AddPoint(points[q[0]][0], points[q[0]][1], points[q[0]][2])
		face.AddPoint(points[q[1]][0], points[q[1]][1], points[q[1]][2])
		face.AddPoint(points[q[2]][0], points[q[2]][1], points[q[2]][2])
		face.AddPoint(points[q[3]][0], points[q[3]][1], points[q[3]][2])

		// Finalize the face's geometry. Use FACE_NORMAL because the winding
		// order is defined correctly to produce an outward normal.
		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	obj.Finished(false, false)

	return obj
}

// NewSubdividedRectangle creates a new Object_3d in the shape of a cuboid where each face
// is subdivided into a grid of smaller triangles.
// width, height, length: The dimensions of the cuboid.
// clr: The color for all faces of the cuboid.
// subdivisions: The number of divisions along each edge of a face. For example, a value
//
//	of 2 will split a face into a 2x2 grid of quads (8 triangles). A value of 1
//	will result in one quad per face (2 triangles).
func NewSubdividedRectangle(width, height, length float64, clr color.RGBA, subdivisions int) *Model {
	obj := NewObject_3d()
	w2, h2, l2 := width/2.0, height/2.0, length/2.0

	if subdivisions < 1 {
		subdivisions = 1 // Ensure at least one subdivision.
	}

	// generateFace is a helper function that constructs one of the six faces of the cuboid.
	// It takes an origin point and two vectors (u, v) that define the plane and dimensions
	// of the face. It then creates a grid of vertices and generates triangles.
	generateFace := func(origin, u, v [3]float64) {
		// Create a grid of vertices for the current face.
		vertices := make([][][3]float64, subdivisions+1)
		for i := range vertices {
			vertices[i] = make([][3]float64, subdivisions+1)
			for j := range vertices[i] {
				// Calculate the position of vertex (i, j) on the grid plane.
				ui := float64(i) / float64(subdivisions)
				vj := float64(j) / float64(subdivisions)
				vertices[i][j] = [3]float64{
					origin[0] + ui*u[0] + vj*v[0],
					origin[1] + ui*u[1] + vj*v[1],
					origin[2] + ui*u[2] + vj*v[2],
				}
			}
		}

		// Create two triangles for each quad in the subdivision grid.
		for i := 0; i < subdivisions; i++ {
			for j := 0; j < subdivisions; j++ {
				// Get the four corner vertices of the current quad.
				p1 := vertices[i][j]
				p2 := vertices[i+1][j]
				p3 := vertices[i+1][j+1]
				p4 := vertices[i][j+1]

				// Create the first triangle for the quad (p1, p2, p3).
				face1 := NewFace(nil, clr, nil)
				face1.AddPoint(p1[0], p1[1], p1[2])
				face1.AddPoint(p2[0], p2[1], p2[2])
				face1.AddPoint(p3[0], p3[1], p3[2])
				// The vertices are wound counter-clockwise to produce an outward-facing normal.
				face1.Finished(FACE_NORMAL)
				obj.theFaces.AddFace(face1)

				// Create the second triangle for the quad (p1, p3, p4).
				face2 := NewFace(nil, clr, nil)
				face2.AddPoint(p1[0], p1[1], p1[2])
				face2.AddPoint(p3[0], p3[1], p3[2])
				face2.AddPoint(p4[0], p4[1], p4[2])
				face2.Finished(FACE_NORMAL)
				obj.theFaces.AddFace(face2)
			}
		}
	}

	// Define the 6 faces of the cuboid by specifying their origin, u-vector, and v-vector.
	// The u and v vectors are chosen so their cross-product (which determines the normal)
	// points outwards from the center of the cuboid.

	// Back face (-Z direction)
	generateFace([3]float64{w2, -h2, -l2}, [3]float64{-width, 0, 0}, [3]float64{0, height, 0})

	// Front face (+Z direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{width, 0, 0}, [3]float64{0, height, 0})

	// Left face (-X direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{0, 0, -length}, [3]float64{0, height, 0})

	// Right face (+X direction)
	generateFace([3]float64{w2, -h2, -l2}, [3]float64{0, 0, length}, [3]float64{0, height, 0})

	// Bottom face (-Y direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{width, 0, 0}, [3]float64{0, 0, -length})

	// Top face (+Y direction)
	generateFace([3]float64{-w2, h2, -l2}, [3]float64{width, 0, 0}, [3]float64{0, 0, length})

	// Finalize the object by building its BSP tree.
	obj.Finished(false, false)
	return obj
}

func gen(xWidth, yLength float64, clr color.RGBA, subdivisions int, useTriangles bool) *Model {
	obj := NewObject_3d()
	w2, l2 := xWidth/2.0, yLength/2.0

	if subdivisions < 1 {
		subdivisions = 1 // Ensure at least one subdivision.
	}

	// generateFace is a helper function that constructs one of the six faces of the cuboid.
	// It takes an origin point and two vectors (u, v) that define the plane and dimensions
	// of the face. It then creates a grid of vertices and generates triangles.
	generateFace := func(origin, u, v [3]float64) {
		// Create a grid of vertices for the current face.
		vertices := make([][][3]float64, subdivisions+1)
		for i := range vertices {
			vertices[i] = make([][3]float64, subdivisions+1)
			for j := range vertices[i] {
				// Calculate the position of vertex (i, j) on the grid plane.
				ui := float64(i) / float64(subdivisions)
				vj := float64(j) / float64(subdivisions)
				vertices[i][j] = [3]float64{
					origin[0] + ui*u[0] + vj*v[0],
					origin[1] + ui*u[1] + vj*v[1],
					origin[2] + ui*u[2] + vj*v[2],
				}
			}
		}

		// Create two triangles for each quad in the subdivision grid.
		for i := 0; i < subdivisions; i++ {
			for j := 0; j < subdivisions; j++ {
				// Get the four corner vertices of the current quad.
				p1 := vertices[i][j]
				p2 := vertices[i+1][j]
				p3 := vertices[i+1][j+1]
				p4 := vertices[i][j+1]

				if useTriangles {
					// Create the first triangle for the quad (p1, p2, p3).
					face1 := NewFace(nil, clr, nil)
					face1.AddPoint(p1[0], p1[1], p1[2])
					face1.AddPoint(p2[0], p2[1], p2[2])
					face1.AddPoint(p3[0], p3[1], p3[2])
					// The vertices are wound counter-clockwise to produce an outward-facing normal.
					face1.Finished(FACE_REVERSE)
					obj.theFaces.AddFace(face1)

					// Create the second triangle for the quad (p1, p3, p4).
					face2 := NewFace(nil, clr, nil)
					face2.AddPoint(p1[0], p1[1], p1[2])
					face2.AddPoint(p3[0], p3[1], p3[2])
					face2.AddPoint(p4[0], p4[1], p4[2])
					face2.Finished(FACE_REVERSE)
					obj.theFaces.AddFace(face2)

				} else {
					// Create the first triangle for the quad (p1, p2, p3).
					face1 := NewFace(nil, clr, nil)
					face1.AddPoint(p1[0], p1[1], p1[2])
					face1.AddPoint(p2[0], p2[1], p2[2])
					face1.AddPoint(p3[0], p3[1], p3[2])
					face1.AddPoint(p4[0], p4[1], p4[2])
					// The vertices are wound counter-clockwise to produce an outward-facing normal.
					face1.Finished(FACE_REVERSE)
					obj.theFaces.AddFace(face1)
				}
			}
		}
	}

	// Top face (+Y direction)
	generateFace([3]float64{-w2, 0, -l2}, [3]float64{xWidth, 0, 0}, [3]float64{0, 0, yLength})

	return obj
}

func NewSubdividedPlane(xWidth, yLength float64, clr color.RGBA, subdivisions int, useTriangles bool) *Model {
	obj := gen(xWidth, yLength, clr, subdivisions, useTriangles)

	obj.SetDrawAllFaces(true)

	// Finalize the object by building its BSP tree.
	obj.Finished(false, false)
	return obj
}

func NewSubdividedPlaneHeightMap(xWidth,
	yLength float64,
	clr color.RGBA,
	subdivisions int,
	flatInX float64,
	flatInY float64,
) *Model {
	obj := NewObject_3d()
	w2, l2 := xWidth/2.0, yLength/2.0

	if subdivisions < 1 {
		subdivisions = 1
	}

	// generateFace is a helper function that constructs one of the six faces of the cuboid.
	// It takes an origin point and two vectors (u, v) that define the plane and dimensions
	// of the face. It then creates a grid of vertices and generates triangles.
	generateFace := func(origin, u, v [3]float64) {

		// Create a grid of vertices for the current face.
		vertices := make([][][3]float64, subdivisions+1)
		for i := range vertices {
			vertices[i] = make([][3]float64, subdivisions+1)
			for j := range vertices[i] {
				heightAdjust := 0.0
				// atEdgeOfPlane := (i == 0 || i == subdivisions || j == 0 || j == subdivisions)
				// if !atEdgeOfPlane {

				ui := float64(i) / float64(subdivisions)
				vj := float64(j) / float64(subdivisions)
				x := origin[0] + ui*u[0] + vj*v[0]
				z := origin[2] + ui*u[2] + vj*v[2]

				if x > 0 && x < flatInX && z > 0 && z < flatInY {
					heightAdjust = 0.0 // Flat area in the middle
				} else {

					// rndHeight := rand.Float64()
					// heightAdjust = rndHeight * 50.0

					dist := math.Sqrt((x*x + z*z))
					// // heightAdjust += (math.Sin(dist/100.0) * 20.0) // Add some wave-like variation
					// heightAdjust = heightAdjust * (dist / 400.0)
					// heightAdjust = (math.Sin(dist/100.0) * (dist / 10.0))
					heightAdjust = ((math.Sin(x/200.0) + math.Sin(z/200.0)) * (dist / 25.0))
					// + (rndHeight * 20)
				}

				vertices[i][j] = [3]float64{
					x,
					(origin[1] + ui*u[1] + vj*v[1]) - heightAdjust,
					z,
				}

			}
		}

		// // Create two triangles for each quad in the subdivision grid.
		// for i := 0; i < subdivisions; i++ {
		// 	for j := 0; j < subdivisions; j++ {
		// 		// Get the four corner vertices of the current quad.
		// 		p1 := vertices[i][j]
		// 		p2 := vertices[i+1][j]
		// 		p3 := vertices[i+1][j+1]
		// 		p4 := vertices[i][j+1]

		// 		heightAdjust := 0.0

		// 		// Create the first triangle for the quad (p1, p2, p3).
		// 		face1 := NewFace(nil, clr, nil)
		// 		face1.AddPoint(p1[0], p1[1]+heightAdjust, p1[2])
		// 		face1.AddPoint(p2[0], p2[1]+heightAdjust, p2[2])
		// 		face1.AddPoint(p3[0], p3[1]+heightAdjust, p3[2])
		// 		// The vertices are wound counter-clockwise to produce an outward-facing normal.
		// 		face1.Finished(FACE_REVERSE)
		// 		obj.theFaces.AddFace(face1)

		// 		// Create the second triangle for the quad (p1, p3, p4).
		// 		face2 := NewFace(nil, clr, nil)
		// 		face2.AddPoint(p1[0], p1[1]+heightAdjust, p1[2])
		// 		face2.AddPoint(p3[0], p3[1]+heightAdjust, p3[2])
		// 		face2.AddPoint(p4[0], p4[1]+heightAdjust, p4[2])
		// 		face2.Finished(FACE_REVERSE)
		// 		obj.theFaces.AddFace(face2)
		// 	}
		// }

		for i := 0; i < subdivisions; i++ {
			for j := 0; j < subdivisions; j++ {
				// Get the four corner vertices of the current quad.
				p1 := vertices[i][j]
				p2 := vertices[i+1][j]
				p3 := vertices[i+1][j+1]
				p4 := vertices[i][j+1]

				// Create the first triangle for the quad (p1, p2, p3).
				face1 := NewFace(nil, clr, nil)

				face1.AddPoint(p1[0], p1[1], p1[2])

				face1.AddPoint(p2[0], p2[1], p2[2])

				face1.AddPoint(p3[0], p3[1], p3[2])
				face1.AddPoint(p4[0], p4[1], p4[2])

				face1.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(face1)

			}
		}

		obj.theFaces.SortFacesByDistance(NewVector3(0, 0, 0))

	}

	// Top face (+Y direction)
	generateFace([3]float64{-w2, 0, -l2}, [3]float64{xWidth, 0, 0}, [3]float64{0, 0, yLength})

	// Finalize the object by building its BSP tree.
	obj.Finished(false, false)
	return obj
}

func NewSubdividedPlaneHeightMapPerlin(xWidth,
	yLength float64,
	clr color.RGBA,
	subdivisions int,
	flatInX float64,
	flatInY float64,
) *Model {
	obj := NewObject_3d()
	w2, l2 := xWidth/2.0, yLength/2.0

	if subdivisions < 1 {
		subdivisions = 1
	}
	// Seed the random number generator for reproducibility.
	randomSeed := rand.Int63()

	per := perlin.NewPerlin(2, 3, 1, randomSeed)

	// generateFace is a helper function that constructs one of the six faces of the cuboid.
	// It takes an origin point and two vectors (u, v) that define the plane and dimensions
	// of the face. It then creates a grid of vertices and generates triangles.
	generateFace := func(origin, u, v [3]float64) {

		// Create a grid of vertices for the current face.
		vertices := make([][][3]float64, subdivisions+1)
		for i := range vertices {
			vertices[i] = make([][3]float64, subdivisions+1)
			for j := range vertices[i] {
				heightAdjust := 0.0
				// atEdgeOfPlane := (i == 0 || i == subdivisions || j == 0 || j == subdivisions)
				// if !atEdgeOfPlane {

				ui := float64(i) / float64(subdivisions)
				vj := float64(j) / float64(subdivisions)
				x := origin[0] + ui*u[0] + vj*v[0]
				z := origin[2] + ui*u[2] + vj*v[2]

				// if x > -flatInX && x < flatInX && z > -flatInY && z < flatInY {
				// 	// heightAdjust = 0.0 // Flat area in the middle

				// 	heightAdjust = per.Noise2D(x/1000.0, z/1000.0) * 200.0
				// } else {

				distFromZero := math.Abs(z)

				if distFromZero > 700.0 {
					distFromZero = 700.0
				}

				// heightAdjust = per.Noise2D(x/1000.0, z/1000.0) * 700.0 // Perlin noise for height variation

				heightAdjust = per.Noise2D(x/1000.0, z/1000.0) * distFromZero

				// }

				vertices[i][j] = [3]float64{
					x,
					(origin[1] + ui*u[1] + vj*v[1]) - heightAdjust,
					z,
				}

			}
		}

		// Create two triangles for each quad in the subdivision grid.
		for i := 0; i < subdivisions; i++ {
			for j := 0; j < subdivisions; j++ {
				// Get the four corner vertices of the current quad.
				p1 := vertices[i][j]
				p2 := vertices[i+1][j]
				p3 := vertices[i+1][j+1]
				p4 := vertices[i][j+1]

				heightAdjust := 0.0

				// Create the first triangle for the quad (p1, p2, p3).
				face1 := NewFace(nil, clr, nil)
				face1.AddPoint(p1[0], p1[1]+heightAdjust, p1[2])
				face1.AddPoint(p2[0], p2[1]+heightAdjust, p2[2])
				face1.AddPoint(p3[0], p3[1]+heightAdjust, p3[2])
				// The vertices are wound counter-clockwise to produce an outward-facing normal.
				face1.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(face1)

				// Create the second triangle for the quad (p1, p3, p4).
				face2 := NewFace(nil, clr, nil)
				face2.AddPoint(p1[0], p1[1]+heightAdjust, p1[2])
				face2.AddPoint(p3[0], p3[1]+heightAdjust, p3[2])
				face2.AddPoint(p4[0], p4[1]+heightAdjust, p4[2])
				face2.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(face2)
			}
		}

		// for i := 0; i < subdivisions; i++ {
		// 	for j := 0; j < subdivisions; j++ {
		// 		// Get the four corner vertices of the current quad.
		// 		p1 := vertices[i][j]
		// 		p2 := vertices[i+1][j]
		// 		p3 := vertices[i+1][j+1]
		// 		p4 := vertices[i][j+1]

		// 		// Create the first triangle for the quad (p1, p2, p3).
		// 		face1 := NewFace(nil, clr, nil)

		// 		face1.AddPoint(p1[0], p1[1], p1[2])

		// 		face1.AddPoint(p2[0], p2[1], p2[2])

		// 		face1.AddPoint(p3[0], p3[1], p3[2])
		// 		face1.AddPoint(p4[0], p4[1], p4[2])

		// 		face1.Finished(FACE_REVERSE)
		// 		obj.theFaces.AddFace(face1)

		// 	}
		// }

		obj.theFaces.SortFacesByDistance(NewVector3(0, 0, 0))

	}

	// Top face (+Y direction)
	generateFace([3]float64{-w2, 0, -l2}, [3]float64{xWidth, 0, 0}, [3]float64{0, 0, yLength})

	// Finalize the object by building its BSP tree.
	obj.Finished(false, false)
	return obj
}

// NewSphere creates a new Object_3d in the shape of an icosphere.
// An icosphere is a sphere made of a mesh of triangles, which is more
// uniform than a traditional UV sphere.
func NewSphere(radius float64, subdivisions int, clr color.RGBA, finish bool) *Model {
	obj := NewObject_3d()

	// Define the 12 vertices of an Icosahedron.
	// An icosahedron is a 20-sided polyhedron that forms the base of our sphere.
	// The 't' value is the golden ratio, which helps define the vertex positions.
	t := (1.0 + math.Sqrt(5.0)) / 2.0

	vertices := [][3]float64{
		{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0},
		{0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
		{t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1},
	}

	// Define the 20 triangular faces of the Icosahedron using indices
	// into the vertex list. The order is important for correct normals.
	faces := [][]int{
		{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
		{1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
		{3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
		{4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1},
	}

	// Create initial list of faces for subdivision
	subdivisionFaces := make([][3][3]float64, len(faces))
	for i, faceIndices := range faces {
		subdivisionFaces[i] = [3][3]float64{
			vertices[faceIndices[0]],
			vertices[faceIndices[1]],
			vertices[faceIndices[2]],
		}
	}

	// Subdivide the faces recursively to make the sphere smoother.
	for i := 0; i < subdivisions; i++ {
		newFaces := make([][3][3]float64, 0)
		for _, face := range subdivisionFaces {
			v1 := face[0]
			v2 := face[1]
			v3 := face[2]

			// Calculate the midpoint of each edge of the triangle.
			a := [3]float64{(v1[0] + v2[0]) / 2, (v1[1] + v2[1]) / 2, (v1[2] + v2[2]) / 2}
			b := [3]float64{(v2[0] + v3[0]) / 2, (v2[1] + v3[1]) / 2, (v2[2] + v3[2]) / 2}
			c := [3]float64{(v3[0] + v1[0]) / 2, (v3[1] + v1[1]) / 2, (v3[2] + v1[2]) / 2}

			// Split the original triangle into 4 new ones.
			newFaces = append(newFaces, [3][3]float64{v1, a, c})
			newFaces = append(newFaces, [3][3]float64{v2, b, a})
			newFaces = append(newFaces, [3][3]float64{v3, c, b})
			newFaces = append(newFaces, [3][3]float64{a, b, c})
		}
		subdivisionFaces = newFaces
	}

	// Create the final faces for the Object_3d
	for _, faceVerts := range subdivisionFaces {
		face := NewFace(nil, clr, nil)

		for _, v := range faceVerts {
			// Normalize the vertex to push it onto the surface of a unit sphere.
			length := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
			normX := v[0] / length
			normY := v[1] / length
			normZ := v[2] / length

			// Scale by the desired radius and add the point.
			face.AddPoint(normX*radius, normY*radius, normZ*radius)
		}

		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	if finish {
		obj.Finished(true, false)
	}
	return obj
}

// NewUVSphere creates a sphere based on latitude/longitude rings (sectors and stacks).
// This structure allows for a perfectly straight horizontal stripe.
func NewUVSphere(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Model {
	obj := NewObject_3d()

	// We loop through stacks (latitude) and sectors (longitude).
	vertices := make([][3]float64, 0)
	for i := 0; i <= stacks; i++ {
		stackAngle := math.Pi/2 - float64(i)*math.Pi/float64(stacks) // phi
		xy := radius * math.Cos(stackAngle)
		z := radius * math.Sin(stackAngle)

		for j := 0; j <= sectors; j++ {
			sectorAngle := float64(j) * 2 * math.Pi / float64(sectors) // theta
			x := xy * math.Cos(sectorAngle)
			y := xy * math.Sin(sectorAngle)
			vertices = append(vertices, [3]float64{x, y, z})
		}
	}

	// Determine the start and end stacks for the stripe.
	// The stripe is centered around the equator (the middle stack).
	middleStack := stacks / 2
	stripeStart := middleStack - (stripeStacks / 2)
	stripeEnd := middleStack + (stripeStacks / 2)

	for i := 0; i < stacks; i++ {
		// Determine the color for this entire ring of faces.
		var faceColor color.RGBA
		if i >= stripeStart && i < stripeEnd {
			faceColor = stripeClr
		} else {
			faceColor = bodyClr
		}

		k1 := i * (sectors + 1)
		k2 := k1 + sectors + 1

		for j := 0; j < sectors; j++ {
			k1j := k1 + j
			k1j1 := k1j + 1
			k2j := k2 + j
			k2j1 := k2j + 1

			// For each quad, we create two triangles.
			// Special handling for the poles.
			if i != 0 {
				// First triangle of the quad
				f1 := NewFace(nil, faceColor, nil)
				f1.AddPoint(vertices[k1j][0], vertices[k1j][1], vertices[k1j][2])
				f1.AddPoint(vertices[k2j][0], vertices[k2j][1], vertices[k2j][2])
				f1.AddPoint(vertices[k1j1][0], vertices[k1j1][1], vertices[k1j1][2])
				f1.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(f1)
			}

			if i != (stacks - 1) {
				// Second triangle of the quad
				f2 := NewFace(nil, faceColor, nil)
				f2.AddPoint(vertices[k1j1][0], vertices[k1j1][1], vertices[k1j1][2])
				f2.AddPoint(vertices[k2j][0], vertices[k2j][1], vertices[k2j][2])
				f2.AddPoint(vertices[k2j1][0], vertices[k2j1][1], vertices[k2j1][2])
				f2.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(f2)
			}
		}
	}

	// 3. Finalize the object by building its BSP tree.
	obj.Finished(false, false)
	return obj
}

func NewUVSphere2(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Model {
	obj := NewObject_3d()

	vertices := make([][3]float64, 0)
	for i := 0; i <= stacks; i++ {
		stackAngle := math.Pi/2 - float64(i)*math.Pi/float64(stacks)
		xy := radius * math.Cos(stackAngle)
		z := radius * math.Sin(stackAngle)

		for j := 0; j <= sectors; j++ {
			sectorAngle := float64(j) * 2 * math.Pi / float64(sectors)
			x := xy * math.Cos(sectorAngle)
			y := xy * math.Sin(sectorAngle)
			vertices = append(vertices, [3]float64{x, y, z})
		}
	}

	middleStack := stacks / 2
	stripeStart := middleStack - (stripeStacks / 2)
	stripeEnd := middleStack + (stripeStacks / 2)

	// A helper function to create a face and set its normal correctly
	createFace := func(p1, p2, p3 [3]float64, clr color.RGBA) *Face {
		face := NewFace(nil, clr, nil)
		face.AddPoint(p1[0], p1[1], p1[2])
		face.AddPoint(p2[0], p2[1], p2[2])
		face.AddPoint(p3[0], p3[1], p3[2])

		// Find the center of the triangle face.
		centerX := (p1[0] + p2[0] + p3[0]) / 3.0
		centerY := (p1[1] + p2[1] + p3[1]) / 3.0
		centerZ := (p1[2] + p2[2] + p3[2]) / 3.0

		// The normal is the vector from the sphere's center (0,0,0) to the face's center, normalized.
		normalVec := []float64{centerX, centerY, centerZ, 0.0}
		length := GetLength(normalVec)
		if length > 0 {
			normalVec[0] /= length
			normalVec[1] /= length
			normalVec[2] /= length
		}

		normalVec[0] *= -1
		normalVec[1] *= -1
		normalVec[2] *= -1

		// Explicitly set this perfect normal on the face.
		face.SetNormal(NewVector3(normalVec[0], normalVec[1], normalVec[2]))
		// We no longer need the FACE_NORMAL/FACE_REVERSE flag.
		face.Finished(FACE_REVERSE)
		return face
	}

	for i := 0; i < stacks; i++ {
		var faceColor color.RGBA
		if i >= stripeStart && i < stripeEnd {
			faceColor = stripeClr
		} else {
			faceColor = bodyClr
		}

		k1 := i * (sectors + 1)
		k2 := k1 + sectors + 1

		for j := 0; j < sectors; j++ {
			// Get the four vertices for the quad on the sphere surface
			v1 := vertices[k1+j]
			v2 := vertices[k2+j]
			v3 := vertices[k1+j+1]
			v4 := vertices[k2+j+1]

			if i != 0 {
				// First triangle of the quad: v1, v2, v3
				obj.theFaces.AddFace(createFace(v1, v2, v3, faceColor))
			}

			if i != (stacks - 1) {
				// Second triangle of the quad: v3, v2, v4
				obj.theFaces.AddFace(createFace(v3, v2, v4, faceColor))
			}
		}
	}

	obj.Finished(false, false)
	return obj
}

// NewCylinder creates a new Object3d in the shape of a cylinder.
// The cylinder is centered at the origin, and its height extends along the Y-axis.
// 'segments' controls how many rectangular faces are used to approximate the curved surface;
// a higher number results in a smoother cylinder.
func NewCylinder(radius, height float64, segments int, clr color.RGBA) *Model {
	// A cylinder must have at least 3 segments to form a valid polygonal base.
	if segments < 3 {
		return NewObject_3d() // Return an empty object if segments are insufficient.
	}

	obj := NewObject_3d()
	topFace := NewFace(nil, clr, nil)
	baseFace := NewFace(nil, clr, nil)

	// --- 1. Generate Vertices and Build Cap Faces ---

	// Create slices to hold the ordered vertices for the top and bottom caps.
	topVertices := make([][3]float64, segments)
	baseVertices := make([][3]float64, segments)

	// Calculate the vertex positions for one full circle.
	for i := 0; i < segments; i++ {
		angle := (2.0 * math.Pi / float64(segments)) * float64(i)
		x := radius * math.Cos(angle)
		z := radius * math.Sin(angle)

		// Store vertices for later use when building the side walls.
		topVertices[i] = [3]float64{x, height, z}
		baseVertices[i] = [3]float64{x, 0, z}

		// Add the calculated vertex to the top face polygon.
		// The points are added in counter-clockwise (CCW) order, which will
		// result in a normal vector pointing up (+Y) after finalization.
		topFace.AddPoint(x, height, z)
	}

	// For the base face, vertices must be in clockwise (CW) order to produce a
	// downward (-Y) normal. We achieve this by adding the pre-calculated
	// base vertices to the face in reverse order.
	for i := segments - 1; i >= 0; i-- {
		v := baseVertices[i]
		baseFace.AddPoint(v[0], v[1], v[2])
	}

	// Finalize the cap faces and add them to the object.
	topFace.Finished(FACE_NORMAL)
	obj.theFaces.AddFace(topFace)

	baseFace.Finished(FACE_NORMAL)
	obj.theFaces.AddFace(baseFace)

	// --- 2. Build the Side Wall Faces ---

	// Create a quad for each segment connecting the top and bottom caps.
	for i := 0; i < segments; i++ {
		// Get the four corners of the quad. The modulo operator (%) ensures
		// that the last vertex connects back to the first one.
		p1_base := baseVertices[i]
		p2_base := baseVertices[(i+1)%segments]
		p1_top := topVertices[i]
		p2_top := topVertices[(i+1)%segments]

		// Create the side face with vertices in CCW order for an outward normal.
		sideFace := NewFace(nil, clr, nil)
		sideFace.AddPoint(p1_base[0], p1_base[1], p1_base[2]) // Bottom-start
		sideFace.AddPoint(p2_base[0], p2_base[1], p2_base[2]) // Bottom-end
		sideFace.AddPoint(p2_top[0], p2_top[1], p2_top[2])    // Top-end
		sideFace.AddPoint(p1_top[0], p1_top[1], p1_top[2])    // Top-start

		sideFace.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(sideFace)
	}

	obj.Finished(true, false)
	obj.CentreObject()
	return obj
}

// Assume the following types from your project are available:
// type Object3d struct { /* ... */ }
// func NewObject_3d(bool) *Object3d { /* ... */ }
// func (o *Object3d) Finished(bool) { /* ... */ }
// func (o *Object3d) CentreObject() { /* ... */ }
//
// type Face struct { /* ... */ }
// func NewFace(nil, color.RGBA, nil) *Face { /* ... */ }
// func (f *Face) AddPoint(x, y, z float64) { /* ... */ }
// func (f *Face) Finished(string) { /* ... */ }
//
// const FACE_NORMAL = "FACE_NORMAL" // Or whatever the constant is.

// NewRing creates a new Object3d in the shape of a ring or pipe.
// The ring is initially built with its base at Y=0 and is then centered at the origin.
// 'segments' controls the smoothness of the curved surfaces.
// 'outerRadius' and 'innerRadius' define the ring's thickness.
func NewRing(outerRadius, innerRadius, height float64, segments int, clr color.RGBA, finish bool) *Model {
	// A ring must have at least 3 segments and the outer radius must be larger than the inner.
	if segments < 3 || outerRadius <= innerRadius {
		return NewObject_3d() // Return an empty object if parameters are invalid.
	}

	obj := NewObject_3d()

	// --- 1. Generate Vertices ---

	// Create slices to hold the ordered vertices for the four circular edges of the ring.
	outerTopVertices := make([][3]float64, segments)
	innerTopVertices := make([][3]float64, segments)
	outerBaseVertices := make([][3]float64, segments)
	innerBaseVertices := make([][3]float64, segments)

	// Calculate vertex positions for one full circle at both radii.
	for i := 0; i < segments; i++ {
		angle := (2.0 * math.Pi / float64(segments)) * float64(i)
		cosAngle := math.Cos(angle)
		sinAngle := math.Sin(angle)

		// Outer vertices
		ox, oz := outerRadius*cosAngle, outerRadius*sinAngle
		outerTopVertices[i] = [3]float64{ox, height, oz}
		outerBaseVertices[i] = [3]float64{ox, 0, oz}

		// Inner vertices
		ix, iz := innerRadius*cosAngle, innerRadius*sinAngle
		innerTopVertices[i] = [3]float64{ix, height, iz}
		innerBaseVertices[i] = [3]float64{ix, 0, iz}
	}

	// --- 2. Build Faces ---

	// Create faces for each segment connecting the inner and outer vertices.
	for i := 0; i < segments; i++ {
		// Indices for the current segment and the next, wrapping around.
		i2 := (i + 1) % segments

		// Get the four corners of the quad for the top surface of this segment.
		ot1 := outerTopVertices[i]  // Outer Top, point 1
		ot2 := outerTopVertices[i2] // Outer Top, point 2
		it1 := innerTopVertices[i]  // Inner Top, point 1
		it2 := innerTopVertices[i2] // Inner Top, point 2

		// Get the four corners for the bottom surface.
		ob1 := outerBaseVertices[i]  // Outer Base, point 1
		ob2 := outerBaseVertices[i2] // Outer Base, point 2
		ib1 := innerBaseVertices[i]  // Inner Base, point 1
		ib2 := innerBaseVertices[i2] // Inner Base, point 2

		// --- Top Face ---
		// Create the top face with vertices in CCW order for an upward (+Y) normal.
		topFace := NewFace(nil, clr, nil)
		topFace.AddPoint(ot1[0], ot1[1], ot1[2])
		topFace.AddPoint(ot2[0], ot2[1], ot2[2])
		topFace.AddPoint(it2[0], it2[1], it2[2])
		topFace.AddPoint(it1[0], it1[1], it1[2])
		topFace.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(topFace)

		// --- Bottom Face ---
		// Create the bottom face with vertices in CW order for a downward (-Y) normal.
		bottomFace := NewFace(nil, clr, nil)
		bottomFace.AddPoint(ob1[0], ob1[1], ob1[2])
		bottomFace.AddPoint(ib1[0], ib1[1], ib1[2])
		bottomFace.AddPoint(ib2[0], ib2[1], ib2[2])
		bottomFace.AddPoint(ob2[0], ob2[1], ob2[2])
		bottomFace.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(bottomFace)

		// --- Outer Wall Face ---
		// Create the outer wall with CCW order for an outward-pointing normal.
		outerWallFace := NewFace(nil, clr, nil)
		outerWallFace.AddPoint(ob1[0], ob1[1], ob1[2])
		outerWallFace.AddPoint(ob2[0], ob2[1], ob2[2])
		outerWallFace.AddPoint(ot2[0], ot2[1], ot2[2])
		outerWallFace.AddPoint(ot1[0], ot1[1], ot1[2])
		outerWallFace.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(outerWallFace)

		// --- Inner Wall Face ---
		// Create the inner wall. The vertex order is reversed (CW from outside,
		// but CCW from inside) to make the normal point inward.
		innerWallFace := NewFace(nil, clr, nil)
		innerWallFace.AddPoint(ib1[0], ib1[1], ib1[2])
		innerWallFace.AddPoint(it1[0], it1[1], it1[2])
		innerWallFace.AddPoint(it2[0], it2[1], it2[2])
		innerWallFace.AddPoint(ib2[0], ib2[1], ib2[2])
		innerWallFace.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(innerWallFace)
	}

	if finish {
		obj.Finished(true, false)
	}
	return obj
}
