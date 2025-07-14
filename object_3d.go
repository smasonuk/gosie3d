package gosie3d

import (
	"bufio"
	"fmt"
	"image/color"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"github.com/hajimehoshi/ebiten/v2"
)

func (o *Object3d) PaintSolid(screen *ebiten.Image, x, y int, lightingChange bool) {
	if o.root != nil {
		o.root.PaintWithShading(screen, x, y, o.transFaceMesh.Points, o.transNormalMesh.Points, lightingChange)
	}
}

type Object3d struct {
	faceMesh        *FaceMesh
	normalMesh      *NormalMesh
	transFaceMesh   *FaceMesh
	transNormalMesh *NormalMesh
	theFaces        *FaceStore
	root            *BspNode
	rotMatrix       *Matrix
	position        *Point3d
}

func NewObject_3d() *Object3d {
	return &Object3d{
		transFaceMesh:   NewFaceMesh(),
		transNormalMesh: NewNormalMesh(),
		theFaces:        NewFaceStore(),
		rotMatrix:       IdentMatrix(),
	}
}

func (o *Object3d) Clone() *Object3d {
	clone := &Object3d{
		// shared
		faceMesh:   o.faceMesh,
		normalMesh: o.normalMesh,
		theFaces:   o.theFaces,
		root:       o.root,

		// instance-specific
		transFaceMesh:   o.transFaceMesh.Copy(),
		transNormalMesh: o.transNormalMesh.Copy(),
		rotMatrix:       IdentMatrix(),
	}
	return clone
}

func (o *Object3d) Finished() {
	log.Println("Creating BSP Tree...")
	if o.theFaces.FaceCount() > 0 {
		o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
	}
	log.Println("BSP Tree Created.")

	o.faceMesh = &FaceMesh{Mesh: *o.transFaceMesh.Mesh.Copy()}
	o.normalMesh = &NormalMesh{Mesh: *o.transNormalMesh.Mesh.Copy()}

	log.Printf("Points: %d", len(o.faceMesh.Points.ThisMatrix))
	log.Printf("Normals: %d", len(o.normalMesh.Points.ThisMatrix))
}

func (o *Object3d) ApplyMatrixBatch(m *Matrix) {
	o.rotMatrix = m.MultiplyBy(o.rotMatrix)
}

func (o *Object3d) ApplyMatrixTemp(aMatrix *Matrix) {
	rotMatrixTemp := aMatrix.MultiplyBy(o.rotMatrix)

	// Use the new, correct method to transform the normals (rotation only).
	rotMatrixTemp.TransformNormals(o.normalMesh.Points, o.transNormalMesh.Points)

	// Use the original method to transform the vertex positions (rotation and translation).
	rotMatrixTemp.TransformObj(o.faceMesh.Points, o.transFaceMesh.Points)
}

func (o *Object3d) ApplyMatrixPermanent(aMatrix *Matrix) {
	o.ApplyMatrixTemp(aMatrix)

	o.normalMesh.Points = o.transNormalMesh.Points.Copy()
	o.faceMesh.Points = o.transFaceMesh.Points.Copy()
}

func (o *Object3d) createBspTree(faces *FaceStore, newFaces *FaceMesh, newNormMesh *NormalMesh) *BspNode {
	if faces.FaceCount() == 0 {
		return nil
	}

	parentFace := o.choosePlane(faces)
	originalNormal, normalIndex := newNormMesh.AddNormal(parentFace.GetNormal())
	parentFace.SetNormal(NewVector3(originalNormal[0], originalNormal[1], originalNormal[2]))
	newFace, parentIndices := newFaces.AddFace(parentFace)
	parent := NewBspNode(newFace.Points, newFace.GetNormal(), newFace.Col, parentIndices, normalIndex)
	pPlane := NewPlane(newFace, newFace.GetNormal())

	// Create two new lists to hold the faces that fall on either side of the plane.
	fvLeft := NewFaceStore()
	fvRight := NewFaceStore()

	// Partition the *remaining* faces against the splitting plane.
	for a := 0; a < faces.FaceCount(); a++ {
		currentFace := faces.GetFace(a)
		if pPlane.FaceIntersect(currentFace) {
			// If the face is split by the plane...
			split := pPlane.SplitFace(currentFace)
			if split == nil {
				continue
			}
			for _, facePart := range split {
				if facePart != nil && len(facePart.Points) > 0 {
					if pPlane.Where(facePart) <= 0 {
						fvLeft.AddFace(NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal()))
					} else {
						fvRight.AddFace(NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal()))
					}
				}
			}
		} else {
			// If the face is entirely on one side...
			w := pPlane.Where(currentFace)
			if w <= 0 {
				fvLeft.AddFace(currentFace)
			} else {
				fvRight.AddFace(currentFace)
			}
		}
	}

	// Build the left and right sub-trees from the new lists.
	if fvLeft.FaceCount() > 0 {
		parent.Left = o.createBspTree(fvLeft, newFaces, newNormMesh)
	}
	if fvRight.FaceCount() > 0 {
		parent.Right = o.createBspTree(fvRight, newFaces, newNormMesh)
	}

	return parent
}

func (o *Object3d) choosePlane(fs *FaceStore) *Face {
	leastFace, leastFaceTotal := 0, fs.FaceCount()

	for chosen := 0; chosen < fs.FaceCount(); chosen++ {
		total := 0
		p := fs.GetFace(chosen).GetPlane()
		for i := 0; i < fs.FaceCount(); i++ {
			if i == chosen {
				continue
			}
			f := fs.GetFace(i)
			if p.FaceIntersect(f) {
				total++
			}
		}
		if total < leastFaceTotal {
			leastFaceTotal = total
			leastFace = chosen
			if total == 0 {
				break
			}
		}
	}
	return fs.RemoveFaceAt(leastFace)
}

func NewCube() *Object3d {
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

	obj.Finished()
	return obj
}

// NewRectangle creates a new Object_3d in the shape of a cuboid (a 3D rectangle).
// It's centered at the origin and has the specified dimensions and color.
func NewRectangle(width, height, length float64, clr color.RGBA) *Object3d {
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

	obj.Finished()

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
func NewSubdividedRectangle(width, height, length float64, clr color.RGBA, subdivisions int) *Object3d {
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
	obj.Finished()
	return obj
}

// NewSphere creates a new Object_3d in the shape of an icosphere.
// An icosphere is a sphere made of a mesh of triangles, which is more
// uniform than a traditional UV sphere.
func NewSphere(radius float64, subdivisions int, clr color.RGBA) *Object3d {
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

	obj.Finished()
	return obj
}

// NewUVSphere creates a sphere based on latitude/longitude rings (sectors and stacks).
// This structure allows for a perfectly straight horizontal stripe.
func NewUVSphere(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Object3d {
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
	obj.Finished()
	return obj
}

func NewUVSphere2(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Object3d {
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

	obj.Finished()
	return obj
}

func LoadObjectFromDXFFile(fileName string, reverse int) (*Object3d, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, fmt.Errorf("could not open DXF file %s: %w", fileName, err)
	}
	defer file.Close()

	obj, err := NewObjectFromDXF(file, reverse)
	if err != nil {
		return nil, fmt.Errorf("error parsing DXF file %s: %w", fileName, err)
	}

	return obj, nil
}

// NewObjectFromDXF creates a new Object_3d by reading a simplified DXF file from
// the given reader. It returns a fully constructed object with its BSP tree
// already built, or an error if the file cannot be parsed.
func NewObjectFromDXF(reader io.Reader, reverse int) (*Object3d, error) {
	// 1. Create a new, empty object to populate.
	obj := NewObject_3d()

	scanner := bufio.NewScanner(reader)

	// Helper function to read the next line and parse it as a float64
	readFloatLine := func() (float64, error) {
		if !scanner.Scan() {
			if err := scanner.Err(); err != nil {
				return 0, err
			}
			return 0, io.EOF // Clean end of file
		}
		val, err := strconv.ParseFloat(strings.TrimSpace(scanner.Text()), 64)
		if err != nil {
			return 0, fmt.Errorf("could not parse float value '%s': %w", scanner.Text(), err)
		}
		return val, nil
	}

	for scanner.Scan() {
		if !strings.HasPrefix(scanner.Text(), "3DFACE") {
			continue
		}

		for i := 0; i < 3; i++ {
			if !scanner.Scan() {
				// On error, return a nil object and the error
				return nil, fmt.Errorf("unexpected end of file while parsing 3DFACE header")
			}
		}

		faceColor := color.RGBA{
			R: uint8(rand.Intn(256)),
			G: uint8(rand.Intn(256)),
			B: uint8(rand.Intn(256)),
			A: 255,
		}
		aFace := NewFace(nil, faceColor, nil)

		for c := 0; c < 4; c++ {
			x, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading X coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			y, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Y coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			z, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Z coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			aFace.AddPoint(x, y, z)
		}

		aFace.Finished(reverse)

		// 2. Add the face to the new object we created at the start.
		obj.theFaces.AddFace(aFace)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading from DXF source: %w", err)
	}

	// 3. Finalize the new object by building its BSP tree.
	obj.Finished()

	// 4. On success, return the fully populated object and a nil error.
	return obj, nil
}

func (o *Object3d) GetPosition() *Point3d {
	if o.position == nil {
		o.position = NewPoint3d(0, 0, 0)
	}
	return o.position
}

func (o *Object3d) SetPosition(x, y, z float64) {
	if o.position == nil {
		o.position = NewPoint3d(x, y, z)
	} else {
		o.position.X = x
		o.position.Y = y
		o.position.Z = z
	}
}
