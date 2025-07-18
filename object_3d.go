package gosie3d

import (
	"bufio"
	"fmt"
	"image/color"
	"io"
	"log"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"github.com/hajimehoshi/ebiten/v2"
)

// const (
// 	PAINT_MODE_BSP PaintMode = iota
// 	PAINT_MODE_IGNORE_BACKFACE
// )

// type PaintMode int

type Object3d struct {
	faceMesh           *FaceMesh
	normalMesh         *NormalMesh
	transFaceMesh      *FaceMesh
	transNormalMesh    *NormalMesh
	theFaces           *FaceStore
	root               *BspNode
	rotMatrix          *Matrix
	position           *Point3d
	canPaintWithoutBSP bool
	xLength            float64
	yLength            float64
	zLength            float64
	objectDirection    *Vector3
}

func (o *Object3d) PaintObject(screen *ebiten.Image, x, y int, lightingChange bool) {
	if o.root != nil {
		o.root.PaintWithShading(screen, x, y, o.transFaceMesh.Points, o.transNormalMesh.Points, lightingChange, true)
	}
}

func (o *Object3d) SetDirectionVector(v *Vector3) {
	// normalize
	o.objectDirection = v.Copy()
	o.objectDirection.Normalize()

	fmt.Println("Object direction vector set to:", o.objectDirection)

}

// apply matrix to the direction vector to transform it.
func (o *Object3d) ApplyDirectionVector(m *Matrix) *Vector3 {
	if o.objectDirection == nil {
		return nil
	}

	// create matrix from the direction vector
	vecMatrix := NewMatrix()
	vecMatrix.AddRow([]float64{o.objectDirection.X, o.objectDirection.Y, o.objectDirection.Z, 1.0})
	// apply the matrix transformation
	transformed := m.MultiplyBy(vecMatrix)

	// return the transformed vector
	return NewVector3(
		transformed.ThisMatrix[0][0],
		transformed.ThisMatrix[0][1],
		transformed.ThisMatrix[0][2],
	)
}

// apply matrix to the direction vector to transform it.
func ApplyDirectionVector(v *Vector3, m *Matrix) *Vector3 {

	// create matrix from the direction vector
	vecMatrix := NewMatrix()
	vecMatrix.AddRow([]float64{v.X, v.Y, v.Z, 1.0})
	// apply the matrix transformation
	transformed := m.MultiplyBy(vecMatrix)

	// return the transformed vector
	return NewVector3(
		transformed.ThisMatrix[0][0],
		transformed.ThisMatrix[0][1],
		transformed.ThisMatrix[0][2],
	)
}

func (o *Object3d) SetRotMatrix(m *Matrix) {
	o.rotMatrix = m
}

func NewObject_3d(canPaintWithoutBSP bool) *Object3d {
	return &Object3d{
		transFaceMesh:      NewFaceMesh(),
		transNormalMesh:    NewNormalMesh(),
		theFaces:           NewFaceStore(),
		rotMatrix:          IdentMatrix(),
		canPaintWithoutBSP: canPaintWithoutBSP,
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

func (o *Object3d) Finished(centerObject bool) {

	log.Println("Creating BSP Tree...")
	if o.theFaces.FaceCount() > 0 {
		// o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
		o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
	}

	log.Println("BSP Tree Created.")

	o.faceMesh = &FaceMesh{Mesh: *o.transFaceMesh.Mesh.Copy()}
	o.normalMesh = &NormalMesh{Mesh: *o.transNormalMesh.Mesh.Copy()}

	if centerObject {
		o.CentreObject()
	}

	log.Printf("Points: %d", len(o.faceMesh.Points.ThisMatrix))
	log.Printf("Normals: %d", len(o.normalMesh.Points.ThisMatrix))

	o.calcSize()
}

// Moves all points so that the 0,0,0 is the center of the object.
func (o *Object3d) CentreObject() {
	if o.faceMesh == nil || len(o.faceMesh.Points.ThisMatrix) == 0 {
		return
	}

	// Calculate the bounding box of the object.
	minX, maxX := o.faceMesh.Points.ThisMatrix[0][0], o.faceMesh.Points.ThisMatrix[0][0]
	minY, maxY := o.faceMesh.Points.ThisMatrix[0][1], o.faceMesh.Points.ThisMatrix[0][1]
	minZ, maxZ := o.faceMesh.Points.ThisMatrix[0][2], o.faceMesh.Points.ThisMatrix[0][2]
	for _, point := range o.faceMesh.Points.ThisMatrix {
		if point[0] < minX {
			minX = point[0]
		} else if point[0] > maxX {
			maxX = point[0]
		}
		if point[1] < minY {
			minY = point[1]
		} else if point[1] > maxY {
			maxY = point[1]
		}

		if point[2] < minZ {
			minZ = point[2]
		} else if point[2] > maxZ {
			maxZ = point[2]
		}
	}

	//calculate the center of the bounding box.
	centerX := (minX + maxX) / 2.0
	centerY := (minY + maxY) / 2.0
	centerZ := (minZ + maxZ) / 2.0

	// Move all points so that the center of the bounding box is at 0,0,0.
	for i := range o.faceMesh.Points.ThisMatrix {
		o.faceMesh.Points.ThisMatrix[i][0] -= centerX
		o.faceMesh.Points.ThisMatrix[i][1] -= centerY
		o.faceMesh.Points.ThisMatrix[i][2] -= centerZ
	}
}

func (o *Object3d) ZLength() float64 {
	return o.zLength
}

func (o *Object3d) YLength() float64 {
	return o.yLength
}

func (o *Object3d) XLength() float64 {
	return o.xLength
}

func (o *Object3d) calcSize() {
	if o.faceMesh == nil || len(o.faceMesh.Points.ThisMatrix) == 0 {
		o.xLength = 0
		o.yLength = 0
		o.zLength = 0
		return
	}

	minX, maxX := o.faceMesh.Points.ThisMatrix[0][0], o.faceMesh.Points.ThisMatrix[0][0]
	minY, maxY := o.faceMesh.Points.ThisMatrix[0][1], o.faceMesh.Points.ThisMatrix[0][1]
	minZ, maxZ := o.faceMesh.Points.ThisMatrix[0][2], o.faceMesh.Points.ThisMatrix[0][2]
	for _, point := range o.faceMesh.Points.ThisMatrix {
		if point[0] < minX {
			minX = point[0]
		} else if point[0] > maxX {
			maxX = point[0]
		}
		if point[1] < minY {
			minY = point[1]
		} else if point[1] > maxY {
			maxY = point[1]
		}

		if point[2] < minZ {
			minZ = point[2]
		} else if point[2] > maxZ {
			maxZ = point[2]
		}
	}
	o.xLength = maxX - minX
	o.yLength = maxY - minY
	o.zLength = maxZ - minZ
	log.Printf("Object size: X: %.2f, Y: %.2f, Z: %.2f", o.xLength, o.yLength, o.zLength)
}

func (o *Object3d) ApplyMatrix(m *Matrix) {
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
	obj := NewObject_3d(false)

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
	obj.Finished(false)

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

func (o *Object3d) RollObject(directionOfRoll float64, amountOfMovement float64) {
	o.rotMatrix = ApplyRollObjectMatrix(directionOfRoll, amountOfMovement, o.rotMatrix)
}

func ApplyRollObjectMatrix(directionOfRoll float64, amountOfMovement float64, existingMatrix *Matrix) *Matrix {
	rotMatrixY := NewRotationMatrix(ROTY, -directionOfRoll)
	rotMatrixYBack := NewRotationMatrix(ROTY, directionOfRoll)
	rotMatrixX := NewRotationMatrix(ROTX, -amountOfMovement)
	all := rotMatrixY.MultiplyBy(rotMatrixX).MultiplyBy(rotMatrixYBack)

	return all.MultiplyBy(existingMatrix)
}

func (o *Object3d) PaintObjectWithBackfaceCullingOnly() {

}
