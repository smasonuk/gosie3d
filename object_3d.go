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
	drawLinesOnly      bool

	// for when not using BSP
	faceIndicies   [][]int // Indices of the faces in the faceMesh
	normalIndicies []int
}

func (o *Object3d) SetDrawLinesOnly(only bool) {
	o.drawLinesOnly = only
}

func (o *Object3d) GetDrawLinesOnly() bool {
	return o.drawLinesOnly
}

func (o *Object3d) PaintObject(batcher *PolygonBatcher, x, y int, lightingChange bool, screenWidth, screenHeight float32) {
	if o.canPaintWithoutBSP {
		o.PaintWithoutBSP(batcher, x, y, screenHeight, screenWidth)

	} else {
		if o.root != nil {
			o.root.PaintWithShading(batcher, x, y, o.transFaceMesh.Points, o.transNormalMesh.Points, lightingChange, o.drawLinesOnly, screenWidth, screenHeight)
		}
	}
}

func (o *Object3d) SetDirectionVector(v *Vector3) {
	// normalize
	o.objectDirection = v.Copy()
	o.objectDirection.Normalize()

	fmt.Println("Object direction vector set to:", o.objectDirection)
}

func (o *Object3d) TranslateAllPoints(x, y, z float64) {
	if o.faceMesh == nil || o.faceMesh.Points == nil {
		return
	}

	for i := range o.faceMesh.Points.ThisMatrix {
		o.faceMesh.Points.ThisMatrix[i][0] += x
		o.faceMesh.Points.ThisMatrix[i][1] += y
		o.faceMesh.Points.ThisMatrix[i][2] += z
	}
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
		faceIndicies:       make([][]int, 0),
		normalIndicies:     make([]int, 0),
	}
}

func (o *Object3d) Clone() *Object3d {

	clone := &Object3d{
		// shared
		faceMesh:       o.faceMesh,
		normalMesh:     o.normalMesh,
		theFaces:       o.theFaces,
		root:           o.root,
		faceIndicies:   o.faceIndicies,
		normalIndicies: o.normalIndicies,

		// instance-specific
		transFaceMesh:      o.transFaceMesh.Copy(),
		transNormalMesh:    o.transNormalMesh.Copy(),
		rotMatrix:          IdentMatrix(),
		canPaintWithoutBSP: o.canPaintWithoutBSP,
	}
	return clone
}

func (o *Object3d) createFaceList(faces *FaceStore, newFaces *FaceMesh, newNormMesh *NormalMesh) {
	for i := 0; i < faces.FaceCount(); i++ {
		originalFace := faces.GetFace(i)
		newFace, ind := newFaces.AddFace(originalFace)
		newFaces.AddFace(newFace)

		normal := originalFace.GetNormal()
		_, normalIndex := newNormMesh.AddNormal(normal)

		o.faceIndicies = append(o.faceIndicies, ind)
		o.normalIndicies = append(o.normalIndicies, normalIndex)
	}
}

func (o *Object3d) Finished(centerObject bool, useBspTree bool) {
	log.Println("Creating BSP Tree...")
	if o.theFaces.FaceCount() > 0 {
		if useBspTree {
			o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
			o.canPaintWithoutBSP = false
		} else {
			o.createFaceList(o.theFaces, o.transFaceMesh, o.transNormalMesh)
			o.canPaintWithoutBSP = true
		}
	}

	log.Println("BSP Tree Created.")

	o.faceMesh = o.transFaceMesh.Copy()
	o.normalMesh = o.transNormalMesh.Copy()
	// &NormalMesh{Mesh: *o.transNormalMesh.Mesh.Copy()}

	if centerObject {
		o.CentreObject()
	}

	log.Printf("Points: %d", len(o.faceMesh.Points.ThisMatrix))
	log.Printf("Normals: %d", len(o.normalMesh.Points.ThisMatrix))

	o.CalcSize()
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

func (o *Object3d) GetExtents() (float64, float64, float64) {
	return o.xLength, o.yLength, o.zLength
}

func (o *Object3d) XLength() float64 {
	return o.xLength
}

func (o *Object3d) YLength() float64 {
	return o.yLength
}

func (o *Object3d) ZLength() float64 {
	return o.zLength
}

func (o *Object3d) CalcSize() {
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
	obj.Finished(false, true)

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

func (o *Object3d) RotateY(amountOfMovementInRads float64) {
	o.rotMatrix = ApplyRotateYObjectMatrix(amountOfMovementInRads, o.rotMatrix)
}

func ApplyRollObjectMatrix(directionOfRoll float64, amountOfMovement float64, existingMatrix *Matrix) *Matrix {
	rotMatrixY := NewRotationMatrix(ROTY, -directionOfRoll)
	rotMatrixYBack := NewRotationMatrix(ROTY, directionOfRoll)
	rotMatrixX := NewRotationMatrix(ROTX, -amountOfMovement)
	all := rotMatrixY.MultiplyBy(rotMatrixX).MultiplyBy(rotMatrixYBack)

	return all.MultiplyBy(existingMatrix)
}

func ApplyRotateYObjectMatrix(amountOfMovementInRads float64, existingMatrix *Matrix) *Matrix {
	rotMatrix := NewRotationMatrix(ROTY, amountOfMovementInRads)
	return rotMatrix.MultiplyBy(existingMatrix)
}

func (o *Object3d) AddFacesFromObject(other *Object3d) {
	if other == nil || other.theFaces == nil || other.theFaces.FaceCount() == 0 {
		return
	}

	for i := 0; i < other.theFaces.FaceCount(); i++ {
		face := other.theFaces.GetFace(i)
		newFace := face.Copy()

		o.theFaces.AddFace(newFace)
	}
}

func (o *Object3d) PaintWithoutBSP(batcher *PolygonBatcher, x, y int, screenHeight, screenWidth float32) {
	points := make([][]float64, 0, 10)
	for i := 0; i < len(o.faceIndicies); i++ {
		faceIndices := o.faceIndicies[i]
		normalIndex := o.normalIndicies[i]
		facePointsInCameraSpace := points[:]
		for _, index := range faceIndices {
			point := o.transFaceMesh.Points.ThisMatrix[index]
			facePointsInCameraSpace = append(facePointsInCameraSpace, point)
		}
		// face := o.transFaceMesh.faces[i]
		face := o.theFaces.faces[i]

		normal := o.transNormalMesh.Points.ThisMatrix[normalIndex]

		o.paintFace(batcher, x, y, facePointsInCameraSpace, normal, screenWidth, screenHeight, face)
	}
}

func (o *Object3d) paintFace(batcher *PolygonBatcher, x, y int, points [][]float64, normal []float64, screenWidth, screenHeight float32, face *Face) {

	firstPoint := points[0]
	where := normal[0]*firstPoint[0] +
		normal[1]*firstPoint[1] +
		normal[2]*firstPoint[2]

	if where > 0 { // Facing the camera
		o.paintFace2(batcher,
			x,
			y,
			firstPoint,
			points,
			normal,
			screenWidth,
			screenHeight,
			face,
		)
	}

}

func (o *Object3d) paintFace2(
	batcher *PolygonBatcher,
	x, y int,
	firstTransformedPoint []float64,
	initial3DPoints [][]float64,
	transformedNormal []float64,
	screenWidth, screenHeight float32,
	face *Face,
) bool {

	pointsToUse := clipPolygonAgainstNearPlane(initial3DPoints)
	if len(pointsToUse) < 3 {
		return false
	}

	initialScreenPoints := make([]Point, len(pointsToUse))
	for i, point := range pointsToUse {
		// At this stage, point[2] (z) is guaranteed to be >= nearPlaneZ,
		// so perspective division is safe.
		z := float32(point[2])
		initialScreenPoints[i] = Point{
			X: float32((conversionFactor*point[0])/float64(z)) + float32(x),
			Y: float32((conversionFactor*point[1])/float64(z)) + float32(y),
		}
	}

	// clip to screen
	clippedPoints := clipPolygon(initialScreenPoints, screenWidth, screenHeight)
	if len(clippedPoints) < 3 {
		return false
	}

	finalScreenPointsX := make([]float32, len(clippedPoints))
	finalScreenPointsY := make([]float32, len(clippedPoints))
	for i, p := range clippedPoints {
		finalScreenPointsX[i] = p.X
		finalScreenPointsY[i] = p.Y
	}

	col := face.Col
	polyColor := col
	if true {
		shadingRefPoint := firstTransformedPoint
		polyColor = getColor(shadingRefPoint, transformedNormal, polyColor)
	}

	// if !linesOnly {
	black := color.RGBA{R: 100, G: 100, B: 100, A: 25}
	batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, polyColor, black, 1.0)

	// } else {
	// 	black := color.RGBA{R: 0, G: 0, B: 0, A: 255}
	// 	// batcher.AddPolygonOutline(finalScreenPointsX, finalScreenPointsY, 1, polyColor)
	// 	batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, black, polyColor, 1.0)
	// }

	return false
}

// GetColor calculates the color of a polygon based on simple lighting.
func getColor(
	firstTransformedPoint []float64,
	transformedNormal []float64,
	polyColor color.RGBA,
) color.RGBA {
	const ambientLight = 0.65
	const spotlightConePower = 10.0
	const spotlightLightAmount = 1.0 - ambientLight

	diffuseFactor := transformedNormal[2]
	if diffuseFactor < 0 {
		diffuseFactor = 0
	}

	var spotlightFactor float64
	lenVecToPoint := GetLength2(firstTransformedPoint)

	if lenVecToPoint > 0 {
		cosAngle := firstTransformedPoint[2] / lenVecToPoint
		if cosAngle < 0 {
			cosAngle = 0
		}
		spotlightFactor = math.Pow(cosAngle, spotlightConePower)
	} else {
		spotlightFactor = 1.0
	}

	spotlightBrightness := diffuseFactor * spotlightFactor * spotlightLightAmount
	finalBrightness := ambientLight + spotlightBrightness

	c := 240 - int(finalBrightness*240)
	min := 7
	r1 := clamp2(int(polyColor.R)-c, min, 255)
	g1 := clamp2(int(polyColor.G)-c, min, 255)
	b1 := clamp2(int(polyColor.B)-c, min, 255)
	polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: polyColor.A}

	return polyColor
}

// func (o *Object3d) drae(batcher *PolygonBatcher, x, y int, transPoints *Matrix, transNormals *Matrix, doShading bool, linesOnly bool, screenWidth, screenHeight float32) {

// 	transformedNormal := transNormals.ThisMatrix[b.normalIndex]
// 	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

// 	// Determine if the polygon is facing the camera
// 	where := transformedNormal[0]*firstTransformedPoint[0] +
// 		transformedNormal[1]*firstTransformedPoint[1] +
// 		transformedNormal[2]*firstTransformedPoint[2]

// 	if where <= 0 { // Facing away from the camera so don't paint it

// 	} else {

// 		o.paintPoly(batcher, x, y,
// 			transPoints,
// 			transNormals,
// 			doShading,
// 			firstTransformedPoint,
// 			transformedNormal,
// 			linesOnly,
// 			screenWidth,
// 			screenHeight,
// 		)

// 	}
// }

// func (o *Object3d) paintPoly(
// 	batcher *PolygonBatcher,
// 	x, y int,
// 	verticesInCameraSpace *Matrix,
// 	normalsInCameraSpace *Matrix,
// 	shadePoly bool,
// 	firstTransformedPoint []float64,
// 	transformedNormal []float64,
// 	linesOnly bool,
// 	screenWidth, screenHeight float32,
// ) bool {

// 	initial3DPoints := b.pointsToUse[:0] // Reuse the slice to avoid allocation
// 	for _, pointIndex := range b.facePointIndices {
// 		initial3DPoints = append(initial3DPoints, verticesInCameraSpace.ThisMatrix[pointIndex])
// 	}

// 	pointsToUse := clipPolygonAgainstNearPlane(initial3DPoints)

// 	// If clipping results in a polygon with too few vertices, don't draw it.
// 	if len(pointsToUse) < 3 {
// 		return false
// 	}

// 	initialScreenPoints := make([]Point, len(pointsToUse))
// 	for i, point := range pointsToUse {
// 		// At this stage, point[2] (z) is guaranteed to be >= nearPlaneZ,
// 		// so perspective division is safe.
// 		z := float32(point[2])
// 		initialScreenPoints[i] = Point{
// 			X: float32((conversionFactor*point[0])/float64(z)) + float32(x),
// 			Y: float32((conversionFactor*point[1])/float64(z)) + float32(y),
// 		}
// 	}

// 	clippedPoints := clipPolygon(initialScreenPoints, screenWidth, screenHeight)

// 	if len(clippedPoints) < 3 {
// 		return false
// 	}

// 	finalScreenPointsX := make([]float32, len(clippedPoints))
// 	finalScreenPointsY := make([]float32, len(clippedPoints))
// 	for i, p := range clippedPoints {
// 		finalScreenPointsX[i] = p.X
// 		finalScreenPointsY[i] = p.Y
// 	}

// 	polyColor := color.RGBA{R: b.colRed, G: b.colGreen, B: b.colBlue, A: b.colAlpha}
// 	if shadePoly {
// 		shadingRefPoint := verticesInCameraSpace.ThisMatrix[b.facePointIndices[0]]
// 		polyColor = b.GetColor(shadingRefPoint, transformedNormal, polyColor)
// 	}

// 	if !linesOnly {
// 		black := color.RGBA{R: 100, G: 100, B: 100, A: 25}
// 		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, polyColor, black, 1.0)

// 	} else {
// 		black := color.RGBA{R: 0, G: 0, B: 0, A: 255}
// 		// batcher.AddPolygonOutline(finalScreenPointsX, finalScreenPointsY, 1, polyColor)
// 		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, black, polyColor, 1.0)
// 	}

// 	return false
// }
