package gosie3d

import (
	"bufio"
	"fmt"
	"image/color"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

type Model struct {
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
	faceIndicies     [][]int // Indices of the faces in the faceMesh
	normalIndicies   []int
	drawAllFaces     bool // If true, draw all faces regardless of visibility
	dontDrawOutlines bool // If true, don't draw outlines of polygons
}

func (o *Model) SetDontDrawOutlines(dontDraw bool) {
	o.dontDrawOutlines = dontDraw
}

func (o *Model) GetDontDrawOutlines() bool {
	return o.dontDrawOutlines
}

func (o *Model) SetDrawAllFaces(draw bool) {
	o.drawAllFaces = draw
}

func (o *Model) SetDrawLinesOnly(only bool) {
	o.drawLinesOnly = only
}

func (o *Model) GetDrawLinesOnly() bool {
	return o.drawLinesOnly
}

func (o *Model) PaintObject(batcher *PolygonBatcher, x, y int, lightingChange bool, screenWidth, screenHeight float32) {
	if o.canPaintWithoutBSP {
		o.paintWithoutBSP(batcher, x, y, screenHeight, screenWidth)

	} else {
		if o.root != nil {
			o.root.PaintWithShading(batcher, x, y, o.transFaceMesh.Points, o.transNormalMesh.Points, lightingChange, o.drawLinesOnly, screenWidth, screenHeight,
				o.dontDrawOutlines)
		}
	}
}

func (o *Model) SetDirectionVector(v *Vector3) {
	// normalize
	o.objectDirection = v.Copy()
	o.objectDirection.Normalize()

	fmt.Println("Object direction vector set to:", o.objectDirection)
}

func (o *Model) TranslateAllPoints(x, y, z float64) {
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
func (o *Model) ApplyDirectionVector(m *Matrix) *Vector3 {
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

func (o *Model) SetRotMatrix(m *Matrix) {
	o.rotMatrix = m
}

func (o *Model) GetRotMatrix() *Matrix {
	if o.rotMatrix == nil {
		o.rotMatrix = IdentMatrix()
	}
	return o.rotMatrix
}

func NewObject_3d() *Model {
	return &Model{
		transFaceMesh:   NewFaceMesh(),
		transNormalMesh: NewNormalMesh(),
		theFaces:        NewFaceStore(),
		rotMatrix:       IdentMatrix(),
		faceIndicies:    make([][]int, 0),
		normalIndicies:  make([]int, 0),
	}
}

func (o *Model) Clone() *Model {

	clone := &Model{
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

func (o *Model) CreateFaceList() {
	faces, newFaces, newNormMesh := o.theFaces, o.transFaceMesh, o.transNormalMesh
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

func (o *Model) Finished(centerObject bool, useBspTree bool) {
	log.Println("Creating BSP Tree...")
	if o.theFaces.FaceCount() > 0 {
		if useBspTree {
			o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
			o.canPaintWithoutBSP = false
		} else {
			o.CreateFaceList()
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
func (o *Model) CentreObject() {
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

func (o *Model) GetExtents() (float64, float64, float64) {
	return o.xLength, o.yLength, o.zLength
}

func (o *Model) XLength() float64 {
	return o.xLength
}

func (o *Model) YLength() float64 {
	return o.yLength
}

func (o *Model) ZLength() float64 {
	return o.zLength
}

func (o *Model) CalcSize() {
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

func (o *Model) ApplyMatrix(m *Matrix) {
	o.rotMatrix = m.MultiplyBy(o.rotMatrix)
}

func (o *Model) ApplyMatrixTemp(aMatrix *Matrix) {
	rotMatrixTemp := aMatrix.MultiplyBy(o.rotMatrix)

	// Use the new, correct method to transform the normals (rotation only).
	rotMatrixTemp.TransformNormals(o.normalMesh.Points, o.transNormalMesh.Points)

	// Use the original method to transform the vertex positions (rotation and translation).
	rotMatrixTemp.TransformObj(o.faceMesh.Points, o.transFaceMesh.Points)
}

func (o *Model) ApplyMatrixPermanent(aMatrix *Matrix) {
	o.ApplyMatrixTemp(aMatrix)

	o.normalMesh.Points = o.transNormalMesh.Points.Copy()
	o.faceMesh.Points = o.transFaceMesh.Points.Copy()
}

func (o *Model) ApplyObjMatrixPermanent() {
	o.ApplyMatrixTemp(o.rotMatrix)

	o.normalMesh.Points = o.transNormalMesh.Points.Copy()
	o.faceMesh.Points = o.transFaceMesh.Points.Copy()

	o.rotMatrix = IdentMatrix()
}

func (o *Model) createBspTree(faces *FaceStore, newFaces *FaceMesh, newNormMesh *NormalMesh) *BspNode {
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
						f1 := NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal())
						fvLeft.AddFace(f1)
					} else {
						f2 := NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal())
						fvRight.AddFace(f2)
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

func (o *Model) choosePlane(fs *FaceStore) *Face {
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

func LoadObjectFromDXFFile(fileName string, reverse int) (*Model, error) {
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
func NewObjectFromDXF(reader io.Reader, reverse int) (*Model, error) {
	// create a new, empty object to populate.
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
			R: 100,
			G: 10,
			B: 58,
			A: 255,
		}
		aFace := NewFace(nil, faceColor, nil)

		for c := 0; c < 4; c++ {
			x, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading X coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan()

			y, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Y coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan()

			z, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Z coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan()

			aFace.AddPoint(x, y, z)
		}

		aFace.Finished(reverse)

		//Add the face to the new object we created at the start.
		obj.theFaces.AddFace(aFace)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading from DXF source: %w", err)
	}

	// Finalize the new object by building its BSP tree.
	obj.Finished(false, true)

	// On success, return the fully populated object and a nil error.
	return obj, nil
}

func (o *Model) GetPosition() *Point3d {
	if o.position == nil {
		o.position = NewPoint3d(0, 0, 0)
	}
	return o.position
}

func (o *Model) SetPosition(x, y, z float64) {
	if o.position == nil {
		o.position = NewPoint3d(x, y, z)
	} else {
		o.position.X = x
		o.position.Y = y
		o.position.Z = z
	}
}

func (o *Model) RollObject(directionOfRoll float64, amountOfMovement float64) {
	o.rotMatrix = ApplyRollObjectMatrix(directionOfRoll, amountOfMovement, o.rotMatrix)
}

func (o *Model) RotateY(amountOfMovementInRads float64) {
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

func (o *Model) AddFacesFromObject(other *Model) {
	if other == nil || other.theFaces == nil || other.theFaces.FaceCount() == 0 {
		return
	}

	for i := 0; i < other.theFaces.FaceCount(); i++ {
		face := other.theFaces.GetFace(i)
		newFace := face.Copy()

		o.theFaces.AddFace(newFace)
	}
}

// get face indices
func (o *Model) GetFaceIndices() [][]int {
	return o.faceIndicies
}

// get transFaceMesh
func (o *Model) GetTransFaceMesh() *FaceMesh {
	return o.transFaceMesh
}

func (o *Model) paintWithoutBSP(batcher *PolygonBatcher, x, y int, screenHeight, screenWidth float32) {
	points := make([][]float64, 0, 10)
	for i := 0; i < len(o.faceIndicies); i++ {
		faceIndices := o.faceIndicies[i]
		normalIndex := o.normalIndicies[i]
		facePointsInCameraSpace := points[:]
		for _, index := range faceIndices {
			point := o.transFaceMesh.Points.ThisMatrix[index]
			facePointsInCameraSpace = append(facePointsInCameraSpace, point)
		}
		face := o.theFaces.faces[i]

		normal := o.transNormalMesh.Points.ThisMatrix[normalIndex]

		o.paintFace(batcher, x, y, facePointsInCameraSpace, normal, screenWidth, screenHeight, face)
	}
}

// func (o *Model) FacesIntersectingLine(startLine, endLine *Point3d) []*Face {
// 	faces := make([]*Face, 0, 5)
// 	points := make([][]float64, 0, 10)
// 	for i := 0; i < len(o.faceIndicies); i++ {
// 		faceIndices := o.faceIndicies[i]
// 		// normalIndex := o.normalIndicies[i]
// 		facePointsInCameraSpace := points[:]
// 		for _, index := range faceIndices {
// 			point := o.transFaceMesh.Points.ThisMatrix[index]
// 			facePointsInCameraSpace = append(facePointsInCameraSpace, point)
// 		}

// 		intersects := LineIntersectsPolygon(startLine, endLine, facePointsInCameraSpace)
// 		if intersects {
// 			// If the line intersects the polygon, add the face to the list.
// 			face := o.theFaces.GetFace(i)
// 			faces = append(faces, face)
// 		}
// 	}

// 	if len(faces) == 0 {
// 		return nil
// 	}

// 	return faces
// }

func (o *Model) BspNodesIntersectingLine(startLine, endLine *Point3d) *BspNode {
	nodes := make([]*BspNode, 0, 5)
	points := make([][]float64, 0, 10)
	// for i := 0; i < len(o.faceIndicies); i++ {
	// 	faceIndices := o.faceIndicies[i]
	// 	// normalIndex := o.normalIndicies[i]
	// 	facePointsInCameraSpace := points[:]
	// 	for _, index := range faceIndices {
	// 		point := o.transFaceMesh.Points.ThisMatrix[index]
	// 		facePointsInCameraSpace = append(facePointsInCameraSpace, point)
	// 	}

	// 	intersects := LineIntersectsPolygon(startLine, endLine, facePointsInCameraSpace)
	// 	if intersects {
	// 		// If the line intersects the polygon, add the face to the list.
	// 		face := o.theFaces.GetFace(i)
	// 		faces = append(faces, face)
	// 	}
	// }

	var traverseNodes func(node *BspNode)
	traverseNodes = func(node *BspNode) {
		if node == nil {
			return
		}
		// Each node in the BSP tree represents a polygon.
		// allFaceIndices = append(allFaceIndices, node.facePointIndices)
		// Recursively traverse the child nodes.

		// for i := 0; i < len(node.facePointIndices); i++ {
		// 	faceIndices := node.facePointIndices[i]
		// normalIndex := o.normalIndicies[i]

		facePointsInCameraSpace := points[:]
		for _, index := range node.facePointIndices {
			point := o.transFaceMesh.Points.ThisMatrix[index]
			facePointsInCameraSpace = append(facePointsInCameraSpace, point)
		}
		intersects := LineIntersectsPolygon(startLine, endLine, facePointsInCameraSpace)
		if intersects {
			// If the line intersects the polygon, add the face to the list.
			// face := o.theFaces.GetFace(i)
			nodes = append(nodes, node)
		}

		// }
		traverseNodes(node.Left)
		traverseNodes(node.Right)
	}

	traverseNodes(o.root)

	if len(nodes) == 0 {
		return nil
	}

	// loop through the nodes and find the one closest to startLine point
	closestNode := nodes[0]
	distance := math.MaxFloat64
	for _, node := range nodes {
		nodex, nodey, nodez := node.GetAveragePoint(o.transFaceMesh.Points)
		dx := nodex - startLine.X
		dy := nodey - startLine.Y
		dz := nodez - startLine.Z
		d := math.Sqrt(dx*dx + dy*dy + dz*dz)
		if d < distance {
			distance = d
			closestNode = node
		}
	}

	return closestNode
}

func (o *Model) paintFace(batcher *PolygonBatcher, x, y int, points [][]float64, normal []float64, screenWidth, screenHeight float32, face *Face) {

	firstPoint := points[0]
	where := 1.0
	if !o.drawAllFaces {
		where = normal[0]*firstPoint[0] +
			normal[1]*firstPoint[1] +
			normal[2]*firstPoint[2]
	}

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

func (o *Model) paintFace2(
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
		// cf := float64(screenWidth)
		// At this stage, point[2] (z) is guaranteed to be >= nearPlaneZ,
		// so perspective division is safe.
		z := float32(point[2])
		initialScreenPoints[i] = Point{
			X: ConvertToScreenX(float64(screenWidth), float64(screenHeight), point[0], float64(z)),
			Y: ConvertToScreenY(float64(screenWidth), float64(screenHeight), point[1], float64(z)),
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

	if !o.drawLinesOnly {
		black := color.RGBA{R: 50, G: 50, B: 50, A: 25}
		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, polyColor, black, 1.0)

	} else {
		black := color.RGBA{R: 0, G: 0, B: 0, A: 255}
		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, black, polyColor, 1.0)
	}

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
func (o *Model) ScaleAllPoints(scale float64) {
	if o.faceMesh == nil || o.faceMesh.Points == nil {
		return
	}

	for i := range o.faceMesh.Points.ThisMatrix {
		o.faceMesh.Points.ThisMatrix[i][0] *= scale
		o.faceMesh.Points.ThisMatrix[i][1] *= scale
		o.faceMesh.Points.ThisMatrix[i][2] *= scale
	}
}

// Represents a vertex with its color.
type Vertex struct {
	X, Y, Z float64
	Color   color.RGBA
}

func LoadObjectFromPLYFile(fileName string, reverse int, useBsp bool) (*Model, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, fmt.Errorf("could not open PLY file %s: %w", fileName, err)
	}
	defer file.Close()

	obj, err := LoadObjectFromPLYReader(file, reverse, useBsp)
	if err != nil {
		return nil, fmt.Errorf("error parsing PLY file %s: %w", fileName, err)
	}

	return obj, nil
}

// func LoadObjectFromPLYReader(reader io.Reader, reverse int) (*Model, error) {
// 	obj := NewObject_3d()
// 	scanner := bufio.NewScanner(reader)

// 	var vertexCount, faceCount int

// 	// 1. Parse the header first
// 	for scanner.Scan() {
// 		line := scanner.Text()
// 		parts := strings.Fields(line)
// 		if len(parts) == 3 && parts[0] == "element" {
// 			if parts[1] == "vertex" {
// 				vertexCount, _ = strconv.Atoi(parts[2])
// 			} else if parts[1] == "face" {
// 				faceCount, _ = strconv.Atoi(parts[2])
// 			}
// 		}
// 		if line == "end_header" {
// 			break
// 		}
// 	}

// 	// 2. Read all vertices and their colors
// 	vertices := make([]Vertex, 0, vertexCount)
// 	for i := 0; i < vertexCount; i++ {
// 		scanner.Scan()
// 		parts := strings.Fields(scanner.Text())
// 		if len(parts) < 6 { // expecting x, y, z, r, g, b
// 			return nil, fmt.Errorf("invalid vertex data on line %d", i)
// 		}

// 		x, _ := strconv.ParseFloat(parts[0], 64)
// 		y, _ := strconv.ParseFloat(parts[1], 64)
// 		z, _ := strconv.ParseFloat(parts[2], 64)
// 		r, _ := strconv.ParseUint(parts[3], 10, 8)
// 		g, _ := strconv.ParseUint(parts[4], 10, 8)
// 		b, _ := strconv.ParseUint(parts[5], 10, 8)

// 		vertices = append(vertices, Vertex{
// 			X: x, Y: y, Z: z,
// 			Color: color.RGBA{R: uint8(r), G: uint8(g), B: uint8(b), A: 255},
// 		})
// 	}

// 	// 3. Read the face definitions and create faces using the vertices
// 	for i := 0; i < faceCount; i++ {
// 		scanner.Scan()
// 		parts := strings.Fields(scanner.Text())

// 		// parts[0] is the number of vertices in this face
// 		numFaceVerts, _ := strconv.Atoi(parts[0])
// 		if len(parts) != numFaceVerts+1 {
// 			return nil, fmt.Errorf("invalid face data on line %d", i)
// 		}

// 		// We can average the vertex colors to get a face color.
// 		var r, g, b uint32

// 		aFace := NewFace(nil, color.RGBA{44, 2, 65, 255}, nil)

// 		for j := 1; j <= numFaceVerts; j++ {
// 			idx, _ := strconv.Atoi(parts[j])
// 			vert := vertices[idx]
// 			aFace.AddPoint(vert.X, vert.Y, vert.Z)
// 			r += uint32(vert.Color.R)
// 			g += uint32(vert.Color.G)
// 			b += uint32(vert.Color.B)
// 		}

// 		// Set the average color for the face
// 		avgColor := color.RGBA{
// 			R: uint8(r / uint32(numFaceVerts)),
// 			G: uint8(g / uint32(numFaceVerts)),
// 			B: uint8(b / uint32(numFaceVerts)),
// 			A: 255,
// 		}
// 		aFace.SetColor(avgColor) // Assuming you have a method to set the color

// 		aFace.Finished(reverse)
// 		obj.theFaces.AddFace(aFace)
// 	}

// 	if err := scanner.Err(); err != nil {
// 		return nil, fmt.Errorf("error reading from PLY source: %w", err)
// 	}

// 	obj.Finished(false, true)
// 	return obj, nil
// }

func LoadObjectFromPLYReader(reader io.Reader, reverse int, useBsp bool) (*Model, error) {
	obj := NewObject_3d()
	scanner := bufio.NewScanner(reader)

	var vertexCount, faceCount int
	var hasVertexColor, hasFaceColor bool
	var currentElement string

	// 1. Parse the header intelligently to detect color format
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Fields(line)
		if len(parts) == 0 {
			continue
		}

		switch parts[0] {
		case "element":
			if len(parts) == 3 {
				currentElement = parts[1] // Track if we are defining a 'vertex' or 'face'
				if parts[1] == "vertex" {
					vertexCount, _ = strconv.Atoi(parts[2])
				} else if parts[1] == "face" {
					faceCount, _ = strconv.Atoi(parts[2])
				}
			}
		case "property":
			// Check for color properties within the current element definition
			if len(parts) > 2 && (parts[2] == "red" || parts[2] == "diffuse_red") {
				if currentElement == "vertex" {
					hasVertexColor = true
				} else if currentElement == "face" {
					hasFaceColor = true
				}
			}
		case "end_header":
			goto endHeaderLoop
		}
	}
endHeaderLoop:

	// Read all vertices
	vertices := make([]Vertex, 0, vertexCount)
	for i := 0; i < vertexCount; i++ {
		if !scanner.Scan() {
			return nil, fmt.Errorf("unexpected end of file while reading vertices")
		}
		parts := strings.Fields(scanner.Text())

		// Conditional Vertex Parsing
		if hasVertexColor {
			if len(parts) < 6 {
				return nil, fmt.Errorf("invalid vertex-color data on line %d", i)
			}
			x, _ := strconv.ParseFloat(parts[0], 64)
			y, _ := strconv.ParseFloat(parts[1], 64)
			z, _ := strconv.ParseFloat(parts[2], 64)
			r, _ := strconv.ParseUint(parts[3], 10, 8)
			g, _ := strconv.ParseUint(parts[4], 10, 8)
			b, _ := strconv.ParseUint(parts[5], 10, 8)
			vertices = append(vertices, Vertex{X: x, Y: y, Z: z, Color: color.RGBA{R: uint8(r), G: uint8(g), B: uint8(b), A: 255}})
		} else { // No vertex color
			if len(parts) < 3 {
				return nil, fmt.Errorf("invalid vertex data on line %d", i)
			}
			x, _ := strconv.ParseFloat(parts[0], 64)
			y, _ := strconv.ParseFloat(parts[1], 64)
			z, _ := strconv.ParseFloat(parts[2], 64)
			// Use a default white color if none is provided
			vertices = append(vertices, Vertex{X: x, Y: y, Z: z, Color: color.RGBA{R: 255, G: 255, B: 255, A: 255}})
		}
	}

	// Read face definitions
	for i := 0; i < faceCount; i++ {
		if !scanner.Scan() {
			return nil, fmt.Errorf("unexpected end of file while reading faces")
		}
		parts := strings.Fields(scanner.Text())
		numFaceVerts, _ := strconv.Atoi(parts[0])

		var faceColor color.RGBA

		// Conditional Face Parsing
		if hasFaceColor {
			if len(parts) != numFaceVerts+1+3 {
				return nil, fmt.Errorf("invalid face-color data on line %d", i)
			}
			// Color is at the end of the line
			r, _ := strconv.ParseUint(parts[numFaceVerts+1], 10, 8)
			g, _ := strconv.ParseUint(parts[numFaceVerts+2], 10, 8)
			b, _ := strconv.ParseUint(parts[numFaceVerts+3], 10, 8)
			faceColor = color.RGBA{R: uint8(r), G: uint8(g), B: uint8(b), A: 255}
		} else if hasVertexColor {
			// No face color, but vertices have color -> average them
			var r, g, b uint32
			for j := 1; j <= numFaceVerts; j++ {
				idx, _ := strconv.Atoi(parts[j])
				vert := vertices[idx]
				r += uint32(vert.Color.R)
				g += uint32(vert.Color.G)
				b += uint32(vert.Color.B)
			}
			faceColor = color.RGBA{
				R: uint8(r / uint32(numFaceVerts)),
				G: uint8(g / uint32(numFaceVerts)),
				B: uint8(b / uint32(numFaceVerts)),
				A: 255,
			}
		} else {
			// No color information at all, use a default color
			faceColor = color.RGBA{R: 128, G: 128, B: 128, A: 255}
		}

		aFace := NewFace(nil, faceColor, nil)
		for j := 1; j <= numFaceVerts; j++ {
			idx, _ := strconv.Atoi(parts[j])
			vert := vertices[idx]
			aFace.AddPoint(vert.X, vert.Y, vert.Z)
		}

		aFace.Finished(reverse)
		obj.theFaces.AddFace(aFace)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading from PLY source: %w", err)
	}

	obj.Finished(false, useBsp)
	return obj, nil
}

// A small epsilon value for floating-point comparisons to avoid precision errors.
const epsilon = 1e-6

// subtract calculates the vector from p2 to p1 (p1 - p2).
func subtract(p1, p2 *Point3d) *Point3d {
	return &Point3d{X: p1.X - p2.X, Y: p1.Y - p2.Y, Z: p1.Z - p2.Z}
}

// crossProduct calculates the cross product of two vectors.
func crossProduct(v1, v2 *Point3d) *Point3d {
	return &Point3d{
		X: v1.Y*v2.Z - v1.Z*v2.Y,
		Y: v1.Z*v2.X - v1.X*v2.Z,
		Z: v1.X*v2.Y - v1.Y*v2.X,
	}
}

// dotProduct calculates the dot product of two vectors.
func dotProduct(v1, v2 *Point3d) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}

// LineIntersectsPolygon determines if a line segment intersects with a 3D polygon.
// The polygon is assumed to be planar and convex.
func LineIntersectsPolygon(lineStart, lineEnd *Point3d, polygonPoints [][]float64) bool {
	if len(polygonPoints) < 3 {
		// A polygon must have at least 3 vertices.
		return false
	}

	// 1. Define the plane of the polygon using the first three points.
	// We need a point on the plane (p0) and the plane's normal vector.
	p0 := NewPoint3d(polygonPoints[0][0], polygonPoints[0][1], polygonPoints[0][2])
	p1 := NewPoint3d(polygonPoints[1][0], polygonPoints[1][1], polygonPoints[1][2])
	p2 := NewPoint3d(polygonPoints[2][0], polygonPoints[2][1], polygonPoints[2][2])

	// Create two vectors on the plane.
	v1 := subtract(p1, p0)
	v2 := subtract(p2, p0)

	// The plane normal is the cross product of the two vectors.
	planeNormal := crossProduct(v1, v2)

	// 2. Calculate the intersection of the line and the plane.
	// The line is represented as P = lineStart + t * lineDir
	lineDir := subtract(lineEnd, lineStart)

	// Check if the line is parallel to the plane.
	dotNormalDir := dotProduct(planeNormal, lineDir)
	if math.Abs(dotNormalDir) < epsilon {
		return false // Line is parallel, no intersection.
	}

	// Calculate the 't' parameter for the line equation.
	w := subtract(lineStart, p0)
	t := -dotProduct(planeNormal, w) / dotNormalDir

	// 3. Check if the intersection point is within the line segment.
	// If t is not between 0 and 1, the intersection is outside the segment.
	if t < 0.0-epsilon || t > 1.0+epsilon {
		return false
	}

	// 4. Calculate the actual intersection point.
	intersectionPoint := NewPoint3d(
		lineStart.X+t*lineDir.X,
		lineStart.Y+t*lineDir.Y,
		lineStart.Z+t*lineDir.Z,
	)

	// 5. Check if the intersection point is inside the polygon.
	return isPointInPolygon(intersectionPoint, polygonPoints, planeNormal)
}

// isPointInPolygon checks if a 3D point (known to be on the polygon's plane)
// is inside the polygon's boundaries using a 2D projection and ray casting.
func isPointInPolygon(point *Point3d, polygonPoints [][]float64, normal *Point3d) bool {
	// Project the 3D polygon and the point to a 2D plane.
	// We choose the plane based on the largest component of the normal vector
	// to avoid a degenerate projection (i.e., the polygon projecting to a line).
	absX := math.Abs(normal.X)
	absY := math.Abs(normal.Y)
	absZ := math.Abs(normal.Z)

	var u, v int // Indices for the 2D coordinates (0=X, 1=Y, 2=Z)

	if absX > absY && absX > absZ {
		// Project to YZ plane (discarding X)
		u, v = 1, 2
	} else if absY > absX && absY > absZ {
		// Project to XZ plane (discarding Y)
		u, v = 0, 2
	} else {
		// Project to XY plane (discarding Z)
		u, v = 0, 1
	}

	// Get the 2D coordinates of the intersection point.
	point2D := []float64{getCoord(point, u), getCoord(point, v)}

	// Apply the Ray Casting algorithm in 2D.
	intersections := 0
	numVertices := len(polygonPoints)
	for i := 0; i < numVertices; i++ {
		p1_3D := polygonPoints[i]
		p2_3D := polygonPoints[(i+1)%numVertices]

		// Get 2D coordinates of the edge's vertices.
		p1_2D := []float64{p1_3D[u], p1_3D[v]}
		p2_2D := []float64{p2_3D[u], p2_3D[v]}

		// Check if the horizontal ray from the point intersects with the edge.
		if (p1_2D[1] > point2D[1]) != (p2_2D[1] > point2D[1]) {
			// Calculate the x-intersection of the line.
			x_intersection := (p2_2D[0]-p1_2D[0])*(point2D[1]-p1_2D[1])/(p2_2D[1]-p1_2D[1]) + p1_2D[0]
			if point2D[0] < x_intersection {
				intersections++
			}
		}
	}

	// If the number of intersections is odd, the point is inside the polygon.
	return intersections%2 == 1
}

// getCoord is a small helper to get a coordinate by its index (0=X, 1=Y, 2=Z).
func getCoord(p *Point3d, index int) float64 {
	switch index {
	case 0:
		return p.X
	case 1:
		return p.Y
	case 2:
		return p.Z
	}
	return 0
}

// func (o *Model) SaveDXF(batcher *PolygonBatcher, x, y int, screenHeight, screenWidth float32) {
// 	points := make([][]float64, 0, 10)
// 	for i := 0; i < len(o.faceIndicies); i++ {
// 		faceIndices := o.faceIndicies[i]
// 		normalIndex := o.normalIndicies[i]
// 		facePointsInObjectSpace := points[:]
// 		for _, index := range faceIndices {
// 			point := o.faceMesh.Points.ThisMatrix[index]
// 			facePointsInObjectSpace = append(facePointsInObjectSpace, point)
// 		}
// 		face := o.theFaces.faces[i]
// 		normal := o.transNormalMesh.Points.ThisMatrix[normalIndex]

// 		//finish
// 	}
// }

// SaveDXF writes the model's base geometry to a file in DXF format.
// It correctly retrieves face data whether the model uses a BSP tree or a simple face list.
// Note: This function replaces the original stub, which had a signature and logic
// unsuited for file-saving operations.
func (o *Model) SaveDXF(fileName string) error {
	file, err := os.Create(fileName)
	if err != nil {
		return fmt.Errorf("could not create DXF file %s: %w", fileName, err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// writePair is a helper function to write a DXF group code and its value.
	writePair := func(code int, value interface{}) {
		// This will handle potential errors during write, though we'll check the final flush.
		_, _ = fmt.Fprintf(writer, "%d\n%v\n", code, value)
	}

	// --- 1. Collect all face indices from the model ---
	var allFaceIndices [][]int

	if o.canPaintWithoutBSP {
		// If no BSP tree is used, the indices are stored in a simple list.
		allFaceIndices = o.faceIndicies
	} else if o.root != nil {
		// If a BSP tree is used, we must traverse it to collect the face indices from each node.
		// This assumes BspNode has fields: Left *BspNode, Right *BspNode, and pointIndices []int.
		allFaceIndices = make([][]int, 0)
		var collectIndices func(node *BspNode)
		collectIndices = func(node *BspNode) {
			if node == nil {
				return
			}
			// Each node in the BSP tree represents a polygon.
			allFaceIndices = append(allFaceIndices, node.facePointIndices)
			// Recursively traverse the child nodes.
			collectIndices(node.Left)
			collectIndices(node.Right)
		}
		collectIndices(o.root)
	}

	// --- 2. Write the standard DXF file structure ---

	// Header Section (minimal)
	writePair(0, "SECTION")
	writePair(2, "HEADER")
	writePair(0, "ENDSEC")

	// Entities Section
	writePair(0, "SECTION")
	writePair(2, "ENTITIES")

	// --- 3. Write each face from the collected indices as a 3DFACE entity ---
	for _, faceIdxs := range allFaceIndices {
		if len(faceIdxs) < 3 {
			continue // A face must have at least 3 vertices.
		}

		writePair(0, "3DFACE")
		writePair(8, "0") // Layer name

		// A 3DFACE entity requires 4 vertices. For triangles, the 4th vertex repeats the 3rd.
		p1 := o.faceMesh.Points.ThisMatrix[faceIdxs[0]]
		p2 := o.faceMesh.Points.ThisMatrix[faceIdxs[1]]
		p3 := o.faceMesh.Points.ThisMatrix[faceIdxs[2]]

		var p4 []float64
		if len(faceIdxs) > 3 {
			p4 = o.faceMesh.Points.ThisMatrix[faceIdxs[3]]
		} else {
			p4 = p3 // For a triangle, repeat the third vertex.
		}

		// Write vertex coordinates with their corresponding group codes.
		writePair(10, p1[0]) // Vertex 1: X
		writePair(20, p1[1]) // Vertex 1: Y
		writePair(30, p1[2]) // Vertex 1: Z

		writePair(11, p2[0]) // Vertex 2: X
		writePair(21, p2[1]) // Vertex 2: Y
		writePair(31, p2[2]) // Vertex 2: Z

		writePair(12, p3[0]) // Vertex 3: X
		writePair(22, p3[1]) // Vertex 3: Y
		writePair(32, p3[2]) // Vertex 3: Z

		writePair(13, p4[0]) // Vertex 4: X
		writePair(23, p4[1]) // Vertex 4: Y
		writePair(33, p4[2]) // Vertex 4: Z
	}

	// End of Entities Section
	writePair(0, "ENDSEC")

	// End of File marker
	writePair(0, "EOF")

	// Flush the writer to ensure all data is written to the file and return any error.
	return writer.Flush()
}

// func (o *Model) SavePLY(fileName string) error {
// 	file, err := os.Create(fileName)
// 	if err != nil {
// 		return fmt.Errorf("could not create PLY file %s: %w", fileName, err)
// 	}
// 	defer file.Close()

// 	writer := bufio.NewWriter(file)

// 	// Helper struct to associate face indices with their color.
// 	type coloredFace struct {
// 		indices []int
// 		color   color.RGBA
// 	}

// 	// --- 1. Collect all face indices and their corresponding colors ---
// 	allFaces := make([]coloredFace, 0)

// 	if o.canPaintWithoutBSP {
// 		// If no BSP tree, get faces from the simple lists.
// 		for i, faceIdxs := range o.faceIndicies {
// 			if i < len(o.theFaces.faces) {
// 				faceColor := o.theFaces.faces[i].Col
// 				allFaces = append(allFaces, coloredFace{indices: faceIdxs, color: faceColor})
// 			}
// 		}
// 	} else if o.root != nil {
// 		// If a BSP tree is used, traverse it to collect faces.
// 		var collectFromBSP func(node *BspNode)
// 		collectFromBSP = func(node *BspNode) {
// 			if node == nil {
// 				return
// 			}
// 			// Get color and indices from the node.
// 			nodeColor := color.RGBA{R: node.colRed, G: node.colGreen, B: node.colBlue, A: node.colAlpha}
// 			allFaces = append(allFaces, coloredFace{indices: node.facePointIndices, color: nodeColor})

// 			// Recurse on children.
// 			collectFromBSP(node.Left)
// 			collectFromBSP(node.Right)
// 		}
// 		collectFromBSP(o.root)
// 	}

// 	// --- 2. Create a color map for each unique vertex ---
// 	numVertices := len(o.faceMesh.Points.ThisMatrix)
// 	vertexColors := make([]color.RGBA, numVertices)
// 	isColorSet := make([]bool, numVertices) // To track if a vertex color has been assigned.

// 	for _, face := range allFaces {
// 		for _, vertexIndex := range face.indices {
// 			if vertexIndex < numVertices && !isColorSet[vertexIndex] {
// 				vertexColors[vertexIndex] = face.color
// 				isColorSet[vertexIndex] = true
// 			}
// 		}
// 	}

// 	// --- 3. Write the PLY file header ---
// 	_, _ = fmt.Fprintln(writer, "ply")
// 	_, _ = fmt.Fprintln(writer, "format ascii 1.0")
// 	_, _ = fmt.Fprintf(writer, "comment Generated by gosie3d for model with %d vertices and %d faces\n", numVertices, len(allFaces))
// 	_, _ = fmt.Fprintf(writer, "element vertex %d\n", numVertices)
// 	_, _ = fmt.Fprintln(writer, "property float x")
// 	_, _ = fmt.Fprintln(writer, "property float y")
// 	_, _ = fmt.Fprintln(writer, "property float z")
// 	_, _ = fmt.Fprintln(writer, "property uchar red")
// 	_, _ = fmt.Fprintln(writer, "property uchar green")
// 	_, _ = fmt.Fprintln(writer, "property uchar blue")
// 	_, _ = fmt.Fprintf(writer, "element face %d\n", len(allFaces))
// 	_, _ = fmt.Fprintln(writer, "property list uchar int vertex_indices")
// 	_, _ = fmt.Fprintln(writer, "end_header")

// 	// Write vertex data (coordinates and color)
// 	for i, vertex := range o.faceMesh.Points.ThisMatrix {
// 		col := vertexColors[i]
// 		_, _ = fmt.Fprintf(writer, "%f %f %f %d %d %d\n", vertex[0], vertex[1], vertex[2], col.R, col.G, col.B)
// 	}

// 	// write face data (vertex count and indices)
// 	for _, face := range allFaces {
// 		numFaceVerts := len(face.indices)
// 		if numFaceVerts < 3 {
// 			continue
// 		}
// 		// Start line with number of vertices for this face.
// 		_, _ = fmt.Fprintf(writer, "%d", numFaceVerts)

// 		// Append the vertex indices.
// 		for _, vIndex := range face.indices {
// 			_, _ = fmt.Fprintf(writer, " %d", vIndex)
// 		}
// 		_, _ = fmt.Fprintln(writer)
// 	}

// 	return writer.Flush()
// }

func (o *Model) SavePLYWithFaceColors(fileName string) error {
	file, err := os.Create(fileName)
	if err != nil {
		return fmt.Errorf("could not create PLY file %s: %w", fileName, err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// Helper struct to associate face indices with their color.
	type coloredFace struct {
		indices []int
		color   color.RGBA
	}

	// Collect all face indices and their corresponding colors
	allFaces := make([]coloredFace, 0)
	if o.canPaintWithoutBSP {
		for i, faceIdxs := range o.faceIndicies {
			if i < len(o.theFaces.faces) {
				faceColor := o.theFaces.faces[i].Col
				allFaces = append(allFaces, coloredFace{indices: faceIdxs, color: faceColor})
			}
		}
	} else if o.root != nil {
		var collectFromBSP func(node *BspNode)
		collectFromBSP = func(node *BspNode) {
			if node == nil {
				return
			}
			nodeColor := color.RGBA{R: node.colRed, G: node.colGreen, B: node.colBlue, A: node.colAlpha}
			allFaces = append(allFaces, coloredFace{indices: node.facePointIndices, color: nodeColor})
			collectFromBSP(node.Left)
			collectFromBSP(node.Right)
		}
		collectFromBSP(o.root)
	}

	numVertices := len(o.faceMesh.Points.ThisMatrix)

	// Write the modified PLY file header
	// Note that the color properties are now part of the "face" element.
	_, _ = fmt.Fprintln(writer, "ply")
	_, _ = fmt.Fprintln(writer, "format ascii 1.0")
	_, _ = fmt.Fprintf(writer, "comment Generated by gosie3d with face colors\n")
	_, _ = fmt.Fprintf(writer, "element vertex %d\n", numVertices)
	_, _ = fmt.Fprintln(writer, "property float x")
	_, _ = fmt.Fprintln(writer, "property float y")
	_, _ = fmt.Fprintln(writer, "property float z")
	_, _ = fmt.Fprintf(writer, "element face %d\n", len(allFaces))
	_, _ = fmt.Fprintln(writer, "property list uchar int vertex_indices")
	_, _ = fmt.Fprintln(writer, "property uchar red")   // <-- Moved to face element
	_, _ = fmt.Fprintln(writer, "property uchar green") // <-- Moved to face element
	_, _ = fmt.Fprintln(writer, "property uchar blue")  // <-- Moved to face element
	_, _ = fmt.Fprintln(writer, "end_header")

	// --- 3. Write vertex data (coordinates ONLY) ---
	for _, vertex := range o.faceMesh.Points.ThisMatrix {
		_, _ = fmt.Fprintf(writer, "%f %f %f\n", vertex[0], vertex[1], vertex[2])
	}

	// --- 4. Write face data (indices AND color) ---
	for _, face := range allFaces {
		numFaceVerts := len(face.indices)
		if numFaceVerts < 3 {
			continue
		}

		// Start line with number of vertices for this face.
		_, _ = fmt.Fprintf(writer, "%d", numFaceVerts)

		// Append the vertex indices.
		for _, vIndex := range face.indices {
			_, _ = fmt.Fprintf(writer, " %d", vIndex)
		}

		// Append the face color at the end of the line.
		col := face.color
		_, _ = fmt.Fprintf(writer, " %d %d %d\n", col.R, col.G, col.B)
	}

	return writer.Flush()
}
