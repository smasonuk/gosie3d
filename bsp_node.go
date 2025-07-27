package gosie3d

import (
	"image/color"
	"math"
)

// BspNode represents a node in a Binary Space Partitioning tree.
type BspNode struct {
	normal           *Vector3
	Left             *BspNode
	Right            *BspNode
	colRed           uint8
	colGreen         uint8
	colBlue          uint8
	colAlpha         uint8
	facePointIndices []int
	xp               []float32
	yp               []float32
	normalIndex      int
	pointsToUse      [][]float64 //temp buffer
	screenPointsX    []float32
	screenPointsY    []float32
}

// Point represents a 2D point in screen space.
type Point struct {
	X, Y float32
}

const nearPlaneZ = 25
const conversionFactor = 700

// NewBspNode creates a new BSP node.
func NewBspNode(facePoints [][]float64, faceNormal *Vector3, faceColor color.RGBA, pointIndices []int, normalIdx int) *BspNode {
	b := &BspNode{
		normal:           faceNormal,
		colRed:           faceColor.R,
		colGreen:         faceColor.G,
		colBlue:          faceColor.B,
		colAlpha:         faceColor.A,
		facePointIndices: pointIndices,
		normalIndex:      normalIdx,
		xp:               make([]float32, len(facePoints)),
		yp:               make([]float32, len(facePoints)),
		pointsToUse:      make([][]float64, 0, len(facePoints)*2),
		screenPointsX:    make([]float32, 10),
		screenPointsY:    make([]float32, 10),
	}
	return b
}

// PaintWithoutShading paints the BSP tree without lighting effects.
func (b *BspNode) PaintWithoutShading(batcher *PolygonBatcher, x, y int, transPoints *Matrix, transNormals *Matrix, linesOnly bool, screenWidth, screenHeight float32) {
	b.PaintWithShading(batcher, x, y, transPoints, transNormals, false, linesOnly, screenWidth, screenHeight)
}

// PaintWithShading recursively traverses the BSP tree and paints the polygons.
func (b *BspNode) PaintWithShading(batcher *PolygonBatcher, x, y int, transPoints *Matrix, transNormals *Matrix, doShading bool, linesOnly bool, screenWidth, screenHeight float32) {
	if len(b.facePointIndices) == 0 {
		return
	}

	transformedNormal := transNormals.ThisMatrix[b.normalIndex]
	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

	// Determine if the polygon is facing the camera
	where := transformedNormal[0]*firstTransformedPoint[0] +
		transformedNormal[1]*firstTransformedPoint[1] +
		transformedNormal[2]*firstTransformedPoint[2]

	if where <= 0 { // Facing away from the camera
		if b.Left != nil {
			b.Left.PaintWithShading(batcher, x, y, transPoints, transNormals, doShading, linesOnly, screenWidth, screenHeight)
		}
		if b.Right != nil {
			b.Right.PaintWithShading(batcher, x, y, transPoints, transNormals, doShading, linesOnly, screenWidth, screenHeight)
		}
	} else { // Facing towards the camera
		if b.Right != nil {
			b.Right.PaintWithShading(batcher, x, y, transPoints, transNormals, doShading, linesOnly, screenWidth, screenHeight)
		}

		shouldReturn := b.paintPoly(batcher, x, y, transPoints, transNormals, doShading, firstTransformedPoint, transformedNormal, linesOnly, screenWidth, screenHeight)
		if shouldReturn {
			return // Z-clipping occurred, no need to paint left side
		}

		if b.Left != nil {
			b.Left.PaintWithShading(batcher, x, y, transPoints, transNormals, doShading, linesOnly, screenWidth, screenHeight)
		}
	}
}

// clipPolygonAgainstNearPlane clips a 3D polygon against the near Z plane.
// This is crucial for preventing rendering artifacts when polygons are partially behind the camera.
func clipPolygonAgainstNearPlane(polygon [][]float64) [][]float64 {
	if len(polygon) == 0 {
		return polygon
	}

	clippedPolygon := make([][]float64, 0, len(polygon)*2)
	startPoint := polygon[len(polygon)-1]

	for _, endPoint := range polygon {
		startInside := startPoint[2] >= nearPlaneZ
		endInside := endPoint[2] >= nearPlaneZ

		if endInside {
			if !startInside {
				// Edge crosses from outside to inside, calculate intersection
				intersection := intersectNearPlane(startPoint, endPoint)
				clippedPolygon = append(clippedPolygon, intersection)
			}
			// End point is inside, add it
			clippedPolygon = append(clippedPolygon, endPoint)
		} else if startInside {
			// Edge crosses from inside to outside, calculate intersection
			intersection := intersectNearPlane(startPoint, endPoint)
			clippedPolygon = append(clippedPolygon, intersection)
		}
		// If both points are outside, do nothing

		startPoint = endPoint
	}

	return clippedPolygon
}

// intersectNearPlane calculates the intersection point of a line segment with the near Z plane.
func intersectNearPlane(p1, p2 []float64) []float64 {
	// Using linear interpolation to find the intersection point
	// t = (nearPlaneZ - p1.z) / (p2.z - p1.z)
	deltaZ := p2[2] - p1[2]
	if math.Abs(deltaZ) < 1e-6 { // Avoid division by zero for horizontal lines
		return p1
	}
	t := (nearPlaneZ - p1[2]) / deltaZ

	ix := p1[0] + t*(p2[0]-p1[0])
	iy := p1[1] + t*(p2[1]-p1[1])

	return []float64{ix, iy, nearPlaneZ}
}

// paintPoly handles Z-clipping, screen-space clipping, and drawing of a single polygon.
func (b *BspNode) paintPoly(
	batcher *PolygonBatcher,
	x, y int,
	verticesInCameraSpace *Matrix,
	normalsInCameraSpace *Matrix,
	shadePoly bool,
	firstTransformedPoint []float64,
	transformedNormal []float64,
	linesOnly bool,
	screenWidth, screenHeight float32,
) bool {

	initial3DPoints := b.pointsToUse[:0] // Reuse the slice to avoid allocation
	for _, pointIndex := range b.facePointIndices {
		initial3DPoints = append(initial3DPoints, verticesInCameraSpace.ThisMatrix[pointIndex])
	}

	pointsToUse := clipPolygonAgainstNearPlane(initial3DPoints)

	// If clipping results in a polygon with too few vertices, don't draw it.
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

	polyColor := color.RGBA{R: b.colRed, G: b.colGreen, B: b.colBlue, A: b.colAlpha}
	if shadePoly {
		shadingRefPoint := verticesInCameraSpace.ThisMatrix[b.facePointIndices[0]]
		polyColor = b.GetColor(shadingRefPoint, transformedNormal, polyColor)
	}

	if !linesOnly {
		black := color.RGBA{R: 100, G: 100, B: 100, A: 25}
		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, polyColor, black, 1.0)

	} else {
		black := color.RGBA{R: 0, G: 0, B: 0, A: 255}
		// batcher.AddPolygonOutline(finalScreenPointsX, finalScreenPointsY, 1, polyColor)
		batcher.AddPolygonAndOutline(finalScreenPointsX, finalScreenPointsY, black, polyColor, 1.0)
	}

	return false
}

// clipPolygon applies the Sutherland-Hodgman algorithm to clip a polygon against the screen boundaries.
func clipPolygon(subjectPolygon []Point, screenWidth, screenHeight float32) []Point {
	// Clip against the 4 screen edges sequentially.
	// The output of one clipping stage becomes the input for the next.
	clipped := clipAgainstEdge(subjectPolygon, func(p Point) bool { return p.X >= 0 }, func(s, e Point) Point { // Left edge
		if e.X == s.X {
			return Point{X: 0, Y: s.Y}
		} // Avoid division by zero for vertical lines
		return Point{X: 0, Y: s.Y + (e.Y-s.Y)*(0-s.X)/(e.X-s.X)}
	})
	clipped = clipAgainstEdge(clipped, func(p Point) bool { return p.X <= screenWidth }, func(s, e Point) Point { // Right edge
		if e.X == s.X {
			return Point{X: screenWidth, Y: s.Y}
		}
		return Point{X: screenWidth, Y: s.Y + (e.Y-s.Y)*(screenWidth-s.X)/(e.X-s.X)}
	})
	clipped = clipAgainstEdge(clipped, func(p Point) bool { return p.Y >= 0 }, func(s, e Point) Point { // Top edge
		if e.Y == s.Y {
			return Point{X: s.X, Y: 0}
		} // Avoid division by zero for horizontal lines
		return Point{X: s.X + (e.X-s.X)*(0-s.Y)/(e.Y-s.Y), Y: 0}
	})
	clipped = clipAgainstEdge(clipped, func(p Point) bool { return p.Y <= screenHeight }, func(s, e Point) Point { // Bottom edge
		if e.Y == s.Y {
			return Point{X: s.X, Y: screenHeight}
		}
		return Point{X: s.X + (e.X-s.X)*(screenHeight-s.Y)/(e.Y-s.Y), Y: screenHeight}
	})

	return clipped
}

// clipAgainstEdge clips a polygon against a single arbitrary edge.
// The edge is defined by the `inside` function, which returns true if a point is on the visible side.
// The `intersection` function calculates where a line segment crosses the edge boundary.
func clipAgainstEdge(subjectPolygon []Point, inside func(Point) bool, intersection func(Point, Point) Point) []Point {
	if len(subjectPolygon) == 0 {
		return subjectPolygon
	}

	outputList := make([]Point, 0, len(subjectPolygon)*2)
	s := subjectPolygon[len(subjectPolygon)-1] // Start with the last vertex to process the edge from last to first vertex

	for _, e := range subjectPolygon {
		sInside := inside(s)
		eInside := inside(e)

		if eInside {
			if !sInside {
				// Was outside, is now inside: crossed the boundary, so add intersection point.
				outputList = append(outputList, intersection(s, e))
			}
			// Add the vertex 'e' as it is inside.
			outputList = append(outputList, e)
		} else if sInside {
			// Was inside, is now outside: crossed the boundary, so add intersection point.
			outputList = append(outputList, intersection(s, e))
		}
		// If both were outside, do nothing.

		s = e // Advance to the next edge
	}

	return outputList
}

// GetColor calculates the color of a polygon based on simple lighting.
func (b *BspNode) GetColor(
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
	r1 := clamp2(int(b.colRed)-c, min, 255)
	g1 := clamp2(int(b.colGreen)-c, min, 255)
	b1 := clamp2(int(b.colBlue)-c, min, 255)
	polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: polyColor.A}

	return polyColor
}

// --- Helper functions (assumed to exist elsewhere, added for completeness) ---

// clamp2 restricts an integer value to a given range.
func clamp2(value, min, max int) int {
	if value < min {
		return min
	}
	if value > max {
		return max
	}
	return value
}

// GetLength2 calculates the magnitude of a 3D vector.
func GetLength2(v []float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}
