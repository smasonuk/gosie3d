package gosie3d

import (
	"image/color"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
)

type BspNode struct {
	normal           *Vector3
	Left             *BspNode
	Right            *BspNode
	colRed           uint8
	colGreen         uint8
	colBlue          uint8
	facePointIndices []int
	xp               []float32
	yp               []float32
	normalIndex      int
	pointsToUse      [][]float64 //temp buffer
}

const nearPlaneZ = 25
const conversionFactor = 700
const antiAlias = true
const antiAliasLines = true

func NewBspNode(facePoints [][]float64, faceNormal *Vector3, faceColor color.RGBA, pointIndices []int, normalIdx int) *BspNode {
	b := &BspNode{
		normal:           faceNormal,
		colRed:           faceColor.R,
		colGreen:         faceColor.G,
		colBlue:          faceColor.B,
		facePointIndices: pointIndices,
		normalIndex:      normalIdx,
		xp:               make([]float32, len(facePoints)),
		yp:               make([]float32, len(facePoints)),
		pointsToUse:      make([][]float64, 0, len(facePoints)*2),
	}
	return b
}

func (b *BspNode) PaintWithoutShading(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, linesOnly bool) {
	b.PaintWithShading(screen, x, y, transPoints, transNormals, false, linesOnly)
}

func (b *BspNode) PaintWithShading(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, doShading bool, linesOnly bool) {
	if len(b.facePointIndices) == 0 {
		return
	}

	transformedNormal := transNormals.ThisMatrix[b.normalIndex]
	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

	where := transformedNormal[0]*firstTransformedPoint[0] +
		transformedNormal[1]*firstTransformedPoint[1] +
		transformedNormal[2]*firstTransformedPoint[2]

	if where <= 0 {
		if b.Left != nil {
			b.Left.PaintWithShading(screen, x, y, transPoints, transNormals, doShading, linesOnly)
		}
		if b.Right != nil {
			b.Right.PaintWithShading(screen, x, y, transPoints, transNormals, doShading, linesOnly)
		}
	} else {
		if b.Right != nil {
			b.Right.PaintWithShading(screen, x, y, transPoints, transNormals, doShading, linesOnly)
		}

		shouldReturn := b.paintPoly(screen, x, y, transPoints, transNormals, doShading, firstTransformedPoint, transformedNormal, linesOnly)
		if shouldReturn {
			return // Z-clipping occurred, no need to paint left side
		}

		if b.Left != nil {
			b.Left.PaintWithShading(screen, x, y, transPoints, transNormals, doShading, linesOnly)
		}
	}
}

func (b *BspNode) createFaceFromVertices(verticesInCameraSpace *Matrix) *Face {
	f := NewFaceEmpty(color.RGBA{R: b.colRed, G: b.colGreen, B: b.colBlue, A: 255}, b.normal)
	for _, pointIndex := range b.facePointIndices {
		pnt := verticesInCameraSpace.ThisMatrix[pointIndex]
		x := pnt[0]
		y := pnt[1]
		z := pnt[2]

		f.AddPoint(x, y, z)
	}
	f.Finished(FACE_NORMAL)

	return f
}

func (b *BspNode) getFaceInFront(face *Face) [][]float64 {
	//from 0,0,1 pointing forward
	plane := NewPlaneFromPoint(NewPoint3d(0, 0, nearPlaneZ), NewVector3(0, 0, 1))

	newFaces := plane.SplitFace(face)

	// find the face that is in front of the plane
	if len(newFaces) > 0 {
		for _, f := range newFaces {
			if f == nil {
				continue
			}
			for _, point := range f.Points {
				z := point[2]
				if z > nearPlaneZ+0.2 {
					return f.Points
				}
			}
		}
	}

	return nil
}

func getMidpoint(points [][]float64) []float64 {
	if len(points) == 0 {
		return nil
	}

	midpoint := make([]float64, 3)
	for _, point := range points {
		midpoint[0] += point[0]
		midpoint[1] += point[1]
		midpoint[2] += point[2]
	}

	midpoint[0] /= float64(len(points))
	midpoint[1] /= float64(len(points))
	midpoint[2] /= float64(len(points))

	return midpoint
}

func (b *BspNode) paintPoly(screen *ebiten.Image,
	x, y int,
	verticesInCameraSpace *Matrix,
	normalsInCameraSpace *Matrix,
	shadePoly bool,
	firstTransformedPoint []float64,
	transformedNormal []float64,
	linesOnly bool,
) bool {

	// are any of the points behind the camera?
	// TODO: this make shouldn't be done every time
	// pointsToUse := make([][]float64, 0, len(b.facePointIndices))
	pointsToUse := b.pointsToUse[:0] // Reuse the slice to avoid allocation

	numPointsBehind := 0
	for _, pointIndex := range b.facePointIndices {
		pnt := verticesInCameraSpace.ThisMatrix[pointIndex]
		if pnt[2] < nearPlaneZ {
			numPointsBehind++
		}
		pointsToUse = append(pointsToUse, pnt)
	}
	somePointsBehindAndInFront := numPointsBehind > 0 && numPointsBehind < len(b.facePointIndices)
	allPointsBehind := numPointsBehind == len(b.facePointIndices)
	midPoint := getMidpoint(pointsToUse)

	firstTransformedPoint = midPoint

	// need to split
	if somePointsBehindAndInFront {
		createFaceFromVertices := b.createFaceFromVertices(verticesInCameraSpace)
		faceInFront := b.getFaceInFront(createFaceFromVertices)
		pointsToUse = faceInFront
	}

	// If all points are behind the camera, we don't need to draw this polygon
	if allPointsBehind {
		if b.Left != nil {
			b.Left.PaintWithShading(screen, x, y,
				verticesInCameraSpace,
				normalsInCameraSpace,
				shadePoly,
				linesOnly)
		}
		return true
	}

	screenPointsX := make([]float32, len(pointsToUse))
	screenPointsY := make([]float32, len(pointsToUse))

	for i, points := range pointsToUse {
		screenPointsX[i] = float32((conversionFactor*points[0])/points[2]) + float32(x)
		screenPointsY[i] = float32((conversionFactor*points[1])/points[2]) + float32(y)
	}

	polyColor := color.RGBA{R: uint8(b.colRed), G: uint8(b.colGreen), B: uint8(b.colBlue), A: 255}

	if shadePoly {
		polyColor = b.GetColor(
			verticesInCameraSpace,
			normalsInCameraSpace,
			firstTransformedPoint,
			transformedNormal,
			polyColor)
	}

	if !linesOnly {
		fillConvexPolygon(screen, screenPointsX, screenPointsY, polyColor)
	} else {
		// fillConvexPolygon(screen, screenPointsX, screenPointsY, polyColor)
		black := color.RGBA{R: 100, G: 100, B: 100, A: 20}

		fillConvexPolygon(screen, screenPointsX, screenPointsY, polyColor)
		drawPolygonOutline(screen, screenPointsX, screenPointsY, 1.0, black)

		// drawPolygonOutline(screen, screenPointsX, screenPointsY, 1.0, polyColor)
	}

	// black := color.RGBA{R: 0, G: 0, B: 0, A: 200}
	// fillConvexPolygon(screen, screenPointsX, screenPointsY, black)
	// drawPolygonOutline(screen, screenPointsX, screenPointsY, 1.0, polyColor)

	// black := color.RGBA{R: 0, G: 0, B: 0, A: 255}
	// fillConvexPolygon(screen, screenPointsX, screenPointsY, black)
	// drawPolygonOutline(screen, screenPointsX, screenPointsY, 1.0, polyColor)

	// transPolyColor := color.RGBA{R: polyColor.R, G: polyColor.G, B: polyColor.B, A: 200}
	// fillConvexPolygon(screen, screenPointsX, screenPointsY, transPolyColor)

	return false
}

func (b *BspNode) GetColor(
	transPoints *Matrix,
	transNormals *Matrix,
	firstTransformedPoint []float64,
	transformedNormal []float64,
	polyColor color.RGBA,
) color.RGBA {

	// Define Lighting Parameters
	// The minimum brightness for any surface (e.g., 20% light).
	const ambientLight = 0.65
	// Higher values create a sharper, more focused spotlight cone.
	const spotlightConePower = 10.0
	// The amount of light available for the spotlight effect is the total (1.0)
	// minus the ambient light we've already added.
	const spotlightLightAmount = 1.0 - ambientLight

	// Diffuse Factor (based on face orientation)
	diffuseFactor := transformedNormal[2]
	if diffuseFactor < 0 {
		diffuseFactor = 0
	}

	// potlight Factor (based on face position in view)
	var spotlightFactor float64
	lenVecToPoint := GetLength(firstTransformedPoint)

	if lenVecToPoint > 0 {
		cosAngle := firstTransformedPoint[2] / lenVecToPoint
		if cosAngle < 0 {
			cosAngle = 0
		}
		spotlightFactor = math.Pow(cosAngle, spotlightConePower)
	} else {
		spotlightFactor = 1.0
	}

	// The spotlight's contribution to brightness is scaled by the amount
	// of light available for it.
	spotlightBrightness := diffuseFactor * spotlightFactor * spotlightLightAmount

	// Add the ambient light to the calculated spotlight brightness.
	// This ensures a minimum light level for the entire scene.
	finalBrightness := ambientLight + spotlightBrightness

	// Use the final brightness with the original color modification formula.
	// A brightness of 1.0 means no color change (subtract 0).
	// A brightness of 0.0 means maximum darkness (subtract 240).
	c := 240 - int(finalBrightness*240)

	min := 7
	r1 := clamp(int(b.colRed)-c, min, 255)
	g1 := clamp(int(b.colGreen)-c, min, 255)
	b1 := clamp(int(b.colBlue)-c, min, 255)
	polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}

	return polyColor
}
