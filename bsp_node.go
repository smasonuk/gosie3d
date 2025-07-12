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
}

// todo: change Normal to a Vector3 type
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
	}
	return b
}

func (b *BspNode) PaintWithoutColorChange(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix) {
	b.PaintWithColor(screen, x, y, transPoints, transNormals, false)
}

func (b *BspNode) PaintWithColor(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, changeColor bool) {
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
			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
		if b.Right != nil {
			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
	} else {
		if b.Right != nil {
			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}

		// Project the points into screen space
		for i, pointIndex := range b.facePointIndices {
			pnt := transPoints.ThisMatrix[pointIndex]
			if pnt[2] < 0 { // Z-clipping
				if b.Left != nil {
					b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
				}
				return
			}
			b.xp[i] = float32((400*pnt[0])/pnt[2]) + float32(x)
			b.yp[i] = float32((400*pnt[1])/pnt[2]) + float32(y)
		}

		polyColor := color.RGBA{R: uint8(b.colRed), G: uint8(b.colGreen), B: uint8(b.colBlue), A: 255}

		if changeColor {
			// Get the color based on lighting and position
			polyColor = b.GetColor(transPoints, transNormals, firstTransformedPoint, transformedNormal, polyColor)
		}

		// draw this polygon
		fillConvexPolygon(screen, b.xp, b.yp, polyColor)

		if b.Left != nil {
			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
	}
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
