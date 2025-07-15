package gosie3d

import (
	"image"
	"image/color"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

var (
	whiteImage = ebiten.NewImage(3, 3)
	whiteSub   *ebiten.Image
)

func init() {
	whiteImage.Fill(color.White)
	whiteSub = whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image)
}

func fillConvexPolygon(screen *ebiten.Image, xp, yp []float32, clr color.RGBA) {
	if len(xp) < 3 {
		return
	}

	indices := make([]uint16, 0, (len(xp)-2)*3)
	for i := 2; i < len(xp); i++ {
		indices = append(indices, 0, uint16(i-1), uint16(i))
	}

	vertices := make([]ebiten.Vertex, len(xp))
	cr := float32(clr.R) / 255.0
	cg := float32(clr.G) / 255.0
	cb := float32(clr.B) / 255.0
	ca := float32(clr.A) / 255.0

	for i := range xp {
		vertices[i] = ebiten.Vertex{
			DstX:   xp[i],
			DstY:   yp[i],
			SrcX:   1,
			SrcY:   1,
			ColorR: cr,
			ColorG: cg,
			ColorB: cb,
			ColorA: ca,
		}
	}

	op := &ebiten.DrawTrianglesOptions{}
	// op.FillRule = ebiten.FillAll
	// op.FillRule = ebiten.FillRuleEvenOdd
	op.AntiAlias = true
	// screen.DrawTriangles(vertices, indices, whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image), op)
	screen.DrawTriangles(vertices, indices, whiteSub, op)

}

// drawPolygonOutline draws the outline of a polygon defined by the given points.
// It uses the vector package to create a path and stroke it.
func drawPolygonOutline(screen *ebiten.Image, xp, yp []float32, strokeWidth float32, clr color.RGBA) {
	// We need at least 2 points to draw a line.
	if len(xp) < 2 {
		return
	}

	// Create a new vector path.
	var path vector.Path

	// Move to the first point to start the path.
	path.MoveTo(xp[0], yp[0])

	// Create line segments to the subsequent points.
	for i := 1; i < len(xp); i++ {
		path.LineTo(xp[i], yp[i])
	}

	// Close the path to connect the last point back to the first.
	path.Close()

	// Stroke options define how the line looks (width, join style, etc.).
	strokeOp := &vector.StrokeOptions{
		Width: strokeWidth,
		// LineJoin: vector.LineJoinMiter,
	}

	// AppendVerticesAndIndicesForStroke generates the vertices and indices needed to render the path's outline.
	// We pass nil for the destination slices to have new ones created.
	vertices, indices := path.AppendVerticesAndIndicesForStroke(nil, nil, strokeOp)

	// Convert the RGBA color to float32 values (0.0-1.0) for the vertices.
	cr := float32(clr.R) / 255.0
	cg := float32(clr.G) / 255.0
	cb := float32(clr.B) / 255.0
	ca := float32(clr.A) / 255.0

	// Apply the color to all the generated vertices.
	// The SrcX and SrcY are set to 1 to use the solid color from our whiteSubImage.
	for i := range vertices {
		vertices[i].ColorR = cr
		vertices[i].ColorG = cg
		vertices[i].ColorB = cb
		vertices[i].ColorA = ca
		vertices[i].SrcX = 1
		vertices[i].SrcY = 1
	}

	// Define the options for drawing the triangles.
	drawOp := &ebiten.DrawTrianglesOptions{
		AntiAlias: true,
	}

	// white := whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image)

	// Draw the triangles that form the line stroke.
	screen.DrawTriangles(vertices, indices, whiteSub, drawOp)
}
