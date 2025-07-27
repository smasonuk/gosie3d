package gosie3d

import (
	"image"
	"image/color"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

var polygonBuffer struct {
	vertices []ebiten.Vertex
	indices  []uint16
}

var (
	whiteImage = ebiten.NewImage(3, 3)
	whiteSub   *ebiten.Image
)

func init() {
	whiteImage.Fill(color.White)
	whiteSub = whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image)
}

type PolygonBatcher struct {
	vertices []ebiten.Vertex
	indices  []uint16
}

func NewPolygonBatcher(initialCap int) *PolygonBatcher {
	return &PolygonBatcher{
		vertices: make([]ebiten.Vertex, 0, initialCap*4), // Guessing 4 vertices per polygon on avg
		indices:  make([]uint16, 0, initialCap*6),        // Guessing 6 indices per polygon on avg
	}
}

// AddPolygon adds a single polygon's geometry to the batch.
func (b *PolygonBatcher) AddPolygon(xp, yp []float32, clr color.RGBA) {
	if len(xp) < 3 {
		return
	}

	// This offset is crucial! It's the starting index for the new vertices.
	vertexOffset := uint16(len(b.vertices))

	cr := float32(clr.R) / 255
	cg := float32(clr.G) / 255
	cb := float32(clr.B) / 255
	ca := float32(clr.A) / 255

	// Add the new vertices to the main slice
	for i := range xp {
		v := ebiten.Vertex{
			DstX:   xp[i],
			DstY:   yp[i],
			SrcX:   1, // From your white texture
			SrcY:   1,
			ColorR: cr,
			ColorG: cg,
			ColorB: cb,
			ColorA: ca,
		}
		b.vertices = append(b.vertices, v)
	}

	// Add the new indices, making sure to apply the vertexOffset
	for i := 2; i < len(xp); i++ {
		b.indices = append(b.indices,
			vertexOffset,             // First vertex of this polygon
			vertexOffset+uint16(i-1), // Previous vertex
			vertexOffset+uint16(i),   // Current vertex
		)
	}
}

// Draw sends the entire batch of polygons to the GPU in a single draw call.
func (b *PolygonBatcher) Draw(screen *ebiten.Image, whiteSub *ebiten.Image) {
	// Don't do anything if there's no geometry to draw.
	if len(b.vertices) == 0 {
		return
	}

	op := &ebiten.DrawTrianglesOptions{
		AntiAlias: true,
	}

	screen.DrawTriangles(b.vertices, b.indices, whiteSub, op)

	// Reset slices for the next frame, but keep the allocated memory. ✨
	b.vertices = b.vertices[:0]
	b.indices = b.indices[:0]
}

func (b *PolygonBatcher) AddPolygonOutline(xp, yp []float32, strokeWidth float32, clr color.RGBA) {
	if len(xp) < 2 {
		return
	}

	// This offset is the starting index for the vertices we are about to add.
	vertexOffset := len(b.vertices)

	var path vector.Path
	path.MoveTo(xp[0], yp[0])
	for i := 1; i < len(xp); i++ {
		path.LineTo(xp[i], yp[i])
	}
	path.Close()

	strokeOp := &vector.StrokeOptions{
		Width:    strokeWidth,
		LineJoin: vector.LineJoinMiter, // Or another join style
	}

	// Append the new vertices and indices directly to the batcher's slices.
	// This is the core of the operation.
	b.vertices, b.indices = path.AppendVerticesAndIndicesForStroke(b.vertices, b.indices, strokeOp)

	// Now, color only the vertices we just added.
	cr, cg, cb, ca := float32(clr.R)/255, float32(clr.G)/255, float32(clr.B)/255, float32(clr.A)/255
	for i := vertexOffset; i < len(b.vertices); i++ {
		b.vertices[i].ColorR = cr
		b.vertices[i].ColorG = cg
		b.vertices[i].ColorB = cb
		b.vertices[i].ColorA = ca
		b.vertices[i].SrcX = 1
		b.vertices[i].SrcY = 1
	}
}

// AddPolygonWithOutline adds a filled polygon and its outline to the batch.
// This is more efficient than calling separate fill and stroke functions.
func (b *PolygonBatcher) AddPolygonAndOutline(xp, yp []float32, fillClr, strokeClr color.RGBA, strokeWidth float32) {
	if len(xp) < 3 {
		return // Need at least 3 vertices for a polygon.
	}

	// --- 1. Add the fill geometry ---
	fillVertexOffset := uint16(len(b.vertices))

	// Pre-calculate fill color components to avoid division in the loop.
	fr, fg, fb, fa := float32(fillClr.R)/255, float32(fillClr.G)/255, float32(fillClr.B)/255, float32(fillClr.A)/255

	for i := range xp {
		v := ebiten.Vertex{
			DstX:   xp[i],
			DstY:   yp[i],
			SrcX:   1, // Corresponds to the 1x1 whiteSub image
			SrcY:   1,
			ColorR: fr,
			ColorG: fg,
			ColorB: fb,
			ColorA: fa,
		}
		b.vertices = append(b.vertices, v)
	}

	// Triangulate the polygon using a simple fan triangulation.
	for i := 2; i < len(xp); i++ {
		b.indices = append(b.indices,
			fillVertexOffset,             // First vertex
			fillVertexOffset+uint16(i-1), // Previous vertex
			fillVertexOffset+uint16(i),   // Current vertex
		)
	}

	// --- 2. Add the stroke geometry ---
	strokeVertexOffset := len(b.vertices)

	var path vector.Path
	path.MoveTo(xp[0], yp[0])
	for i := 1; i < len(xp); i++ {
		path.LineTo(xp[i], yp[i])
	}
	path.Close()

	strokeOp := &vector.StrokeOptions{
		Width:    strokeWidth,
		LineJoin: vector.LineJoinMiter,
	}

	// `AppendVerticesAndIndicesForStroke` efficiently generates the outline mesh.
	b.vertices, b.indices = path.AppendVerticesAndIndicesForStroke(b.vertices, b.indices, strokeOp)

	// Pre-calculate stroke color components.
	sr, sg, sb, sa := float32(strokeClr.R)/255, float32(strokeClr.G)/255, float32(strokeClr.B)/255, float32(strokeClr.A)/255

	// Color only the new vertices added for the stroke.
	for i := strokeVertexOffset; i < len(b.vertices); i++ {
		b.vertices[i].ColorR = sr
		b.vertices[i].ColorG = sg
		b.vertices[i].ColorB = sb
		b.vertices[i].ColorA = sa
		b.vertices[i].SrcX = 1
		b.vertices[i].SrcY = 1
	}
}
