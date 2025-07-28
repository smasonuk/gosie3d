package gosie3d

import (
	"image"
	"image/color"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

const (
	MaxPolygonVertices = 50000
	antiAlias          = false
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
	// vertices []ebiten.Vertex
	// indices  []uint16
	batches      []*PolygonBatch
	currentBatch *PolygonBatch
}

type PolygonBatch struct {
	vertices []ebiten.Vertex
	indices  []uint16
}

func NewPolygonBatcher(initialCap int) *PolygonBatcher {
	currentBatch := PolygonBatch{
		vertices: make([]ebiten.Vertex, 0, MaxPolygonVertices),
		indices:  make([]uint16, 0, MaxPolygonVertices*6),
	}

	batches := make([]*PolygonBatch, 0, 10) // Start with a few batches
	batches = append(batches, &currentBatch)
	return &PolygonBatcher{
		// vertices: make([]ebiten.Vertex, 0, initialCap*4), // Guessing 4 vertices per polygon on avg
		// indices:  make([]uint16, 0, initialCap*6),        // Guessing 6 indices per polygon on avg
		batches:      batches,
		currentBatch: &currentBatch,
	}
}

// add
func (b *PolygonBatcher) AddBatchIfNeeded() {
	if len(b.currentBatch.vertices) > MaxPolygonVertices {
		// If the current batch is full, start a new one.
		b.currentBatch = &PolygonBatch{
			vertices: make([]ebiten.Vertex, 0, MaxPolygonVertices),
			indices:  make([]uint16, 0, MaxPolygonVertices*6),
		}
		b.batches = append(b.batches, b.currentBatch)
	}
}

// AddPolygon adds a single polygon's geometry to the batch.
func (b *PolygonBatcher) AddPolygon(xp, yp []float32, clr color.RGBA) {
	if len(xp) < 3 {
		return
	}

	b.AddBatchIfNeeded()

	// This offset is crucial! It's the starting index for the new vertices.
	vertexOffset := uint16(len(b.currentBatch.vertices))

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
		b.currentBatch.vertices = append(b.currentBatch.vertices, v)
	}

	// Add the new indices, making sure to apply the vertexOffset
	for i := 2; i < len(xp); i++ {
		b.currentBatch.indices = append(b.currentBatch.indices,
			vertexOffset,             // First vertex of this polygon
			vertexOffset+uint16(i-1), // Previous vertex
			vertexOffset+uint16(i),   // Current vertex
		)
	}
}

// Draw sends the entire batch of polygons to the GPU in a single draw call.
func (b *PolygonBatcher) Draw(screen *ebiten.Image, whiteSub *ebiten.Image) {

	for _, batch := range b.batches {
		// Don't do anything if there's no geometry to draw.
		if len(batch.vertices) == 0 {
			return
		}

		op := &ebiten.DrawTrianglesOptions{
			AntiAlias: antiAlias,
		}

		screen.DrawTriangles(batch.vertices, batch.indices, whiteSub, op)

		// Reset slices for the next frame, but keep the allocated memory.
		batch.vertices = batch.vertices[:0]
		batch.indices = batch.indices[:0]
	}

	b.currentBatch = &PolygonBatch{
		vertices: make([]ebiten.Vertex, 0, MaxPolygonVertices),
		indices:  make([]uint16, 0, MaxPolygonVertices*6),
	}
	b.batches = []*PolygonBatch{b.currentBatch}

}

func (b *PolygonBatcher) AddPolygonOutline(xp, yp []float32, strokeWidth float32, clr color.RGBA) {
	if len(xp) < 2 {
		return
	}

	b.AddBatchIfNeeded()

	// This offset is the starting index for the vertices we are about to add.
	vertexOffset := len(b.currentBatch.vertices)

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
	b.currentBatch.vertices, b.currentBatch.indices = path.AppendVerticesAndIndicesForStroke(b.currentBatch.vertices, b.currentBatch.indices, strokeOp)

	// Now, color only the vertices we just added.
	cr, cg, cb, ca := float32(clr.R)/255, float32(clr.G)/255, float32(clr.B)/255, float32(clr.A)/255
	for i := vertexOffset; i < len(b.currentBatch.vertices); i++ {
		b.currentBatch.vertices[i].ColorR = cr
		b.currentBatch.vertices[i].ColorG = cg
		b.currentBatch.vertices[i].ColorB = cb
		b.currentBatch.vertices[i].ColorA = ca
		b.currentBatch.vertices[i].SrcX = 1
		b.currentBatch.vertices[i].SrcY = 1
	}
}

// AddPolygonWithOutline adds a filled polygon and its outline to the batch.
// This is more efficient than calling separate fill and stroke functions.
func (b *PolygonBatcher) AddPolygonAndOutline(xp, yp []float32, fillClr, strokeClr color.RGBA, strokeWidth float32) {
	if len(xp) < 3 {
		return // Need at least 3 vertices for a polygon.
	}

	b.AddBatchIfNeeded()

	// --- 1. Add the fill geometry ---
	fillVertexOffset := uint16(len(b.currentBatch.vertices))

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
		b.currentBatch.vertices = append(b.currentBatch.vertices, v)
	}

	// Triangulate the polygon using a simple fan triangulation.
	for i := 2; i < len(xp); i++ {
		b.currentBatch.indices = append(b.currentBatch.indices,
			fillVertexOffset,             // First vertex
			fillVertexOffset+uint16(i-1), // Previous vertex
			fillVertexOffset+uint16(i),   // Current vertex
		)
	}

	// --- 2. Add the stroke geometry ---
	strokeVertexOffset := len(b.currentBatch.vertices)

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
	b.currentBatch.vertices, b.currentBatch.indices = path.AppendVerticesAndIndicesForStroke(b.currentBatch.vertices, b.currentBatch.indices, strokeOp)

	// Pre-calculate stroke color components.
	sr, sg, sb, sa := float32(strokeClr.R)/255, float32(strokeClr.G)/255, float32(strokeClr.B)/255, float32(strokeClr.A)/255

	// Color only the new vertices added for the stroke.
	for i := strokeVertexOffset; i < len(b.currentBatch.vertices); i++ {
		b.currentBatch.vertices[i].ColorR = sr
		b.currentBatch.vertices[i].ColorG = sg
		b.currentBatch.vertices[i].ColorB = sb
		b.currentBatch.vertices[i].ColorA = sa
		b.currentBatch.vertices[i].SrcX = 1
		b.currentBatch.vertices[i].SrcY = 1
	}
}
