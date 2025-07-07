package gosie3d

import (
	"image"
	"image/color"

	"github.com/hajimehoshi/ebiten/v2"
)

var (
	whiteImage = ebiten.NewImage(3, 3)
)

func init() {
	whiteImage.Fill(color.White)
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
	op.FillRule = ebiten.FillAll
	screen.DrawTriangles(vertices, indices, whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image), op)
}
