package gosie3d

import (
	"image/color"
	"math"
)

type Plane struct {
	A, B, C, D float64
}

const planeThickness = 0.1

func NewPlane(f *Face, normal []float64) *Plane {
	p := &Plane{
		A: normal[0],
		B: normal[1],
		C: normal[2],
	}
	p.D = -(p.A*f.Points[0][0] + p.B*f.Points[0][1] + p.C*f.Points[0][2])
	return p
}

func (p *Plane) PointOnPlane(x, y, z float64) float64 {
	num := p.A*x + p.B*y + p.C*z + p.D
	if math.Abs(num) > 0 && math.Abs(num) < planeThickness {
		return 0.0
	}
	return num
}

func (p *Plane) LIntersect(p1, p2 *Point3d) bool {
	a := p.PointOnPlane(p1.GetX(), p1.GetY(), p1.GetZ())
	b := p.PointOnPlane(p2.GetX(), p2.GetY(), p2.GetZ())
	if a == 0 || b == 0 {
		return false
	}
	return (a > 0 && b < 0) || (a < 0 && b > 0)
}

func (p *Plane) LineIntersect(p1, p2 *Point3d) *Point3d {
	x1, y1, z1 := p1.GetX(), p1.GetY(), p1.GetZ()
	x2, y2, z2 := p2.GetX(), p2.GetY(), p2.GetZ()

	if !p.LIntersect(p1, p2) {
		return nil
	}
	denom := (p.A*(x2-x1) + p.B*(y2-y1) + p.C*(z2-z1))
	if denom == 0 {
		return nil
	}
	t := -(p.A*x1 + p.B*y1 + p.C*z1 + p.D) / denom
	x := x1 + (x2-x1)*t
	y := y1 + (y2-y1)*t
	z := z1 + (z2-z1)*t
	return NewPoint3d(x, y, z)
}

func (p *Plane) FaceIntersect(f *Face) bool {
	var d float64
	initialized := false
	for a := 0; a < f.Cnum; a++ {
		n := p.PointOnPlane(f.Points[a][0], f.Points[a][1], f.Points[a][2])
		if !initialized {
			d = n
			initialized = true
			continue
		}
		if !((d >= 0 && n >= 0) || (d <= 0 && n <= 0)) {
			return true
		}
	}
	return false
}

func (p *Plane) SplitFace(aFace *Face) []*Face {
	faces := make([]*Face, 2)
	faces[0] = NewFace(nil, color.RGBA{}, nil)
	faces[1] = NewFace(nil, color.RGBA{}, nil)
	var inter bool

	if !p.FaceIntersect(aFace) {
		faces[0] = aFace
		faces[1] = nil
		return faces
	}

	currentFace := 0
	pnts := NewClist(aFace.Cnum)
	for i := 0; i < aFace.Cnum; i++ {
		pnts.AddPoint(NewPoint3d(aFace.Points[i][0], aFace.Points[i][1], aFace.Points[i][2]))
	}

	for pnt := 0; pnt < aFace.Cnum; pnt++ {
		p3d1 := pnts.NextPoint()
		p3d2 := pnts.NextPoint()
		pnts.Back()

		if p.LIntersect(p3d1, p3d2) {
			pointIntersect := p.LineIntersect(p3d1, p3d2)
			inter = true
			faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			if pointIntersect != nil {
				faces[currentFace].AddPoint(pointIntersect.GetX(), pointIntersect.GetY(), pointIntersect.GetZ())
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(pointIntersect.GetX(), pointIntersect.GetY(), pointIntersect.GetZ())
			}
		} else {
			if p.PointOnPlane(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ()) == 0 {
				inter = true
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			} else {
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			}
		}
	}

	if !inter {
		faces[0] = aFace
		faces[1] = nil
		return faces
	}

	faces[0].Finished(FACE_NORMAL)
	faces[1].Finished(FACE_NORMAL)
	return faces
}

func (p *Plane) Where(f *Face) float64 {
	var inter float64
	for i := 0; i < len(f.Points); i++ {
		inter += p.PointOnPlane(f.Points[i][0], f.Points[i][1], f.Points[i][2])
	}
	return inter
}
