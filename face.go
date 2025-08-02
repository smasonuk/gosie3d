package gosie3d

import "image/color"

type Face struct {
	Points  [][]float64
	Col     color.RGBA
	normal  *Vector3
	plane   *Plane
	vecPnts [][]float64
	meRev   bool
	Cnum    int
}

const (
	FACE_NORMAL  = 0
	FACE_REVERSE = 1
)

func NewFaceEmpty(col color.RGBA, normal *Vector3) *Face {
	return &Face{
		Points: make([][]float64, 0),
		Col:    col,
		normal: normal,
	}
}

// copy
func (f *Face) Copy() *Face {
	newFace := &Face{
		Points:  make([][]float64, len(f.Points)),
		Col:     f.Col,
		normal:  f.GetNormal().Copy(),
		plane:   f.plane,
		vecPnts: make([][]float64, len(f.vecPnts)),
		meRev:   f.meRev,
		Cnum:    f.Cnum,
	}
	for i, p := range f.Points {
		newFace.Points[i] = make([]float64, len(p))
		copy(newFace.Points[i], p)
	}
	for i, p := range f.vecPnts {
		newFace.vecPnts[i] = make([]float64, len(p))
		copy(newFace.vecPnts[i], p)
	}
	return newFace
}

func NewFace(pnts [][]float64, col color.RGBA, normal *Vector3) *Face {
	f := &Face{
		Points: pnts,
		Col:    col,
		normal: normal,
	}
	if pnts != nil {
		f.Cnum = len(pnts)
	}
	return f
}

func (f *Face) SetColor(col color.RGBA) {
	f.Col = col
}

func (f *Face) AddPoint(x, y, z float64) {
	pnts := []float64{x, y, z, 1.0}
	f.vecPnts = append(f.vecPnts, pnts)
}

func (f *Face) GetNormal() *Vector3 {
	if f.normal == nil {
		f.createNormal()
	}
	return f.normal.Copy()
}

func (f *Face) SetNormal(norm *Vector3) {
	f.normal = norm.Copy()
}

func (f *Face) Finished(reverse int) {
	if reverse == FACE_REVERSE {
		f.meRev = true
	}
	f.Cnum = len(f.vecPnts)
	f.Points = make([][]float64, f.Cnum)
	copy(f.Points, f.vecPnts)
	f.vecPnts = nil
}

func (f *Face) GetPlane() *Plane {
	if f.plane != nil {
		return f.plane
	}
	f.plane = NewPlane(f, f.GetNormal())
	return f.plane
}

func (f *Face) createNormal() {
	if len(f.Points) < 3 {
		f.normal = NewVector3(0, 0, 1)
		return
	}

	f.normal = NewVector3(0, 0, 0)
	x1, y1, z1 := f.Points[0][0], f.Points[0][1], f.Points[0][2]
	x2, y2, z2 := f.Points[1][0], f.Points[1][1], f.Points[1][2]
	x3, y3, z3 := f.Points[2][0], f.Points[2][1], f.Points[2][2]

	u1, u2, u3 := x2-x1, y2-y1, z2-z1
	v1, v2, v3 := x3-x2, y3-y2, z3-z2

	if f.meRev {
		f.normal.X = -(u2*v3 - u3*v2)
		f.normal.Y = -(u3*v1 - u1*v3)
		f.normal.Z = -(u1*v2 - u2*v1)
	} else {
		f.normal.X = u2*v3 - u3*v2
		f.normal.Y = u3*v1 - u1*v3
		f.normal.Z = u1*v2 - u2*v1
	}

	nor := f.normal.Copy()
	nor.Normalize()
	f.normal = nor.Copy()
}

// get midpoint of the face
func (f *Face) GetMidPoint() *Vector3 {
	if len(f.Points) == 0 {
		return NewVector3(0, 0, 0)
	}

	sumX, sumY, sumZ := 0.0, 0.0, 0.0
	for _, p := range f.Points {
		sumX += p[0]
		sumY += p[1]
		sumZ += p[2]
	}
	count := float64(len(f.Points))
	return NewVector3(sumX/count, sumY/count, sumZ/count)
}

// get the distance from the face to a point
func (f *Face) GetDistanceToPoint(p *Vector3) float64 {
	midpoint := f.GetMidPoint()
	return midpoint.DistanceTo(p)
}
