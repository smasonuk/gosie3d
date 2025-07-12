package gosie3d

import (
	"image/color"
	"math"
)

type Plane struct {
	A, B, C, D float64
}

const planeThickness = 0.1

func NewPlaneFromPoint(pointOnPlane *Point3d, normal *Vector3) *Plane {
	newPlane := &Plane{
		A: normal.X,
		B: normal.Y,
		C: normal.Z,
	}
	newPlane.D = -(newPlane.A*pointOnPlane.X + newPlane.B*pointOnPlane.Y + newPlane.C*pointOnPlane.Z)
	return newPlane
}

func NewPlane(f *Face, normal *Vector3) *Plane {
	p := &Plane{
		// A: normal[0],
		// B: normal[1],
		// C: normal[2],
		A: normal.X,
		B: normal.Y,
		C: normal.Z,
	}
	p.D = -(p.A*f.Points[0][0] + p.B*f.Points[0][1] + p.C*f.Points[0][2])
	return p
}

// PointOnPlane determines the position of a point (x, y, z) relative to the plane.
// It substitutes the point's coordinates into the plane's equation.
// The result indicates:
// - Positive value: The point is on the "front" side (in the direction of the normal).
// - Negative value: The point is on the "back" side.
// - Zero: The point is on the plane (within the defined planeThickness).
func (p *Plane) PointOnPlane(x, y, z float64) float64 {
	// Calculate the signed distance from the point to the plane.
	num := p.A*x + p.B*y + p.C*z + p.D

	// Check if the point is very close to the plane, but not exactly on it.
	// If it's within the tolerance (planeThickness), we treat it as being on the plane.
	if math.Abs(num) > 0 && math.Abs(num) < planeThickness {
		return 0.0
	}
	return num
}

// lIntersect checks if a line segment defined by two points (p1, p2) intersects the plane.
// This is a prerequisite for calculating the intersection point.
func (p *Plane) lIntersect(p1, p2 *Point3d) bool {
	a := p.PointOnPlane(p1.X, p1.Y, p1.Z)
	b := p.PointOnPlane(p2.X, p2.Y, p2.Z)

	// If either point is on the plane, this method considers it not a "true" intersection
	// for the purpose of splitting. The calling logic handles this case.
	if a == 0 || b == 0 {
		return false
	}

	// An intersection occurs if and only if the two points are on opposite sides
	// of the plane (one's signed distance is positive, the other's is negative).
	return (a > 0 && b < 0) || (a < 0 && b > 0)
}

// LineIntersect calculates the exact 3D point where a line (defined by p1 and p2)
// intersects with the plane.
func (p *Plane) LineIntersect(p1, p2 *Point3d) *Point3d {
	// Determine which side of the plane each endpoint lies on.
	x1, y1, z1 := p1.X, p1.Y, p1.Z
	x2, y2, z2 := p2.X, p2.Y, p2.Z

	if !p.lIntersect(p1, p2) {
		return nil
	}

	// Calculate the denominator for the intersection formula. This is the dot product
	// of the plane's normal and the line's direction vector.
	denom := (p.A*(x2-x1) + p.B*(y2-y1) + p.C*(z2-z1))

	// If the denominator is zero, the line is parallel to the plane.
	// Since we already checked that the endpoints are on opposite sides, this
	// case is mathematically impossible, but it's a good safeguard.
	if denom == 0 {
		return nil
	}

	// Calculate the parameter 't' which represents the position of the intersection
	// point along the line segment (from p1 to p2).
	t := -(p.A*x1 + p.B*y1 + p.C*z1 + p.D) / denom

	// Use 't' to find the actual coordinates of the intersection point.
	x := x1 + (x2-x1)*t
	y := y1 + (y2-y1)*t
	z := z1 + (z2-z1)*t
	return NewPoint3d(x, y, z)
}

// FaceIntersect checks if a given face (polygon) intersects with the plane.
// It does this by checking if the face's vertices lie on opposite sides of the plane.
func (p *Plane) FaceIntersect(f *Face) bool {
	var d float64
	initialized := false
	for a := 0; a < f.Cnum; a++ {
		// Calculate the signed distance for the current vertex.
		n := p.PointOnPlane(f.Points[a][0], f.Points[a][1], f.Points[a][2])

		// Initialize 'd' with the signed distance of the first vertex.
		if !initialized {
			d = n
			initialized = true
			continue
		}

		// Check if the current vertex 'n' is on a different side of the plane
		// compared to the first vertex 'd'.
		// If the signs are different, the face must cross the plane.
		if !((d >= 0 && n >= 0) || (d <= 0 && n <= 0)) {
			return true
		}
	}

	// If all vertices are on the same side, the face does not intersect the plane.
	return false
}

// SplitFace divides a single face into two new faces along the intersection line
// with the plane.
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

		if p.lIntersect(p3d1, p3d2) {
			pointIntersect := p.LineIntersect(p3d1, p3d2)
			inter = true
			faces[currentFace].AddPoint(p3d1.X, p3d1.Y, p3d1.Z)
			if pointIntersect != nil {
				faces[currentFace].AddPoint(pointIntersect.X, pointIntersect.Y, pointIntersect.Z)
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(pointIntersect.X, pointIntersect.Y, pointIntersect.Z)
			}
		} else {
			if p.PointOnPlane(p3d1.X, p3d1.Y, p3d1.Z) == 0 {
				inter = true
				faces[currentFace].AddPoint(p3d1.X, p3d1.Y, p3d1.Z)
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(p3d1.X, p3d1.Y, p3d1.Z)
			} else {
				faces[currentFace].AddPoint(p3d1.X, p3d1.Y, p3d1.Z)
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
