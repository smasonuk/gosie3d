package gosie3d

type Point3d struct {
	// Points []float64
	X float64
	Y float64
	Z float64
}

func NewPoint3d(x, y, z float64) *Point3d {
	// return &Point3d{
	// 	Points: []float64{x, y, z, 1.0},
	// }
	return &Point3d{
		X: x,
		Y: y,
		Z: z,
	}
}

// func (p *Point3d) GetX() float64 { return p.Points[0] }
// func (p *Point3d) GetY() float64 { return p.Points[1] }
// func (p *Point3d) GetZ() float64 { return p.Points[2] }
