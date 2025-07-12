package gosie3d

type Point3d struct {
	X float64
	Y float64
	Z float64
}

func NewPoint3d(x, y, z float64) *Point3d {
	return &Point3d{
		X: x,
		Y: y,
		Z: z,
	}
}
