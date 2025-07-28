package gosie3d

import "math"

type Vector3 struct {
	X float64
	Y float64
	Z float64
	W float64
}

func (v *Vector3) Add(x, y, z float64) {
	v.X += x
	v.Y += y
	v.Z += z
}

func NewVector3(x, y, z float64) *Vector3 {
	return &Vector3{
		X: x,
		Y: y,
		Z: z,
		W: 1.0, // W component is typically 1.0 for vectors
	}
}

func NewVector3Full(x, y, z, w float64) *Vector3 {
	return &Vector3{
		X: x,
		Y: y,
		Z: z,
		W: 1.0, // W component is typically 1.0 for vectors
	}
}

func NewVector3dFromArray(normal []float64) *Vector3 {
	v := &Vector3{}
	v.X = normal[0]
	v.Y = normal[1]
	v.Z = normal[2]

	return v
}

func (v *Vector3) Normalize() {
	length := math.Sqrt(math.Abs(v.X*v.X + v.Y*v.Y + v.Z*v.Z))
	if length == 0 {
		return
	}
	v.X /= length
	v.Y /= length
	v.Z /= length
}

func (v *Vector3) Copy() *Vector3 {
	return &Vector3{
		X: v.X,
		Y: v.Y,
		Z: v.Z,
		W: v.W,
	}
}

// DistanceTo
func (v *Vector3) DistanceTo(other *Vector3) float64 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}
