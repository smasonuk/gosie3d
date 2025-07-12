package gosie3d

import "math"

type Vector3 struct {
	// Normal [4]float64
	X float64
	Y float64
	Z float64
	W float64 // W component, typically 1.0 for vectors
}

// func (v *Vector3) GetX() float64 { return v.Normal[0] }
// func (v *Vector3) GetY() float64 { return v.Normal[1] }
// func (v *Vector3) GetZ() float64 { return v.Normal[2] }
func (v *Vector3) Add(x, y, z float64) {
	// v.Normal[0] += x
	// v.Normal[1] += y
	// v.Normal[2] += z
	v.X += x
	v.Y += y
	v.Z += z
}

func NewVector3(x, y, z float64) *Vector3 {
	return &Vector3{
		// Normal: [4]float64{x, y, z, 1.0},
		X: x,
		Y: y,
		Z: z,
		W: 1.0, // W component is typically 1.0 for vectors
	}
}

func NewVector3Full(x, y, z, w float64) *Vector3 {
	return &Vector3{
		// Normal: [4]float64{x, y, z, 1.0},
		X: x,
		Y: y,
		Z: z,
		W: 1.0, // W component is typically 1.0 for vectors
	}
}

func NewVector3dFromArray(normal []float64) *Vector3 {
	v := &Vector3{}
	// copy(v.Normal[:], normal)

	v.X = normal[0]
	v.Y = normal[1]
	v.Z = normal[2]

	return v
}

func (v *Vector3) Normalize() {
	// length := math.Sqrt(math.Abs(v.Normal[0]*v.Normal[0] + v.Normal[1]*v.Normal[1] + v.Normal[2]*v.Normal[2]))
	// if length == 0 {
	// 	return
	// }
	// v.Normal[0] /= length
	// v.Normal[1] /= length
	// v.Normal[2] /= length
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
		// Normal: [4]float64{v.Normal[0], v.Normal[1], v.Normal[2], 1.0},
		X: v.X,
		Y: v.Y,
		Z: v.Z,
		W: v.W,
	}
}
