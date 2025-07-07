package gosie3d

import "math"

type Vector3 struct {
	Normal [4]float64
}

func (v *Vector3) GetX() float64 { return v.Normal[0] }
func (v *Vector3) GetY() float64 { return v.Normal[1] }
func (v *Vector3) GetZ() float64 { return v.Normal[2] }
func (v *Vector3) Add(x, y, z float64) {
	v.Normal[0] += x
	v.Normal[1] += y
	v.Normal[2] += z
}

func NewVector3(x, y, z float64) *Vector3 {
	return &Vector3{
		Normal: [4]float64{x, y, z, 1.0},
	}
}

func NewVector3dFromArray(normal []float64) *Vector3 {
	v := &Vector3{}
	copy(v.Normal[:], normal)
	return v
}

func (v *Vector3) Normalize() {
	length := math.Sqrt(math.Abs(v.Normal[0]*v.Normal[0] + v.Normal[1]*v.Normal[1] + v.Normal[2]*v.Normal[2]))
	if length == 0 {
		return
	}
	v.Normal[0] /= length
	v.Normal[1] /= length
	v.Normal[2] /= length
}
