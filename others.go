package gosie3d

// Subtract returns a new Vector3 that is the difference of v1 and v2.
func Subtract(v1, v2 *Vector3) *Vector3 {
	return NewVector3(
		v1.X-v2.X,
		v1.Y-v2.Y,
		v1.Z-v2.Z,
	)
}

// Cross2 computes the Cross2 product of two vectors and returns a new Vector3.
func Cross2(v1, v2 *Vector3) *Vector3 {
	return NewVector3(
		v1.Y*v2.Z-v1.Z*v2.Y,
		v1.Z*v2.X-v1.X*v2.Z,
		v1.X*v2.Y-v1.Y*v2.X,
	)
}

// Dot computes the dot product of two vectors.
func Dot(v1, v2 *Vector3) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}
