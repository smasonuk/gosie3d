package gosie3d

import "math"

type Vector2 struct {
	X float64
	Y float64
}

func NewVectorFromAngle(angle float64) Vector2 {
	return Vector2{
		X: math.Cos(angle),
		Y: math.Sin(angle),
	}
}

func (v Vector2) Normalize() Vector2 {
	magnitude := math.Sqrt(v.X*v.X + v.Y*v.Y)

	if magnitude == 0 {
		return Vector2{X: 0, Y: 0}
	}

	return Vector2{X: v.X / magnitude, Y: v.Y / magnitude}
}

func RotateVector2(v Vector2, angle float64) Vector2 {
	cosAngle := math.Cos(angle)
	sinAngle := math.Sin(angle)

	newX := v.X*cosAngle - v.Y*sinAngle
	newY := v.X*sinAngle + v.Y*cosAngle

	return Vector2{X: newX, Y: newY}
}

func VectorAngle2(v Vector2) float64 {
	return math.Atan2(v.Y, v.X)
}

func XYAngle2(x float64, y float64) float64 {
	return math.Atan2(y, x)
}

// AngleBetweenVectors2 returns the angle between two Vector2s in radians
func AngleBetweenVectors2(a, b Vector2) float64 {
	dot := a.X*b.X + a.Y*b.Y
	magA := math.Sqrt(a.X*a.X + a.Y*a.Y)
	magB := math.Sqrt(b.X*b.X + b.Y*b.Y)

	if magA == 0 || magB == 0 {
		return 0 // or panic/error if you prefer
	}

	cosTheta := dot / (magA * magB)
	// Clamp cosTheta to [-1, 1] to avoid NaN due to floating point error
	if cosTheta > 1 {
		cosTheta = 1
	} else if cosTheta < -1 {
		cosTheta = -1
	}

	return math.Acos(cosTheta)
}

// mult by scalar
func (v Vector2) Mult(scalar float64) Vector2 {
	return Vector2{
		X: v.X * scalar,
		Y: v.Y * scalar,
	}
}
