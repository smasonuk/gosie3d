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
