package gosie3d

import "math"

func clamp(value, min, max int) int {
	if value < min {
		return min
	}
	if value > max {
		return max
	}
	return value
}

const (
	screenWidth  = 640
	screenHeight = 480
)

func GetLength(vec []float64) float64 {
	return math.Sqrt(math.Abs(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))
}
