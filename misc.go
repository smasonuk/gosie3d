package gosie3d

import (
	"image/color"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

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

func DrawLine(screen *ebiten.Image, startX, startY, endX, endY int, col color.Color) {
	vector.StrokeLine(screen, float32(startX), float32(startY), float32(endX), float32(endY), 1, col, false)
}
