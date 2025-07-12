package gosie3d

import (
	"fmt"
	"math"
	"strings"

	"github.com/go-gl/mathgl/mgl64"
)

func NewMatrix() *Matrix {
	return &Matrix{
		ThisMatrix: make([][]float64, 0, 100),
	}
}

type Matrix struct {
	ThisMatrix [][]float64
}

const (
	ROTX = 0
	ROTY = 1
	ROTZ = 2
)

func NewMatrixFromData(aMatrix [][]float64) *Matrix {
	m := &Matrix{
		ThisMatrix: make([][]float64, len(aMatrix)),
	}
	for i := range aMatrix {
		m.ThisMatrix[i] = make([]float64, len(aMatrix[i]))
		copy(m.ThisMatrix[i], aMatrix[i])
	}
	return m
}

func NewRotationMatrix(aRotation int, theta float64) *Matrix {
	m := make([][]float64, 4)
	for i := range m {
		m[i] = make([]float64, 4)
	}
	switch aRotation {
	case ROTX:
		m[0][0] = 1.0
		m[1][1] = math.Cos(theta)
		m[2][1] = -math.Sin(theta)
		m[1][2] = math.Sin(theta)
		m[2][2] = math.Cos(theta)
		m[3][3] = 1.0
	case ROTY:
		m[0][0] = math.Cos(theta)
		m[2][0] = math.Sin(theta)
		m[0][2] = -math.Sin(theta)
		m[2][2] = math.Cos(theta)
		m[1][1] = 1.0
		m[3][3] = 1.0
	case ROTZ:
		m[2][2] = 1.0
		m[3][3] = 1.0
		m[0][0] = math.Cos(theta)
		m[1][0] = -math.Sin(theta)
		m[0][1] = math.Sin(theta)
		m[1][1] = math.Cos(theta)
	}
	return &Matrix{ThisMatrix: m}
}

func NewForwardMatrix(x, y, z float64) *Matrix {
	m := make([][]float64, 4)
	for i := range m {
		m[i] = make([]float64, 4)
	}
	cz, cx, cy := math.Cos(z), math.Cos(x), math.Cos(y)
	sz, sy, sx := math.Sin(z), math.Sin(y), math.Sin(x)

	m[0][0] = cy * cz
	m[1][0] = cy * -sz
	m[2][0] = sy
	m[0][1] = -sx*-sy*cz + cx*sz
	m[1][1] = -sx*-sy*-sz + cx*cz
	m[2][1] = -sx * cy
	m[0][2] = cx*-sy*cz + sx*sz
	m[1][2] = cx*-sy*-sz + sx*cz
	m[2][2] = cx * cy
	m[3][3] = 1.0
	return &Matrix{ThisMatrix: m}
}

func IdentMatrix() *Matrix {
	m := make([][]float64, 4)
	for i := range m {
		m[i] = make([]float64, 4)
	}
	m[0][0], m[1][1], m[2][2], m[3][3] = 1.0, 1.0, 1.0, 1.0
	return &Matrix{ThisMatrix: m}
}

func TransMatrix(x, y, z float64) *Matrix {
	nm := make([][]float64, 4)
	for i := range nm {
		nm[i] = make([]float64, 4)
	}
	nm[3][0] = x
	nm[3][1] = y
	nm[3][2] = z
	nm[0][0], nm[1][1], nm[2][2], nm[3][3] = 1.0, 1.0, 1.0, 1.0
	return &Matrix{ThisMatrix: nm}
}

func (m *Matrix) AddRow(row []float64) {
	m.ThisMatrix = append(m.ThisMatrix, row)
}

func (m *Matrix) MultiplyBy(aMatrix *Matrix) *Matrix {
	newMatrixData := make([][]float64, len(aMatrix.ThisMatrix))
	for i := range newMatrixData {
		newMatrixData[i] = make([]float64, 4)
	}

	for y := 0; y < 4; y++ {
		for x := 0; x < len(aMatrix.ThisMatrix); x++ {
			newMatrixData[x][y] = m.ThisMatrix[0][y]*aMatrix.ThisMatrix[x][0] +
				m.ThisMatrix[1][y]*aMatrix.ThisMatrix[x][1] +
				m.ThisMatrix[2][y]*aMatrix.ThisMatrix[x][2] +
				m.ThisMatrix[3][y]*aMatrix.ThisMatrix[x][3]
		}
	}
	return &Matrix{ThisMatrix: newMatrixData}
}

func (m *Matrix) TransformObj(src, dest *Matrix) {
	for x := 0; x < len(src.ThisMatrix); x++ {
		sx, sy, sz := src.ThisMatrix[x][0], src.ThisMatrix[x][1], src.ThisMatrix[x][2]
		dest.ThisMatrix[x][0] = m.ThisMatrix[0][0]*sx + m.ThisMatrix[1][0]*sy + m.ThisMatrix[2][0]*sz + m.ThisMatrix[3][0]
		dest.ThisMatrix[x][1] = m.ThisMatrix[0][1]*sx + m.ThisMatrix[1][1]*sy + m.ThisMatrix[2][1]*sz + m.ThisMatrix[3][1]
		dest.ThisMatrix[x][2] = m.ThisMatrix[0][2]*sx + m.ThisMatrix[1][2]*sy + m.ThisMatrix[2][2]*sz + m.ThisMatrix[3][2]
	}
}

func (m *Matrix) FindPoints(x, y, z float64) []float64 {
	for _, p := range m.ThisMatrix {
		if p[0] == x && p[1] == y && p[2] == z {
			return p
		}
	}
	return nil
}

func (m *Matrix) Copy() *Matrix {
	return NewMatrixFromData(m.ThisMatrix)
}

func (m *Matrix) String() string {
	var sb strings.Builder
	for i, row := range m.ThisMatrix {
		if i > 0 {
			sb.WriteString("\n")
		}
		for j, val := range row {
			if j > 0 {
				sb.WriteString(" ")
			}
			sb.WriteString(fmt.Sprintf("%f", val))
		}
	}
	return sb.String()
}

// LookAtMatrix generates a complete view matrix that orients the world
// from the perspective of an eye position looking at a target.
func LookAtMatrix(eye, target, up *Vector3) *Matrix {
	// --- Calculate Camera's Local Axes (Right-Handed System) ---

	// 1. zAxis: The "forward" vector of the camera. Points from the target to the eye.
	zAxisVec := NewVector3(
		eye.X-target.X,
		eye.Y-target.Y,
		eye.Z-target.Z,
	)
	zAxisVec.Normalize()
	// zAxis := zAxisVec.Normal
	zAxis := zAxisVec.Copy()

	// 2. xAxis: The "right" vector of the camera.
	// xAxisVec := NewVector3dFromArray(Cross(up.Normal[:], zAxis[:]))
	xAxisVec := Cross(up, zAxis)
	xAxisVec.Normalize()
	// xAxis := xAxisVec.Normal
	xAxis := xAxisVec.Copy()

	// 3. yAxis: The "up" vector of the camera.
	// yAxis := Cross(zAxis[:], xAxis[:])
	yAxis := Cross(zAxis, xAxis)

	// --- Construct the Final View Matrix ---
	// This matrix combines the rotation and translation in the specific
	// format your engine's TransformObj function requires.
	viewMatrix := IdentMatrix()

	// The rotation part is built from the basis vectors in columns.
	// This correctly orients the world to the camera's view.
	// viewMatrix.ThisMatrix[0][0] = xAxis[0]
	// viewMatrix.ThisMatrix[1][0] = xAxis[1]
	// viewMatrix.ThisMatrix[2][0] = xAxis[2]
	viewMatrix.ThisMatrix[0][0] = xAxis.X
	viewMatrix.ThisMatrix[1][0] = xAxis.Y
	viewMatrix.ThisMatrix[2][0] = xAxis.Z

	// viewMatrix.ThisMatrix[0][1] = yAxis[0]
	// viewMatrix.ThisMatrix[1][1] = yAxis[1]
	// viewMatrix.ThisMatrix[2][1] = yAxis[2]
	viewMatrix.ThisMatrix[0][1] = yAxis.X
	viewMatrix.ThisMatrix[1][1] = yAxis.Y
	viewMatrix.ThisMatrix[2][1] = yAxis.Z

	// viewMatrix.ThisMatrix[0][2] = zAxis[0]
	// viewMatrix.ThisMatrix[1][2] = zAxis[1]
	// viewMatrix.ThisMatrix[2][2] = zAxis[2]
	viewMatrix.ThisMatrix[0][2] = zAxis.X
	viewMatrix.ThisMatrix[1][2] = zAxis.Y
	viewMatrix.ThisMatrix[2][2] = zAxis.Z

	// The translation part moves the entire world so the camera is at the origin.
	// It is calculated by taking the dot product of each axis with the eye's position.
	// viewMatrix.ThisMatrix[3][0] = -(xAxis[0]*eye.X + xAxis[1]*eye.Y + xAxis[2]*eye.Z)
	// viewMatrix.ThisMatrix[3][1] = -(yAxis[0]*eye.X + yAxis[1]*eye.Y + yAxis[2]*eye.Z)
	// viewMatrix.ThisMatrix[3][2] = -(zAxis[0]*eye.X + zAxis[1]*eye.Y + zAxis[2]*eye.Z)
	viewMatrix.ThisMatrix[3][0] = -(xAxis.X*eye.X + xAxis.Y*eye.Y + xAxis.Z*eye.Z)
	viewMatrix.ThisMatrix[3][1] = -(yAxis.X*eye.X + yAxis.Y*eye.Y + yAxis.Z*eye.Z)
	viewMatrix.ThisMatrix[3][2] = -(zAxis.X*eye.X + zAxis.Y*eye.Y + zAxis.Z*eye.Z)

	return viewMatrix
}

// Cross calculates the cross product of two 3-element vectors.
func Cross(a, b *Vector3) *Vector3 {
	// return []float64{
	// 	// a[1]*b[2] - a[2]*b[1],
	// 	// a[2]*b[0] - a[0]*b[2],
	// 	// a[0]*b[1] - a[1]*b[0],
	// 	// 0, // W component is 0 for vectors
	// 	a.Y*b.Z - a.Z*b.Y,
	// 	a.Z*b.X - a.X*b.Z,
	// 	a.X*b.Y - a.Y*b.X,
	// 	0, // W component is 0 for vectors
	// }
	return NewVector3(
		a.Y*b.Z-a.Z*b.Y,
		a.Z*b.X-a.X*b.Z,
		a.X*b.Y-a.Y*b.X,
	)
}

func ToGoSieMatrix(m mgl64.Mat4) *Matrix {
	return NewMatrixFromData(
		[][]float64{
			{m[0], m[1], m[2], m[3]},
			{m[4], m[5], m[6], m[7]},
			{m[8], m[9], m[10], m[11]},
			{m[12], m[13], m[14], m[15]},
		},
	)
}

// // rotate vector3 by this matrix
//
//	func (m *Matrix) RotateVector3(v *Vector3) *Vector3 {
//		// x, y, z := v.Normal[0], v.Normal[1], v.Normal[2]
//		x, y, z := v.X, v.Y, v.Z
//		return NewVector3(
//			m.ThisMatrix[0][0]*x+m.ThisMatrix[1][0]*y+m.ThisMatrix[2][0]*z,
//			m.ThisMatrix[0][1]*x+m.ThisMatrix[1][1]*y+m.ThisMatrix[2][1]*z,
//			m.ThisMatrix[0][2]*x+m.ThisMatrix[1][2]*y+m.ThisMatrix[2][2]*z,
//		)
//	}
//

// RotateVector3 rotates a Vector3 by the matrix's 3x3 rotation component.
// It does not apply translation, making it suitable for direction vectors.
func (m *Matrix) RotateVector3(v *Vector3) *Vector3 {
	// Extract the source vector components for clarity.
	vx, vy, vz := v.X, v.Y, v.Z

	// Apply the 3x3 rotation part of the matrix. This is a standard
	// vector-matrix multiplication that ignores the translation part of the matrix.
	newX := m.ThisMatrix[0][0]*vx + m.ThisMatrix[1][0]*vy + m.ThisMatrix[2][0]*vz
	newY := m.ThisMatrix[0][1]*vx + m.ThisMatrix[1][1]*vy + m.ThisMatrix[2][1]*vz
	newZ := m.ThisMatrix[0][2]*vx + m.ThisMatrix[1][2]*vy + m.ThisMatrix[2][2]*vz

	// Return a new Vector3 with the rotated coordinates.
	return NewVector3(newX, newY, newZ)
}

func (m *Matrix) TransformNormals(src, dest *Matrix) {
	for x := 0; x < len(src.ThisMatrix); x++ {
		sx, sy, sz := src.ThisMatrix[x][0], src.ThisMatrix[x][1], src.ThisMatrix[x][2]

		// This is a 3x3 rotation of a vector, it deliberately ignores the
		// translation components of the matrix (m.ThisMatrix[3][...]).
		dest.ThisMatrix[x][0] = m.ThisMatrix[0][0]*sx + m.ThisMatrix[1][0]*sy + m.ThisMatrix[2][0]*sz
		dest.ThisMatrix[x][1] = m.ThisMatrix[0][1]*sx + m.ThisMatrix[1][1]*sy + m.ThisMatrix[2][1]*sz
		dest.ThisMatrix[x][2] = m.ThisMatrix[0][2]*sx + m.ThisMatrix[1][2]*sy + m.ThisMatrix[2][2]*sz
	}
}
