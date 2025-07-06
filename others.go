package gosie3d

// --- Vector Math Helpers (using your Vector3 struct) ---

// Subtract returns a new Vector3 that is the difference of v1 and v2.
func Subtract(v1, v2 *Vector3) *Vector3 {
	return NewVector3(
		v1.GetX()-v2.GetX(),
		v1.GetY()-v2.GetY(),
		v1.GetZ()-v2.GetZ(),
	)
}

// Cross2 computes the Cross2 product of two vectors and returns a new Vector3.
func Cross2(v1, v2 *Vector3) *Vector3 {
	return NewVector3(
		v1.GetY()*v2.GetZ()-v1.GetZ()*v2.GetY(),
		v1.GetZ()*v2.GetX()-v1.GetX()*v2.GetZ(),
		v1.GetX()*v2.GetY()-v1.GetY()*v2.GetX(),
	)
}

// Dot computes the dot product of two vectors.
func Dot(v1, v2 *Vector3) float64 {
	return v1.GetX()*v2.GetX() + v1.GetY()*v2.GetY() + v1.GetZ()*v2.GetZ()
}

// NewLookAtMatrix creates a view matrix for a camera using the provided Vector3 struct.
// It uses a right-handed coordinate system.
//
// position: The 3D world position of the camera (*Vector3).
// target: The 3D world position the camera is looking at (*Vector3).
// up: The world's up direction, typically (0, 1, 0) (*Vector3).
// func NewLookAtMatrix(position, target, up *Vector3) *Matrix {
// 	// 1. Calculate the camera's new basis vectors.

// 	// CORRECTED: The "forward" vector must point from the camera TO the target.
// 	zAxis := Subtract(target, position)
// 	zAxis.Normalize()

// 	// The "right" vector is perpendicular to the world's up vector and the new z-axis.
// 	xAxis := Cross2(up, zAxis)
// 	xAxis.Normalize()

// 	// The camera's "up" vector is perpendicular to the other two axes.
// 	yAxis := Cross2(zAxis, xAxis)

// 	// 2. Create the 4x4 view matrix.
// 	m := make([][]float64, 4)
// 	for i := range m {
// 		m[i] = make([]float64, 4)
// 	}

// 	// Set the rotation part (transposed basis vectors)
// 	m[0][0] = xAxis.GetX()
// 	m[0][1] = yAxis.GetX()
// 	m[0][2] = zAxis.GetX()

// 	m[1][0] = xAxis.GetY()
// 	m[1][1] = yAxis.GetY()
// 	m[1][2] = zAxis.GetY()

// 	m[2][0] = xAxis.GetZ()
// 	m[2][1] = yAxis.GetZ()
// 	m[2][2] = zAxis.GetZ()

// 	// Set the translation part using the dot product
// 	m[3][0] = -Dot(xAxis, position)
// 	m[3][1] = -Dot(yAxis, position)
// 	m[3][2] = -Dot(zAxis, position)
// 	m[3][3] = 1.0

// 	return &Matrix{ThisMatrix: m}
// }

// NewCameraLookAt creates a camera at a position looking at a target.
// This is a robust, one-step implementation that builds the complete view matrix.
func NewCameraLookAt2(x, y, z, lookX, lookY, lookZ float64) *Camera {
	eye := NewVector3(x, y, z)
	target := NewVector3(lookX, lookY, lookZ)
	up := NewVector3(0, 1, 0)

	// --- Calculate Camera's Local Axes (Right-Handed System) ---
	// 1. zAxis: The "forward" direction of the camera.
	zAxisVec := Subtract(target, eye)
	zAxisVec.Normalize()

	// 2. xAxis: The "right" vector, perpendicular to world up and forward.
	xAxisVec := Cross2(up, zAxisVec)
	xAxisVec.Normalize()

	// 3. yAxis: The camera's "up" vector, perpendicular to its right and forward axes.
	yAxisVec := Cross2(zAxisVec, xAxisVec)

	// --- Construct the Final View Matrix (Row-Major) ---
	viewMatrix := IdentMatrix()
	m := viewMatrix.ThisMatrix

	// Set the rotation part with basis vectors as ROWS.
	m[0][0] = xAxisVec.GetX()
	m[0][1] = xAxisVec.GetY()
	m[0][2] = xAxisVec.GetZ()

	m[1][0] = yAxisVec.GetX()
	m[1][1] = yAxisVec.GetY()
	m[1][2] = yAxisVec.GetZ()

	m[2][0] = zAxisVec.GetX()
	m[2][1] = zAxisVec.GetY()
	m[2][2] = zAxisVec.GetZ()

	// Set the translation part in the fourth row.
	m[3][0] = -Dot(xAxisVec, eye)
	m[3][1] = -Dot(yAxisVec, eye)
	m[3][2] = -Dot(zAxisVec, eye)
	m[3][3] = 1.0

	// Create and return the camera
	return &Camera{
		camMatrixRev:   viewMatrix,
		cameraPosition: NewPoint3d(x, y, z),
	}
}

func NewLookAtMatrix2(eye, target, up *Vector3) *Matrix {
	// --- Calculate Camera's Local Axes (Right-Handed System) ---

	// 1. zAxis: The "forward" vector. CORRECTED to be (target - eye).
	zAxisVec := Subtract(target, eye)
	zAxisVec.Normalize()

	// 2. xAxis: The "right" vector.
	xAxisVec := Cross2(up, zAxisVec)
	xAxisVec.Normalize()

	// 3. yAxis: The camera's "up" vector.
	yAxisVec := Cross2(zAxisVec, xAxisVec)
	// No need to normalize yAxisVec, as the other two are orthogonal unit vectors.

	// --- Construct the Final View Matrix ---
	viewMatrix := IdentMatrix()
	m := viewMatrix.ThisMatrix

	// CORRECTED: The rotation part is now built from ROWS, not columns.
	// This matches how your TransformObj function works.
	m[0][0] = xAxisVec.GetX()
	m[0][1] = xAxisVec.GetY()
	m[0][2] = xAxisVec.GetZ()

	m[1][0] = yAxisVec.GetX()
	m[1][1] = yAxisVec.GetY()
	m[1][2] = yAxisVec.GetZ()

	m[2][0] = zAxisVec.GetX()
	m[2][1] = zAxisVec.GetY()
	m[2][2] = zAxisVec.GetZ()

	// The translation part moves the world relative to the camera's new orientation.
	// This calculation remains the same.
	m[3][0] = -Dot(xAxisVec, eye)
	m[3][1] = -Dot(yAxisVec, eye)
	m[3][2] = -Dot(zAxisVec, eye)

	return viewMatrix
}
