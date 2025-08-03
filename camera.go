package gosie3d

import (
	"math"

	"github.com/go-gl/mathgl/mgl64"
)

type Camera struct {
	camMatrixRev   *Matrix
	cameraPosition *Point3d
	cameraAngle    *Vector3
}

func NewCamera(xp, yp, zp, xa, ya, za float64) *Camera {
	c := &Camera{}
	x := NewRotationMatrix(ROTX, -xa)
	y := NewRotationMatrix(ROTY, -ya)
	z := NewRotationMatrix(ROTZ, -za)

	sTransWorldToCamera := TransMatrix(-xp, -yp, -zp)

	c.camMatrixRev = z.MultiplyBy(y)
	c.camMatrixRev = c.camMatrixRev.MultiplyBy(x)
	c.camMatrixRev = c.camMatrixRev.MultiplyBy(c.camMatrixRev.MultiplyBy(sTransWorldToCamera))

	// sMat := sMatX.MultiplyBy(sMatY.MultiplyBy(sTransWorldToCamera))

	c.cameraPosition = NewPoint3d(xp, yp, zp)
	c.cameraAngle = NewVector3(xa, ya, za)

	return c
}

func (c *Camera) GetCameraMatrix() *Matrix {
	return c.camMatrixRev
}

func NewCameraLookAt(camPos Vector3, lookAt Vector3, up Vector3) *Camera {
	lookAtMat := mgl64.LookAt(
		lookAt.X, lookAt.Y, lookAt.Z,
		camPos.X, camPos.Y, camPos.Z,

		0, 1, 0,
	)

	sMat := ToGoSieMatrix(lookAtMat)

	c := &Camera{
		camMatrixRev:   sMat,
		cameraPosition: NewPoint3d(camPos.X, camPos.Y, camPos.Z),
		cameraAngle:    NewVector3(0, 0, 0), //TODO: remove
	}

	return c
}

func NewCameraLookMatrixAt3(cameraLocation Point3d, lookAt Vector3, up Vector3) *Matrix {
	sTransWorldToCamera := TransMatrix(
		-cameraLocation.X,
		-cameraLocation.Y,
		-cameraLocation.Z,
	)

	dirY := lookAt.Y - cameraLocation.Y
	dirX := lookAt.X - cameraLocation.X
	dirZ := lookAt.Z - cameraLocation.Z
	lookADirVec := NewVector3(dirX, dirY, dirZ)

	angleRadY := angleY(lookAt, cameraLocation)
	sMatY := NewRotationMatrix(ROTY, angleRadY)
	newLookAtDirVec := sMatY.RotateVector3(lookADirVec)

	angleDown := angleDown(newLookAtDirVec)

	sMatX := NewRotationMatrix(ROTX, angleDown)

	sMat := sMatX.MultiplyBy(sMatY.MultiplyBy(sTransWorldToCamera))

	return sMat
}

func angleDown(lookADirVec *Vector3) float64 {
	hypot := math.Sqrt(lookADirVec.X*lookADirVec.X + lookADirVec.Y*lookADirVec.Y + lookADirVec.Z*lookADirVec.Z)
	adjacent := lookADirVec.Y

	angleDownRad := math.Acos(adjacent / hypot)        // Angle in radians
	angleDownDegrees := angleDownRad * (180 / math.Pi) // Convert to degrees
	angleDownDegrees = 90 - angleDownDegrees

	return degreesToRadians(angleDownDegrees)

}

func angleY(lookAt Vector3, cameraLocation Point3d) float64 {
	dirY := lookAt.Z - cameraLocation.Z
	dirX := lookAt.X - cameraLocation.X
	angleY := math.Atan2(dirY, dirX)
	angleY = angleY - math.Pi/2 // Adjust to match the camera's forward direction
	return angleY
}

func degreesToRadians(degrees float64) float64 {
	return degrees * (math.Pi / 180)
}

func (c *Camera) LookAt(lookAt Vector3, up Vector3) {

	sMat := NewCameraLookMatrixAt3(*c.cameraPosition, lookAt, up)

	c.camMatrixRev = sMat
}

func (c *Camera) GetPosition() *Point3d {
	if c.cameraPosition == nil {
		return NewPoint3d(0, 0, 0)
	}
	return c.cameraPosition
}

// func (c *Camera) updatedCameraPositionMatrix() {
// 	sTransWorldToCamera := TransMatrix(
// 		-c.cameraPosition.X,
// 		-c.cameraPosition.Y,
// 		-c.cameraPosition.Z,
// 	)

// 	c.camMatrixRev = sTransWorldToCamera.MultiplyBy(c.camMatrixRev)
// }

func (c *Camera) SetCameraPosition(x, y, z float64) {
	c.cameraPosition = NewPoint3d(x, y, z)
}

func (c *Camera) AddXPosition(x float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.X += x
}

func (c *Camera) AddYPosition(y float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.Y += y
}

func (c *Camera) AddZPosition(z float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.Z += z
}

func (c *Camera) GetMatrix() *Matrix {
	return c.camMatrixRev
}

// func (c *Camera) SetMatrix(m *Matrix) {
// 	c.camMatrixRev = m
// }

func (c *Camera) AddAngle(x, y, z float64) {
	// rotY := NewRotationMatrix(ROTY, -y)
	// rotX := NewRotationMatrix(ROTX, -x)
	// c.camMatrixRev = rotY.MultiplyBy(rotX).MultiplyBy(c.camMatrixRev)

	c.cameraAngle.Add(x, y, z)

	rotY := mgl64.HomogRotate3DY(-c.cameraAngle.Y)
	rotX := mgl64.HomogRotate3DX(-c.cameraAngle.X)
	rotZ := mgl64.HomogRotate3DZ(-c.cameraAngle.Z)
	camRev := rotZ.Mul4(rotY).Mul4(rotX)

	c.camMatrixRev = ToGoSieMatrix(camRev)

}
