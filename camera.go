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
	c.camMatrixRev = z.MultiplyBy(y)
	c.camMatrixRev = c.camMatrixRev.MultiplyBy(x)
	c.cameraPosition = NewPoint3d(xp, yp, zp)
	c.cameraAngle = NewVector3(xa, ya, za)

	return c
}

func (c *Camera) GetCameraMatrix() *Matrix {
	return c.camMatrixRev
}

func NewCameraLookAt(camPos Vector3, lookAt Vector3, up Vector3) *Camera {
	lookAtMat := mgl64.LookAt(
		lookAt.GetX(), lookAt.GetY(), lookAt.GetZ(),
		camPos.GetX(), camPos.GetY(), camPos.GetZ(),

		0, 1, 0,
	)

	sMat := ToGoSieMatrix(lookAtMat)

	c := &Camera{
		camMatrixRev:   sMat,
		cameraPosition: NewPoint3d(camPos.GetX(), camPos.GetY(), camPos.GetZ()),
		cameraAngle:    NewVector3(0, 0, 0), //TODO: remove
	}

	return c
}

func NewCameraLookMatrixAt3(cameraLocation Point3d, lookAt Vector3, up Vector3) *Matrix {

	sTransWorldToCamera := TransMatrix(
		-cameraLocation.GetX(),
		-cameraLocation.GetY(),
		-cameraLocation.GetZ(),
	)

	dirY := lookAt.GetZ() - cameraLocation.GetZ()
	dirX := lookAt.GetX() - cameraLocation.GetX()

	// // find angle around the Y axis, and x
	angleY := math.Atan2(dirY, dirX)
	// angleX := math.Atan2(lookAt.GetY()-cameraLocation.GetY(), math.Sqrt(dirX*dirX+dirY*dirY))

	angleDegrees := angleY * (180 / math.Pi)
	angleDegrees = angleDegrees - 90 // Adjust to match the coordinate system
	angleRadY := angleDegrees * (math.Pi / 180)

	sMatY := NewRotationMatrix(ROTY, angleRadY)
	// sMatX := NewRotationMatrix(ROTX, -0.2)

	// sMat := sMatX.MultiplyBy(sMatY.MultiplyBy(sTransWorldToCamera))
	sMat := sMatY.MultiplyBy(sTransWorldToCamera)
	return sMat

	// lookAtMat := mgl64.LookAt(
	// 	c.GetX(),
	// 	c.GetY(),
	// 	c.GetZ(),
	// 	lookAt.GetX(),
	// 	lookAt.GetY(),
	// 	lookAt.GetZ(),

	// 	up.GetX(),
	// 	up.GetY(),
	// 	up.GetZ(),
	// )

	// return ToGoSieMatrix(lookAtMat)

	// pVec := mgl64.Vec3{up.GetX(), up.GetY(), up.GetZ()}
	// lookAtQuat := mgl64.QuatLookAtV(
	// 	mgl64.Vec3{c.GetX(), c.GetY(), c.GetZ()},
	// 	mgl64.Vec3{lookAt.GetX(), lookAt.GetY(), lookAt.GetZ()},
	// 	mgl64.Vec3{pVec[0], pVec[1], pVec[2]},
	// )
	// lookAtMat := lookAtQuat.Mat4()

	// return ToGoSieMatrix(lookAtMat)

}

func (c *Camera) LookAt(lookAt Vector3, up Vector3) {

	sMat := NewCameraLookMatrixAt3(*c.cameraPosition, lookAt, up)

	// .MultiplyBy(sMatX)

	// lookAtMat := mgl64.LookAt(

	// 	lookAt.GetX(),
	// 	lookAt.GetY(),
	// 	lookAt.GetZ(),
	// 	c.cameraPosition.GetX(),
	// 	c.cameraPosition.GetY(),
	// 	c.cameraPosition.GetZ(),
	// 	up.GetX(),
	// 	up.GetY(),
	// 	up.GetZ(),
	// )

	// // Create a quaternion for the up vector to ensure the camera's up direction is respected.
	// upVec := mgl64.Vec3{up.GetX(), up.GetY(), up.GetZ()}
	// lookAtQuat := mgl64.QuatLookAtV(

	// 	mgl64.Vec3{lookAt.GetX(), lookAt.GetY(), lookAt.GetZ()},
	// 	mgl64.Vec3{c.cameraPosition.GetX(), c.cameraPosition.GetY(), c.cameraPosition.GetZ()},

	// 	upVec,
	// )
	// lookAtMat := lookAtQuat.Mat4()

	// sMat := ToGoSieMatrix(lookAtMat)
	c.camMatrixRev = sMat
}

func (c *Camera) GetPosition() *Point3d {
	if c.cameraPosition == nil {
		return NewPoint3d(0, 0, 0)
	}
	return c.cameraPosition
}

func (c *Camera) SetCameraPosition(x, y, z float64) {
	c.cameraPosition = NewPoint3d(x, y, z)
}

func (c *Camera) AddXPosition(x float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.Points[0] += x
}

func (c *Camera) AddYPosition(y float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.Points[1] += y
}

func (c *Camera) AddZPosition(z float64) {
	if c.cameraPosition == nil {
		c.cameraPosition = NewPoint3d(0, 0, 0)
	}
	c.cameraPosition.Points[2] += z
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

	rotY := mgl64.HomogRotate3DY(-c.cameraAngle.GetY())
	rotX := mgl64.HomogRotate3DX(-c.cameraAngle.GetX())
	rotZ := mgl64.HomogRotate3DZ(-c.cameraAngle.GetZ())
	camRev := rotZ.Mul4(rotY).Mul4(rotX)

	c.camMatrixRev = ToGoSieMatrix(camRev)

}
