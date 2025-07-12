package main

import (
	"fmt"
	"math"

	"github.com/smasonuk/gosie3d"
)

func main() {
	// println("Application has started.")

	// mat := gosie3d.NewRotationMatrix(gosie3d.ROTX, 1.0)
	// matz := gosie3d.NewRotationMatrix(gosie3d.ROTZ, 1.0)

	// mat = mat.MultiplyBy(matz)
	// fmt.Println("Rotation Matrix:")
	// fmt.Println(mat)

	// // q := mgl64.QuatLookAtV(mgl64.Vec3{1, 0, 0}, mgl64.Vec3{0, 0, 0}, mgl64.Vec3{0, 1, 0})
	// // fmt.Println("Quaternion:", q)

	// m := mgl64.HomogRotate3DX(1.0)
	// mz := mgl64.HomogRotate3DZ(1.0)
	// m = m.Mul4(mz)
	// fmt.Println("Homogeneous Rotation Matrix:")
	// fmt.Println(m)

	// sieMatrix := gosie3d.ToGoSieMatrix(m)
	// fmt.Println("Converted to GoSie Matrix:")
	// fmt.Println(sieMatrix)

	// 	x1, y1 := 100.0, 100.0
	// 	x2, y2 := 100.0, 150.0

	// 	dirx := x2 - x1
	// 	diry := y2 - y1

	// 	fmt.Printf("dirx: %f, diry: %f\n", dirx, diry)

	// 	// // find angle around the Y axis
	// 	angleY := math.Atan2(
	// 		diry,
	// 		dirx)

	// 	angleDegrees := angleY * (180 / math.Pi)
	// 	fmt.Printf("Camera angle around Y axis: %f degrees", angleDegrees)

	lookAt := gosie3d.NewVector3(20, 0, 20)
	cameraLocation := gosie3d.NewPoint3d(20, -20, 10)

	angleRadY := angleY(lookAt, cameraLocation)
	angleDegreesY := angleRadY * (180 / math.Pi)

	fmt.Println("angledeg:", angleDegreesY)

	dirY := lookAt.Z - cameraLocation.Z
	dirX := lookAt.X - cameraLocation.X
	dirZ := lookAt.Y - cameraLocation.Y
	lookADirVec := gosie3d.NewVector3(dirY, dirX, dirZ)

	sMatY := gosie3d.NewRotationMatrix(gosie3d.ROTY, angleRadY)
	newLookAtDirVec := sMatY.RotateVector3(lookADirVec)
	fmt.Printf("angleDegreesY: %v\n", angleDegreesY)
	fmt.Printf("lookAt: %v\n", lookADirVec)
	fmt.Printf("newLookAtDirVec: %v\n", newLookAtDirVec)

	down := angleDown(newLookAtDirVec)
	fmt.Printf("angle down: %f degrees\n", down)

}

func angleY(lookAt *gosie3d.Vector3, cameraLocation *gosie3d.Point3d) float64 {
	dirY := lookAt.Z - cameraLocation.Z
	dirX := lookAt.X - cameraLocation.X
	angleY := math.Atan2(dirY, dirX)
	angleY = angleY - math.Pi/2 // Adjust to match the camera's forward direction
	return angleY
}

func angleDown(lookADirVec *gosie3d.Vector3) float64 {
	hypot := math.Sqrt(lookADirVec.X*lookADirVec.X + lookADirVec.Y*lookADirVec.Y + lookADirVec.Z*lookADirVec.Z)
	adjacent := lookADirVec.Z

	angleDownRad := math.Acos(adjacent / hypot)        // Angle in radians
	angleDownDegrees := angleDownRad * (180 / math.Pi) // Convert to degrees
	angleDownDegrees = 90 - angleDownDegrees           // Adjust to match the camera's downward angle
	fmt.Printf("Angle down: %f degrees\n", angleDownDegrees)

	return degreesToRadians(angleDownDegrees)

}

func degreesToRadians(degrees float64) float64 {
	return degrees * (math.Pi / 180)
}
