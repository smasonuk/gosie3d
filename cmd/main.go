package main

import (
	"fmt"
	"math"
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

	x1, y1 := 100.0, 100.0
	x2, y2 := 100.0, 150.0

	dirx := x2 - x1
	diry := y2 - y1

	fmt.Printf("dirx: %f, diry: %f\n", dirx, diry)

	// // find angle around the Y axis
	angleY := math.Atan2(
		diry,
		dirx)

	angleDegrees := angleY * (180 / math.Pi)
	fmt.Printf("Camera angle around Y axis: %f degrees", angleDegrees)

}
