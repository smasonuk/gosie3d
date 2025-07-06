package main

import (
	"fmt"

	"github.com/go-gl/mathgl/mgl64"
	"github.com/smasonuk/gosie3d"
)

func main() {
	println("Application has started.")

	mat := gosie3d.NewRotationMatrix(gosie3d.ROTX, 1.0)
	fmt.Println("Rotation Matrix:")
	fmt.Println(mat)

	// q := mgl64.QuatLookAtV(mgl64.Vec3{1, 0, 0}, mgl64.Vec3{0, 0, 0}, mgl64.Vec3{0, 1, 0})
	// fmt.Println("Quaternion:", q)

	m := mgl64.HomogRotate3DX(1.0)
	fmt.Println("Homogeneous Rotation Matrix:")
	fmt.Println(m)

	sieMatrix := gosie3d.ToGoSieMatrix(m)
	fmt.Println("Converted to GoSie Matrix:")
	fmt.Println(sieMatrix)

}
