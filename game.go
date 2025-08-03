package gosie3d

import (
	"image/color"
	"log"
)

type Game struct {
	world        *World
	i, p         float64
	cube         *Model
	lastX, lastY int
	dragged      bool
}

func NewGame() *Game {
	g := &Game{}
	log.Println("Initializing World...")
	g.world = NewWorld3d()
	theCamera := NewCamera(0, 0, 0, 0, 0, 0)
	g.world.AddCamera(theCamera, 0, 0, 0)

	log.Println("Creating Cube...")
	// g.cube = NewCube()
	// fileNAme := "sphere.dxf"
	// reader, err := os.Open(fileNAme)
	// if err != nil {
	// 	log.Fatalf("Error opening DXF file %s: %v", fileNAme, err)
	// }
	// defer reader.Close()

	// sp, err := NewObjectFromDXF(reader, 1)
	// if err != nil {
	// 	log.Fatalf("Error parsing DXF file %s: %v", fileNAme, err)
	// }
	// g.cube = sp

	// g.cube = NewRectangle(100, 100, 200, color.RGBA{R: 255, G: 0, B: 0, A: 255})
	// g.cube = NewSphere(100, 1, color.RGBA{R: 255, G: 0, B: 0, A: 255})
	// g.cube = NewSphereStriped(100, 2, color.RGBA{R: 255, G: 0, B: 0, A: 255}, color.RGBA{R: 0, G: 255, B: 0, A: 255}, 0.5)
	cube := NewUVSphere(100, 14, 8, color.RGBA{R: 255, G: 0, B: 0, A: 255}, color.RGBA{R: 0, G: 255, B: 0, A: 255}, 3)
	// cube2 := NewUVSphere(100, 14, 8, color.RGBA{R: 255, G: 0, B: 0, A: 255}, color.RGBA{R: 0, G: 255, B: 0, A: 255}, 3)

	// cubeClone := NewUVSphere(100, 14, 8, color.RGBA{R: 255, G: 0, B: 0, A: 255}, color.RGBA{R: 0, G: 255, B: 0, A: 255}, 3)
	cubeClone := cube.Clone()

	g.world.AddObject(cube, 0, 0, 1000)
	g.world.AddObject(cubeClone, 200, 0, 1000)

	log.Println("Initialization Complete.")

	return g
}
