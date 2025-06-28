package main

import (
	"fmt"
	"image/color"
	"log"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

const (
	screenWidth  = 640
	screenHeight = 480
	// Constants to mirror Java Matrix statics
	RotX = 0
	RotY = 1
	RotZ = 2
)

// --- Math and Core Types ---

// Vector3 represents a 3D vector or point.
type Vector3 [4]float32

// Matrix4 represents a 4x4 matrix for 3D transformations.
type Matrix4 [4][4]float32

func NewIdentMatrix() *Matrix4 {
	m := &Matrix4{}
	m[0][0] = 1.0
	m[1][1] = 1.0
	m[2][2] = 1.0
	m[3][3] = 1.0
	return m
}

func NewRotationMatrix(axis int, theta float64) *Matrix4 {
	m := &Matrix4{}
	c, s := float32(math.Cos(theta)), float32(math.Sin(theta))
	m[3][3] = 1.0
	switch axis {
	case RotX:
		m[0][0] = 1.0
		m[1][1] = c
		m[1][2] = s
		m[2][1] = -s
		m[2][2] = c
	case RotY:
		m[1][1] = 1.0
		m[0][0] = c
		m[0][2] = -s
		m[2][0] = s
		m[2][2] = c
	case RotZ:
		m[2][2] = 1.0
		m[0][0] = c
		m[0][1] = s
		m[1][0] = -s
		m[1][1] = c
	}
	return m
}

func NewTranslationMatrix(x, y, z float32) *Matrix4 {
	m := NewIdentMatrix()
	m[3][0] = x
	m[3][1] = y
	m[3][2] = z
	return m
}

func (m *Matrix4) MultiplyBy(other *Matrix4) *Matrix4 {
	result := &Matrix4{}
	for y := 0; y < 4; y++ {
		for x := 0; x < 4; x++ {
			result[x][y] = m[0][y]*other[x][0] + m[1][y]*other[x][1] + m[2][y]*other[x][2] + m[3][y]*other[x][3]
		}
	}
	return result
}

func (m *Matrix4) TransformPoints(src []Vector3) []Vector3 {
	dest := make([]Vector3, len(src))
	for i := range src {
		srcW := src[i][3]
		dest[i][0] = m[0][0]*src[i][0] + m[1][0]*src[i][1] + m[2][0]*src[i][2] + m[3][0]*srcW
		dest[i][1] = m[0][1]*src[i][0] + m[1][1]*src[i][1] + m[2][1]*src[i][2] + m[3][1]*srcW
		dest[i][2] = m[0][2]*src[i][0] + m[1][2]*src[i][1] + m[2][2]*src[i][2] + m[3][2]*srcW
		dest[i][3] = m[0][3]*src[i][0] + m[1][3]*src[i][1] + m[2][3]*src[i][2] + m[3][3]*srcW
	}
	return dest
}

// --- Geometry and Scene Structures ---

// Object3D represents a 3D model.
type Object3D struct {
	AllPoints   []Vector3 // Master list of unique vertices in local space
	TransPoints []Vector3 // Transformed vertices in camera space, updated each frame
	AllFaces    []*Face
	BspRoot     *BspNode
	RotMatrix   *Matrix4
}

// Face represents a single polygon, referencing vertices by index.
type Face struct {
	PointIndices []int
	Normal       Vector3
	Color        color.RGBA
	Owner        *Object3D // Back-reference to access vertex data
	cachedPlane  *Plane
}

// BspNode is a node in the BSP tree.
type BspNode struct {
	Face        *Face
	Left, Right *BspNode
}

// Plane represents a plane in 3D space (Ax + By + Cz + D = 0).
type Plane struct {
	A, B, C, D float32
}

func NewObject3D() *Object3D {
	return &Object3D{
		RotMatrix: NewIdentMatrix(),
	}
}

func (o *Object3D) Finished() {
	o.TransPoints = make([]Vector3, len(o.AllPoints))
	log.Printf("Building BSP Tree with %d faces...", len(o.AllFaces))
	o.BspRoot = o.createBspTree(o.AllFaces)
	log.Println("BSP Tree finished.")
}

func (f *Face) createNormal() {
	if len(f.PointIndices) < 3 {
		f.Normal = Vector3{0, 0, 1, 0}
		return
	}
	p1 := f.Owner.AllPoints[f.PointIndices[0]]
	p2 := f.Owner.AllPoints[f.PointIndices[1]]
	p3 := f.Owner.AllPoints[f.PointIndices[2]]

	u_x, u_y, u_z := p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]
	v_x, v_y, v_z := p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]

	nx, ny, nz := u_y*v_z-u_z*v_y, u_z*v_x-u_x*v_z, u_x*v_y-u_y*v_x
	l := float32(math.Sqrt(float64(nx*nx + ny*ny + nz*nz)))
	if l > 0 {
		f.Normal = Vector3{nx / l, ny / l, nz / l, 0}
	}
}

// flip normal
func (f *Face) FlipNormal() {
	if len(f.PointIndices) < 3 {
		return
	}
	f.Normal[0] = -f.Normal[0]
	f.Normal[1] = -f.Normal[1]
	f.Normal[2] = -f.Normal[2]
	if f.cachedPlane != nil {
		f.cachedPlane.A = -f.cachedPlane.A
		f.cachedPlane.B = -f.cachedPlane.B
		f.cachedPlane.C = -f.cachedPlane.C
		f.cachedPlane.D = -f.cachedPlane.D
	}
}

func (f *Face) GetPlane() *Plane {
	if f.cachedPlane == nil {
		p := f.Owner.AllPoints[f.PointIndices[0]]
		n := f.Normal
		d := -(n[0]*p[0] + n[1]*p[1] + n[2]*p[2])
		f.cachedPlane = &Plane{A: n[0], B: n[1], C: n[2], D: d}
	}
	return f.cachedPlane
}

const planeThickness float32 = 0.01

func (p *Plane) PointOnPlane(point Vector3) float32 {
	dist := p.A*point[0] + p.B*point[1] + p.C*point[2] + p.D
	if math.Abs(float64(dist)) < float64(planeThickness) {
		return 0
	}
	return dist
}

func (p *Plane) Where(face *Face) float32 {
	var totalDist float32
	for _, idx := range face.PointIndices {
		totalDist += p.PointOnPlane(face.Owner.AllPoints[idx])
	}
	return totalDist
}

func (o *Object3D) choosePlane(faces []*Face) (*Face, []*Face) {
	if len(faces) == 0 {
		return nil, nil
	}
	return faces[0], faces[1:]
}

func (o *Object3D) createBspTree(faces []*Face) *BspNode {
	if len(faces) == 0 {
		return nil
	}
	parentFace, remainingFaces := o.choosePlane(faces)
	parent := &BspNode{Face: parentFace}
	plane := parentFace.GetPlane()

	var leftFaces, rightFaces []*Face
	for _, currentFace := range remainingFaces {
		if plane.Where(currentFace) > 0 {
			rightFaces = append(rightFaces, currentFace)
		} else {
			leftFaces = append(leftFaces, currentFace)
		}
	}
	parent.Left = o.createBspTree(leftFaces)
	parent.Right = o.createBspTree(rightFaces)
	return parent
}

func (o *Object3D) ApplyMatrixBatch(m *Matrix4) {
	o.RotMatrix = m.MultiplyBy(o.RotMatrix)
}

func (o *Object3D) ApplyMatrixTemp(m *Matrix4) {
	fullMatrix := m.MultiplyBy(o.RotMatrix)
	o.TransPoints = fullMatrix.TransformPoints(o.AllPoints)
}

func (o *Object3D) PaintSolid(screen *ebiten.Image, x, y int) {
	if o.BspRoot == nil {
		return
	}
	vertices := make([]ebiten.Vertex, 0, len(o.AllFaces)*4)
	indices := make([]uint16, 0, len(o.AllFaces)*6)

	vertices, indices = o.BspRoot.Paint(screen, x, y, vertices, indices)

	if len(vertices) > 0 && len(indices) > 0 {
		opts := &ebiten.DrawTrianglesOptions{FillRule: ebiten.FillAll}
		screen.DrawTriangles(vertices, indices, whiteImage, opts)
	}
}

func (n *BspNode) Paint(screen *ebiten.Image, x, y int, vertices []ebiten.Vertex, indices []uint16) ([]ebiten.Vertex, []uint16) {
	if n == nil {
		return vertices, indices
	}

	p0 := n.Face.Owner.TransPoints[n.Face.PointIndices[0]]
	p1 := n.Face.Owner.TransPoints[n.Face.PointIndices[1]]
	p2 := n.Face.Owner.TransPoints[n.Face.PointIndices[2]]

	// Calculate the normal of the transformed face to check its direction relative to the camera.
	u_x, u_y, u_z := p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
	v_x, v_y, v_z := p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]
	normal_x, normal_y, normal_z := u_y*v_z-u_z*v_y, u_z*v_x-u_x*v_z, u_x*v_y-u_y*v_x

	// Test if the face is pointing towards the camera at (0,0,0).
	// A positive result means it's pointing away from the camera.
	dotProductWithOrigin := normal_x*p0[0] + normal_y*p0[1] + normal_z*p0[2]

	// The face is visible if its normal is pointing towards the camera.
	// In a right-handed system with camera looking down -Z, the normal of a visible face
	// should form an obtuse angle with the vector to a point on the face, making the dot product negative.
	if dotProductWithOrigin > 0 { // Front-facing
		vertices, indices = n.Right.Paint(screen, x, y, vertices, indices)
		vertices, indices = n.drawFace(screen, x, y, vertices, indices)
		vertices, indices = n.Left.Paint(screen, x, y, vertices, indices)
	} else { // Back-facing
		vertices, indices = n.Left.Paint(screen, x, y, vertices, indices)
		vertices, indices = n.Right.Paint(screen, x, y, vertices, indices)
	}
	return vertices, indices
}

func (n *BspNode) drawFace(screen *ebiten.Image, x, y int, vertices []ebiten.Vertex, indices []uint16) ([]ebiten.Vertex, []uint16) {
	facePnts := make([]Vector3, len(n.Face.PointIndices))
	for i, idx := range n.Face.PointIndices {
		facePnts[i] = n.Face.Owner.TransPoints[idx]
	}

	projected := make([][2]float32, len(facePnts))
	for i, p := range facePnts {
		if p[2] <= 1.0 { // Near-plane clipping
			return vertices, indices
		}
		projX := (400*p[0])/p[2] + float32(x)
		projY := -(400*p[1])/p[2] + float32(y) // Flipped Y for screen coordinates
		projected[i] = [2]float32{projX, projY}
	}

	// Simple flat lighting
	lightVec := Vector3{0.577, 0.577, -0.577, 0} // Normalized light vector
	dot := n.Face.Normal[0]*lightVec[0] + n.Face.Normal[1]*lightVec[1] + n.Face.Normal[2]*lightVec[2]
	if dot < 0 {
		dot = 0
	}
	intensity := 0.2 + 0.8*dot // Ambient + Diffuse

	r, g, b, _ := n.Face.Color.RGBA()
	colR, colG, colB := float32(r>>8)/255.0, float32(g>>8)/255.0, float32(b>>8)/255.0

	v := ebiten.Vertex{
		SrcX: 1, SrcY: 1,
		ColorR: colR * intensity, ColorG: colG * intensity, ColorB: colB * intensity, ColorA: 1,
	}

	if len(projected) < 3 {
		return vertices, indices
	}

	baseIndex := uint16(len(vertices))
	for i := 0; i < len(projected); i++ {
		v.DstX, v.DstY = projected[i][0], projected[i][1]
		vertices = append(vertices, v)
	}

	// Triangulate (assuming convex polygon)
	for i := 1; i < len(projected)-1; i++ {
		indices = append(indices, baseIndex, baseIndex+uint16(i), baseIndex+uint16(i+1))
	}
	return vertices, indices
}

// Camera represents the viewer in the 3D world.
type Camera struct {
	RevMatrix *Matrix4
}

func NewCamera(xa, ya, za float32) *Camera {
	c := &Camera{}
	x := NewRotationMatrix(RotX, -float64(xa))
	y := NewRotationMatrix(RotY, -float64(ya))
	z := NewRotationMatrix(RotZ, -float64(za))
	c.RevMatrix = z.MultiplyBy(y).MultiplyBy(x)
	return c
}

func (c *Camera) AddAngle(x, y float32) {
	rotX := NewRotationMatrix(RotX, float64(-x))
	rotY := NewRotationMatrix(RotY, float64(-y))
	c.RevMatrix = rotX.MultiplyBy(rotY).MultiplyBy(c.RevMatrix)
}

// World3D manages all objects and the camera in the scene.
type World3D struct {
	objects []*Object3D
	objPos  [][3]float32
	camera  *Camera
	camPos  [3]float32
}

func NewWorld3D() *World3D { return &World3D{} }

func (w *World3D) AddObject(obj *Object3D, x, y, z float32) {
	w.objects = append(w.objects, obj)
	w.objPos = append(w.objPos, [3]float32{x, y, z})
}

func (w *World3D) AddCamera(cam *Camera, x, y, z float32) {
	w.camera = cam
	w.camPos = [3]float32{x, y, z}
}

func (w *World3D) PaintObjects(screen *ebiten.Image, xsize, ysize int) {
	for i, obj := range w.objects {
		pos := w.objPos[i]
		transMatrix := NewTranslationMatrix(
			pos[0]-w.camPos[0],
			pos[1]-w.camPos[1],
			pos[2]-w.camPos[2],
		)
		finalMatrix := w.camera.RevMatrix.MultiplyBy(transMatrix)
		obj.ApplyMatrixTemp(finalMatrix)
		obj.PaintSolid(screen, xsize/2, ysize/2)
	}
}

// --- Game Logic and Ebiten Implementation ---

var whiteImage = ebiten.NewImage(3, 3)

func init() {
	whiteImage.Fill(color.White)
}

type Game struct {
	world        *World3D
	isDragging   bool
	lastX, lastY int
}

func createCube() *Object3D {
	obj := NewObject3D()
	s := float32(30.0) // Increased size for better visibility

	obj.AllPoints = []Vector3{
		{-s, -s, -s, 1}, {s, -s, -s, 1}, {s, s, -s, 1}, {-s, s, -s, 1}, // Z- face points (0-3)
		{-s, -s, s, 1}, {s, -s, s, 1}, {s, s, s, 1}, {-s, s, s, 1}, // Z+ face points (4-7)
	}

	facesData := []struct {
		indices []int
		color   color.RGBA
	}{
		{[]int{4, 5, 6, 7}, color.RGBA{0, 0, 255, 255}},   // Z+ face (blue)
		{[]int{1, 2, 3, 0}, color.RGBA{255, 0, 0, 255}},   // Z- face (red)
		{[]int{5, 1, 0, 4}, color.RGBA{0, 255, 0, 255}},   // Y- face (green)
		{[]int{7, 6, 2, 3}, color.RGBA{255, 255, 0, 255}}, // Y+ face (yellow)
		{[]int{5, 6, 2, 1}, color.RGBA{0, 255, 255, 255}}, // X+ face (cyan)
		{[]int{0, 3, 7, 4}, color.RGBA{255, 0, 255, 255}}, // X- face (magenta)
	}

	for _, fd := range facesData {
		face := &Face{
			PointIndices: fd.indices,
			Color:        fd.color,
			Owner:        obj,
		}
		face.createNormal()
		// face.FlipNormal() // Ensure normals face outward
		obj.AllFaces = append(obj.AllFaces, face)
	}

	obj.Finished()
	return obj
}

func NewGame() *Game {
	g := &Game{world: NewWorld3D()}
	model := createCube()

	// model := load
	g.world.AddObject(model, 0, 0, 200) // Pushed farther back
	g.world.AddCamera(NewCamera(0, 0, 0), 0, 0, 0)
	return g
}

func (g *Game) Update() error {
	// Automatic rotation
	rotX := NewRotationMatrix(RotX, 0.005)
	rotY := NewRotationMatrix(RotY, 0.009)
	g.world.objects[0].ApplyMatrixBatch(rotY.MultiplyBy(rotX))

	// Mouse camera control
	if inpututil.IsMouseButtonJustPressed(ebiten.MouseButtonLeft) {
		g.isDragging = true
		g.lastX, g.lastY = ebiten.CursorPosition()
	}
	if g.isDragging {
		x, y := ebiten.CursorPosition()
		dx := float32(x-g.lastX) / 200.0
		dy := float32(y-g.lastY) / 200.0
		g.world.camera.AddAngle(dy, dx)
		g.lastX, g.lastY = x, y
	}
	if inpututil.IsMouseButtonJustReleased(ebiten.MouseButtonLeft) {
		g.isDragging = false
	}
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	screen.Fill(color.Black)
	g.world.PaintObjects(screen, screenWidth, screenHeight)
	ebitenutil.DebugPrint(screen, fmt.Sprintf("FPS: %0.2f", ebiten.ActualFPS()))
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return screenWidth, screenHeight
}

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Go 3D Engine (Ebiten Conversion)")
	if err := ebiten.RunGame(NewGame()); err != nil {
		log.Fatal(err)
	}
}
