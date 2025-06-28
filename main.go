// To run this code:
// 1. Make sure you have Go installed.
// 2. Set up a new Go module: go mod init <your_module_name>
// 3. Get the Ebiten dependency: go get github.com/hajimehoshi/ebiten/v2
// 4. Save the code as main.go and run: go run .

package main

import (
	"image"
	"image/color"
	"log"
	"math"
	"math/rand"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

const (
	screenWidth  = 640
	screenHeight = 480
)

// =====================================================================================
// Vector3d (from java/Vector3d.java)
// =====================================================================================

type Vector3d struct {
	Normal [4]float64
}

func NewVector3d(x, y, z float64) *Vector3d {
	return &Vector3d{
		Normal: [4]float64{x, y, z, 1.0},
	}
}

func NewVector3dFromArray(normal []float64) *Vector3d {
	v := &Vector3d{}
	copy(v.Normal[:], normal)
	return v
}

func (v *Vector3d) Normalize() {
	length := math.Sqrt(math.Abs(v.Normal[0]*v.Normal[0] + v.Normal[1]*v.Normal[1] + v.Normal[2]*v.Normal[2]))
	if length == 0 {
		return
	}
	v.Normal[0] /= length
	v.Normal[1] /= length
	v.Normal[2] /= length
}

func GetLength(vec []float64) float64 {
	return math.Sqrt(math.Abs(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))
}

// =====================================================================================
// Point3d (from java/Point3d.java)
// =====================================================================================

type Point3d struct {
	Points []float64
}

func NewPoint3d(x, y, z float64) *Point3d {
	return &Point3d{
		Points: []float64{x, y, z, 1.0},
	}
}

func (p *Point3d) GetX() float64 { return p.Points[0] }
func (p *Point3d) GetY() float64 { return p.Points[1] }
func (p *Point3d) GetZ() float64 { return p.Points[2] }

// =====================================================================================
// Clist (from java/Clist.java)
// =====================================================================================

type Clist struct {
	points []*Point3d
	count  int
	max    int
	end    int
}

func NewClist(size int) *Clist {
	return &Clist{
		points: make([]*Point3d, size),
		max:    size,
	}
}
func (c *Clist) AddPoint(p *Point3d) {
	if c.end < c.max {
		c.points[c.end] = p
		c.end++
	}
}
func (c *Clist) Back() {
	c.count--
	if c.count < 0 {
		c.count = c.end - 1
	}
}
func (c *Clist) NextPoint() *Point3d {
	old := c.count
	c.count++
	if c.count >= c.end {
		c.count = 0
	}
	return c.points[old]
}

// =====================================================================================
// Face (from java/Face.java)
// =====================================================================================

type Face struct {
	Points  [][]float64
	Col     color.RGBA
	normal  []float64
	plane   *Plane
	vecPnts [][]float64
	meRev   bool
	Cnum    int
}

const (
	FACE_NORMAL  = 0
	FACE_REVERSE = 1
)

func NewFace(pnts [][]float64, col color.RGBA, normal []float64) *Face {
	f := &Face{
		Points: pnts,
		Col:    col,
		normal: normal,
	}
	if pnts != nil {
		f.Cnum = len(pnts)
	}
	return f
}

func (f *Face) AddPoint(x, y, z float64) {
	pnts := []float64{x, y, z, 1.0}
	f.vecPnts = append(f.vecPnts, pnts)
}

func (f *Face) GetNormal() []float64 {
	if f.normal == nil {
		f.createNormal()
	}
	return f.normal
}

func (f *Face) SetNormal(norm []float64) {
	f.normal = norm
}

func (f *Face) Finished(reverse int) {
	if reverse == FACE_REVERSE {
		f.meRev = true
	}
	f.Cnum = len(f.vecPnts)
	f.Points = make([][]float64, f.Cnum)
	copy(f.Points, f.vecPnts)
	f.vecPnts = nil
}

func (f *Face) GetPlane() *Plane {
	if f.plane != nil {
		return f.plane
	}
	f.plane = NewPlane(f, f.GetNormal())
	return f.plane
}

func (f *Face) createNormal() {
	if len(f.Points) < 3 {
		f.normal = []float64{0, 0, 1, 0}
		return
	}
	f.normal = make([]float64, 4)
	x1, y1, z1 := f.Points[0][0], f.Points[0][1], f.Points[0][2]
	x2, y2, z2 := f.Points[1][0], f.Points[1][1], f.Points[1][2]
	x3, y3, z3 := f.Points[2][0], f.Points[2][1], f.Points[2][2]

	u1, u2, u3 := x2-x1, y2-y1, z2-z1
	v1, v2, v3 := x3-x2, y3-y2, z3-z2

	if f.meRev {
		f.normal[0] = -(u2*v3 - u3*v2)
		f.normal[1] = -(u3*v1 - u1*v3)
		f.normal[2] = -(u1*v2 - u2*v1)
	} else {
		f.normal[0] = u2*v3 - u3*v2
		f.normal[1] = u3*v1 - u1*v3
		f.normal[2] = u1*v2 - u2*v1
	}

	nor := NewVector3dFromArray(f.normal)
	nor.Normalize()
	f.normal = nor.Normal[:]
}

// =====================================================================================
// Matrix (from java/Matrix.java)
// =====================================================================================

type Matrix struct {
	ThisMatrix [][]float64
}

const (
	ROTX = 0
	ROTY = 1
	ROTZ = 2
)

func NewMatrix() *Matrix {
	return &Matrix{
		ThisMatrix: make([][]float64, 0, 100),
	}
}

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
		// The original Java code implicitly used w=1 for multiplication.
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

// =====================================================================================
// Mesh (from java/Mesh.java)
// =====================================================================================

type Mesh struct {
	Points *Matrix
}

func NewMesh() *Mesh {
	return &Mesh{Points: NewMatrix()}
}

func (m *Mesh) AddPoint(point []float64) []float64 {
	ret := m.Points.FindPoints(point[0], point[1], point[2])
	if ret != nil {
		return ret
	}
	pointCopy := make([]float64, len(point))
	copy(pointCopy, point)
	m.Points.AddRow(pointCopy)
	return pointCopy
}

func (m *Mesh) Copy() *Mesh {
	return &Mesh{
		Points: m.Points.Copy(),
	}
}

type FaceMesh struct {
	Mesh
}

func NewFaceMesh() *FaceMesh {
	return &FaceMesh{Mesh: *NewMesh()}
}

func (fm *FaceMesh) AddFace(f *Face) *Face {
	newPoints := make([][]float64, len(f.Points))
	for i, p := range f.Points {
		newPoints[i] = fm.AddPoint(p)
	}
	return NewFace(newPoints, f.Col, f.GetNormal())
}

type NormalMesh struct {
	Mesh
}

func NewNormalMesh() *NormalMesh {
	return &NormalMesh{Mesh: *NewMesh()}
}

func (nm *NormalMesh) AddNormal(pnts []float64) []float64 {
	return nm.AddPoint(pnts)
}

// =====================================================================================
// FaceStore (from java/FaceStore.java)
// =====================================================================================

type FaceStore struct {
	faces []*Face
}

func NewFaceStore() *FaceStore {
	return &FaceStore{faces: make([]*Face, 0, 10)}
}
func (fs *FaceStore) AddFace(f *Face) {
	fs.faces = append(fs.faces, f)
}
func (fs *FaceStore) GetFace(i int) *Face {
	return fs.faces[i]
}
func (fs *FaceStore) FaceCount() int {
	return len(fs.faces)
}
func (fs *FaceStore) RemoveFaceAt(i int) *Face {
	if i < 0 || i >= len(fs.faces) {
		return nil
	}
	f := fs.faces[i]
	fs.faces = append(fs.faces[:i], fs.faces[i+1:]...)
	return f
}

// =====================================================================================
// Plane (from java/Plane.java)
// =====================================================================================

type Plane struct {
	A, B, C, D float64
}

const planeThickness = 0.1

func NewPlane(f *Face, normal []float64) *Plane {
	p := &Plane{
		A: normal[0],
		B: normal[1],
		C: normal[2],
	}
	p.D = -(p.A*f.Points[0][0] + p.B*f.Points[0][1] + p.C*f.Points[0][2])
	return p
}

func (p *Plane) PointOnPlane(x, y, z float64) float64 {
	num := p.A*x + p.B*y + p.C*z + p.D
	if math.Abs(num) > 0 && math.Abs(num) < planeThickness {
		return 0.0
	}
	return num
}

func (p *Plane) LIntersect(p1, p2 *Point3d) bool {
	a := p.PointOnPlane(p1.GetX(), p1.GetY(), p1.GetZ())
	b := p.PointOnPlane(p2.GetX(), p2.GetY(), p2.GetZ())
	if a == 0 || b == 0 {
		return false
	}
	return (a > 0 && b < 0) || (a < 0 && b > 0)
}

func (p *Plane) LineIntersect(p1, p2 *Point3d) *Point3d {
	x1, y1, z1 := p1.GetX(), p1.GetY(), p1.GetZ()
	x2, y2, z2 := p2.GetX(), p2.GetY(), p2.GetZ()

	if !p.LIntersect(p1, p2) {
		return nil
	}
	denom := (p.A*(x2-x1) + p.B*(y2-y1) + p.C*(z2-z1))
	if denom == 0 {
		return nil
	}
	t := -(p.A*x1 + p.B*y1 + p.C*z1 + p.D) / denom
	x := x1 + (x2-x1)*t
	y := y1 + (y2-y1)*t
	z := z1 + (z2-z1)*t
	return NewPoint3d(x, y, z)
}

func (p *Plane) FaceIntersect(f *Face) bool {
	var d float64
	initialized := false
	for a := 0; a < f.Cnum; a++ {
		n := p.PointOnPlane(f.Points[a][0], f.Points[a][1], f.Points[a][2])
		if !initialized {
			d = n
			initialized = true
			continue
		}
		if !((d >= 0 && n >= 0) || (d <= 0 && n <= 0)) {
			return true
		}
	}
	return false
}

func (p *Plane) SplitFace(aFace *Face) []*Face {
	faces := make([]*Face, 2)
	faces[0] = NewFace(nil, color.RGBA{}, nil)
	faces[1] = NewFace(nil, color.RGBA{}, nil)
	var inter bool

	if !p.FaceIntersect(aFace) {
		faces[0] = aFace
		faces[1] = nil
		return faces
	}

	currentFace := 0
	pnts := NewClist(aFace.Cnum)
	for i := 0; i < aFace.Cnum; i++ {
		pnts.AddPoint(NewPoint3d(aFace.Points[i][0], aFace.Points[i][1], aFace.Points[i][2]))
	}

	for pnt := 0; pnt < aFace.Cnum; pnt++ {
		p3d1 := pnts.NextPoint()
		p3d2 := pnts.NextPoint()
		pnts.Back()

		if p.LIntersect(p3d1, p3d2) {
			pointIntersect := p.LineIntersect(p3d1, p3d2)
			inter = true
			faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			if pointIntersect != nil {
				faces[currentFace].AddPoint(pointIntersect.GetX(), pointIntersect.GetY(), pointIntersect.GetZ())
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(pointIntersect.GetX(), pointIntersect.GetY(), pointIntersect.GetZ())
			}
		} else {
			if p.PointOnPlane(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ()) == 0 {
				inter = true
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
				currentFace = 1 - currentFace // flip
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			} else {
				faces[currentFace].AddPoint(p3d1.GetX(), p3d1.GetY(), p3d1.GetZ())
			}
		}
	}

	if !inter {
		faces[0] = aFace
		faces[1] = nil
		return faces
	}

	faces[0].Finished(FACE_NORMAL)
	faces[1].Finished(FACE_NORMAL)
	return faces
}

func (p *Plane) Where(f *Face) float64 {
	var inter float64
	for i := 0; i < len(f.Points); i++ {
		inter += p.PointOnPlane(f.Points[i][0], f.Points[i][1], f.Points[i][2])
	}
	return inter
}

// =====================================================================================
// BspNode (from java/BspNode.java)
// =====================================================================================

type BspNode struct {
	facePnts [][]float64
	normal   []float64
	Left     *BspNode
	Right    *BspNode
	colRed   uint8
	colGreen uint8
	colBlue  uint8
	xp       []float32
	yp       []float32
}

func NewBspNode(f *Face) *BspNode {
	b := &BspNode{
		facePnts: f.Points,
		normal:   f.GetNormal(),
		xp:       make([]float32, len(f.Points)),
		yp:       make([]float32, len(f.Points)),
	}
	b.colRed = f.Col.R
	b.colGreen = f.Col.G
	b.colBlue = f.Col.B
	return b
}

func (b *BspNode) Paint(screen *ebiten.Image, x, y int) {
	if len(b.facePnts) == 0 {
		return
	}
	where := b.normal[0]*b.facePnts[0][0] + b.normal[1]*b.facePnts[0][1] + b.normal[2]*b.facePnts[0][2]

	if where <= 0 {
		if b.Left != nil {
			b.Left.Paint(screen, x, y)
		}
		if b.Right != nil {
			b.Right.Paint(screen, x, y)
		}
	} else {
		if b.Right != nil {
			b.Right.Paint(screen, x, y)
		}

		for i := 0; i < len(b.facePnts); i++ {
			if b.facePnts[i][2] < 0.1 {
				if b.Left != nil {
					b.Left.Paint(screen, x, y)
				}
				return
			}
			b.xp[i] = float32((400*b.facePnts[i][0])/b.facePnts[i][2]) + float32(x)
			b.yp[i] = float32((400*b.facePnts[i][1])/b.facePnts[i][2]) + float32(y)
		}

		cosTheta := (b.normal[0]*b.facePnts[0][0] + b.normal[1]*b.facePnts[0][1] + b.normal[2]*b.facePnts[0][2]) / GetLength(b.facePnts[0])
		c := 240 - int(cosTheta*240)

		b1 := int(b.colBlue) - c
		if b1 < 0 {
			b1 = 0
		}
		if b1 > 255 {
			b1 = 255
		}
		r1 := int(b.colRed) - c
		if r1 < 0 {
			r1 = 0
		}
		if r1 > 255 {
			r1 = 255
		}
		g1 := int(b.colGreen) - c
		if g1 < 0 {
			g1 = 0
		}
		if g1 > 255 {
			g1 = 255
		}
		polyColor := color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}

		fillConvexPolygon(screen, b.xp, b.yp, polyColor)

		if b.Left != nil {
			b.Left.Paint(screen, x, y)
		}
	}
}

// =====================================================================================
// Object_3d (from java/Object_3d.java)
// =====================================================================================
type Object_3d struct {
	faceMesh        *FaceMesh
	normalMesh      *NormalMesh
	transFaceMesh   *FaceMesh
	transNormalMesh *NormalMesh
	theFaces        *FaceStore
	root            *BspNode
	rotMatrix       *Matrix
}

func NewObject_3d() *Object_3d {
	return &Object_3d{
		transFaceMesh:   NewFaceMesh(),
		transNormalMesh: NewNormalMesh(),
		theFaces:        NewFaceStore(),
		rotMatrix:       IdentMatrix(),
	}
}

func (o *Object_3d) Finished() {
	log.Println("Creating BSP Tree...")
	if o.theFaces.FaceCount() > 0 {
		o.root = o.createBspTree(o.theFaces, o.transFaceMesh, o.transNormalMesh)
	}
	log.Println("BSP Tree Created.")

	o.faceMesh = &FaceMesh{Mesh: *o.transFaceMesh.Mesh.Copy()}
	o.normalMesh = &NormalMesh{Mesh: *o.transNormalMesh.Mesh.Copy()}

	log.Printf("Points: %d", len(o.faceMesh.Points.ThisMatrix))
	log.Printf("Normals: %d", len(o.normalMesh.Points.ThisMatrix))
}

func (o *Object_3d) ApplyMatrixBatch(m *Matrix) {
	o.rotMatrix = m.MultiplyBy(o.rotMatrix)
}

func (o *Object_3d) ApplyMatrixTemp(aMatrix *Matrix) {
	rotMatrixTemp := aMatrix.MultiplyBy(o.rotMatrix)
	rotMatrixTemp.TransformObj(o.normalMesh.Points, o.transNormalMesh.Points)

	// FIX: Re-normalize all normals after transformation
	for _, n := range o.transNormalMesh.Points.ThisMatrix {
		v := NewVector3dFromArray(n)
		v.Normalize()
		copy(n, v.Normal[:])
	}

	rotMatrixTemp.TransformObj(o.faceMesh.Points, o.transFaceMesh.Points)
}

func (o *Object_3d) PaintSolid(screen *ebiten.Image, x, y int) {
	if o.root != nil {
		o.root.Paint(screen, x, y)
	}
}

func (o *Object_3d) createBspTree(faces *FaceStore, newFaces *FaceMesh, newNormMesh *NormalMesh) *BspNode {
	if faces.FaceCount() == 0 {
		return nil
	}
	parentFace := o.choosePlane(faces)
	parentFace.SetNormal(newNormMesh.AddNormal(parentFace.GetNormal()))
	newFace := newFaces.AddFace(parentFace)
	parent := NewBspNode(newFace)
	pPlane := NewPlane(newFace, newFace.GetNormal())

	fvLeft := NewFaceStore()
	fvRight := NewFaceStore()

	for a := 0; a < faces.FaceCount(); a++ {
		currentFace := faces.GetFace(a)
		if pPlane.FaceIntersect(currentFace) {
			split := pPlane.SplitFace(currentFace)
			if split == nil || len(split) == 0 {
				continue
			}
			if split[0] != nil && len(split[0].Points) > 0 {
				if pPlane.Where(split[0]) <= 0 {
					fvLeft.AddFace(NewFace(split[0].Points, currentFace.Col, currentFace.GetNormal()))
				} else {
					fvRight.AddFace(NewFace(split[0].Points, currentFace.Col, currentFace.GetNormal()))
				}
			}
			if len(split) > 1 && split[1] != nil && len(split[1].Points) > 0 {
				if pPlane.Where(split[1]) <= 0 {
					fvLeft.AddFace(NewFace(split[1].Points, currentFace.Col, currentFace.GetNormal()))
				} else {
					fvRight.AddFace(NewFace(split[1].Points, currentFace.Col, currentFace.GetNormal()))
				}
			}
		} else {
			w := pPlane.Where(currentFace)
			if w <= 0 {
				fvLeft.AddFace(currentFace)
			} else {
				fvRight.AddFace(currentFace)
			}
		}
	}
	if fvLeft.FaceCount() > 0 {
		parent.Left = o.createBspTree(fvLeft, newFaces, newNormMesh)
	}
	if fvRight.FaceCount() > 0 {
		parent.Right = o.createBspTree(fvRight, newFaces, newNormMesh)
	}
	return parent
}

func (o *Object_3d) choosePlane(fs *FaceStore) *Face {
	leastFace, leastFaceTotal := 0, fs.FaceCount()

	for chosen := 0; chosen < fs.FaceCount(); chosen++ {
		total := 0
		p := fs.GetFace(chosen).GetPlane()
		for i := 0; i < fs.FaceCount(); i++ {
			if i == chosen {
				continue
			}
			f := fs.GetFace(i)
			if p.FaceIntersect(f) {
				total++
			}
		}
		if total < leastFaceTotal {
			leastFaceTotal = total
			leastFace = chosen
			if total == 0 {
				break
			}
		}
	}
	return fs.RemoveFaceAt(leastFace)
}

func NewCube() *Object_3d {
	obj := NewObject_3d()
	s := 40.0 // size

	points := [][3]float64{
		{-s, -s, -s}, {s, -s, -s}, {s, s, -s}, {-s, s, -s}, // Front
		{-s, -s, s}, {s, -s, s}, {s, s, s}, {-s, s, s}, // Back
	}

	colors := []color.RGBA{
		{255, 0, 0, 255}, {0, 255, 0, 255}, {0, 0, 255, 255},
		{255, 255, 0, 255}, {0, 255, 255, 255}, {255, 0, 255, 255},
	}

	quads := [][]int{
		{0, 1, 2, 3}, // Front
		{1, 5, 6, 2}, // Right
		{5, 4, 7, 6}, // Back
		{4, 0, 3, 7}, // Left
		{3, 2, 6, 7}, // Top
		{4, 5, 1, 0}, // Bottom
	}

	for i, q := range quads {
		face := NewFace(nil, colors[i], nil)
		face.AddPoint(points[q[0]][0], points[q[0]][1], points[q[0]][2])
		face.AddPoint(points[q[1]][0], points[q[1]][1], points[q[1]][2])
		face.AddPoint(points[q[2]][0], points[q[2]][1], points[q[2]][2])
		face.AddPoint(points[q[3]][0], points[q[3]][1], points[q[3]][2])
		face.Finished(FACE_NORMAL)
		obj.theFaces.AddFace(face)
	}

	obj.Finished()
	return obj
}

// =====================================================================================
// Camera (from java/Camera.java)
// =====================================================================================

type Camera struct {
	CamMatrixRev *Matrix
}

func NewCamera(xa, ya, za float64) *Camera {
	c := &Camera{}
	x := NewRotationMatrix(ROTX, -xa)
	y := NewRotationMatrix(ROTY, -ya)
	z := NewRotationMatrix(ROTZ, -za)
	c.CamMatrixRev = z.MultiplyBy(y)
	c.CamMatrixRev = c.CamMatrixRev.MultiplyBy(x)
	return c
}

func (c *Camera) AddAngle(x, y float64) {
	rotY := NewRotationMatrix(ROTY, -y)
	rotX := NewRotationMatrix(ROTX, -x)
	c.CamMatrixRev = rotY.MultiplyBy(rotX).MultiplyBy(c.CamMatrixRev)
}

// =====================================================================================
// World_3d (from java/World_3d.java)
// =====================================================================================

type World_3d struct {
	objects       []*Object_3d
	objXpos       []float64
	objYpos       []float64
	objZpos       []float64
	cameras       []*Camera
	camXpos       []float64
	camYpos       []float64
	camZpos       []float64
	currentCamera int
}

func NewWorld_3d() *World_3d {
	return &World_3d{
		currentCamera: -1,
	}
}

func (w *World_3d) AddObject(obj *Object_3d, x, y, z float64) {
	w.objects = append(w.objects, obj)
	w.objXpos = append(w.objXpos, x)
	w.objYpos = append(w.objYpos, y)
	w.objZpos = append(w.objZpos, z)
}

func (w *World_3d) AddCamera(c *Camera, x, y, z float64) {
	w.cameras = append(w.cameras, c)
	w.camXpos = append(w.camXpos, x)
	w.camYpos = append(w.camYpos, y)
	w.camZpos = append(w.camZpos, z)
	w.currentCamera = len(w.cameras) - 1
}

func (w *World_3d) PaintObjects(screen *ebiten.Image, xsize, ysize int) {
	if w.currentCamera == -1 || len(w.cameras) == 0 {
		return
	}
	cam := w.cameras[w.currentCamera]
	camX, camY, camZ := w.camXpos[w.currentCamera], w.camYpos[w.currentCamera], w.camZpos[w.currentCamera]

	for i, obj := range w.objects {
		m := TransMatrix(w.objXpos[i]-camX, w.objYpos[i]-camY, w.objZpos[i]-camZ)
		obj.ApplyMatrixTemp(cam.CamMatrixRev.MultiplyBy(m))
		obj.PaintSolid(screen, xsize/2, ysize/2)
	}
}

// =====================================================================================
// Main Game struct and loop (replaces Start.java)
// =====================================================================================
var (
	whiteImage = ebiten.NewImage(3, 3)
)

func init() {
	whiteImage.Fill(color.White)
}

type Game struct {
	world        *World_3d
	i, p         float64
	cube         *Object_3d
	lastX, lastY int
	dragged      bool
}

func NewGame() *Game {
	g := &Game{}
	log.Println("Initializing World...")
	g.world = NewWorld_3d()
	theCamera := NewCamera(0, 0, 0)
	g.world.AddCamera(theCamera, 0, 0, 0)

	log.Println("Creating Cube...")
	g.cube = NewCube()
	g.world.AddObject(g.cube, 0, 0, 200)
	log.Println("Initialization Complete.")

	return g
}

func (g *Game) Update() error {
	s := NewForwardMatrix(math.Sin(g.i)/8/20, math.Cos(g.p)/8/20, math.Cos(g.i+1.14)/8/20)
	g.cube.ApplyMatrixBatch(s)
	g.i += 0.02
	g.p += 0.05

	x, y := ebiten.CursorPosition()
	if inpututil.IsMouseButtonJustPressed(ebiten.MouseButtonLeft) {
		if !g.dragged {
			g.dragged = true
			g.lastX, g.lastY = x, y
		}
	}
	if inpututil.IsMouseButtonJustReleased(ebiten.MouseButtonLeft) {
		g.dragged = false
	}
	if g.dragged {
		newX := float64(g.lastX - x)
		newX = -newX / 100.0
		newY := float64(g.lastY - y)
		newY /= 100.0

		cam := g.world.cameras[g.world.currentCamera]
		if cam != nil {
			cam.AddAngle(newY, newX)
		}
		g.lastX, g.lastY = x, y
	}
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	screen.Fill(color.Black)
	g.world.PaintObjects(screen, screenWidth, screenHeight)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return screenWidth, screenHeight
}

func main() {
	rand.Seed(time.Now().UnixNano())
	game := NewGame()
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Go 3D Engine (from Java)")
	if err := ebiten.RunGame(game); err != nil {
		log.Fatal(err)
	}
}

// =====================================================================================
// Ebiten Drawing Helpers
// =====================================================================================
func fillConvexPolygon(screen *ebiten.Image, xp, yp []float32, clr color.RGBA) {
	if len(xp) < 3 {
		return
	}

	indices := make([]uint16, 0, (len(xp)-2)*3)
	for i := 2; i < len(xp); i++ {
		indices = append(indices, 0, uint16(i-1), uint16(i))
	}

	vertices := make([]ebiten.Vertex, len(xp))
	// FIX: Use the uint8 RGBA fields directly and normalize to float32 [0,1]
	cr := float32(clr.R) / 255.0
	cg := float32(clr.G) / 255.0
	cb := float32(clr.B) / 255.0
	ca := float32(clr.A) / 255.0

	for i := range xp {
		vertices[i] = ebiten.Vertex{
			DstX:   xp[i],
			DstY:   yp[i],
			SrcX:   1,
			SrcY:   1,
			ColorR: cr,
			ColorG: cg,
			ColorB: cb,
			ColorA: ca,
		}
	}

	op := &ebiten.DrawTrianglesOptions{}
	op.FillRule = ebiten.FillAll
	screen.DrawTriangles(vertices, indices, whiteImage.SubImage(image.Rect(1, 1, 2, 2)).(*ebiten.Image), op)
}
