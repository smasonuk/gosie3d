package gosie3d

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

const (
	screenWidth  = 640
	screenHeight = 480
)

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

// Quaternion represents a rotation in 3D space.
type Quaternion struct {
	W, X, Y, Z float64
}

// NewQuaternion creates a new quaternion from its raw components.
func NewQuaternion(w, x, y, z float64) *Quaternion {
	return &Quaternion{W: w, X: x, Y: y, Z: z}
}

// NewQuaternionFromAxisAngle creates a quaternion representing a rotation
// of `angle` radians around the given `axis`. The axis vector must be normalized.
func NewQuaternionFromAxisAngle(axis [3]float64, angle float64) *Quaternion {
	halfAngle := angle / 2.0
	sinHalfAngle := math.Sin(halfAngle)
	return &Quaternion{
		W: math.Cos(halfAngle),
		X: axis[0] * sinHalfAngle,
		Y: axis[1] * sinHalfAngle,
		Z: axis[2] * sinHalfAngle,
	}
}

// Normalize ensures the quaternion is a unit quaternion, which is necessary
// for it to represent a pure rotation.
func (q *Quaternion) Normalize() {
	length := math.Sqrt(q.W*q.W + q.X*q.X + q.Y*q.Y + q.Z*q.Z)
	if length == 0 {
		return
	}
	q.W /= length
	q.X /= length
	q.Y /= length
	q.Z /= length
}

// ToRotationMatrix converts the quaternion into a 4x4 rotation matrix.
func (q *Quaternion) ToRotationMatrix() *Matrix {
	q.Normalize() // Ensure the quaternion is a unit quaternion

	xx, yy, zz := q.X*q.X, q.Y*q.Y, q.Z*q.Z
	xy, xz, yz := q.X*q.Y, q.X*q.Z, q.Y*q.Z
	wx, wy, wz := q.W*q.X, q.W*q.Y, q.W*q.Z

	m := make([][]float64, 4)
	for i := range m {
		m[i] = make([]float64, 4)
	}

	m[0][0] = 1 - 2*(yy+zz)
	m[0][1] = 2 * (xy - wz)
	m[0][2] = 2 * (xz + wy)

	m[1][0] = 2 * (xy + wz)
	m[1][1] = 1 - 2*(xx+zz)
	m[1][2] = 2 * (yz - wx)

	m[2][0] = 2 * (xz - wy)
	m[2][1] = 2 * (yz + wx)
	m[2][2] = 1 - 2*(xx+yy)

	m[3][3] = 1.0

	return &Matrix{ThisMatrix: m}
}

// LookAtMatrix generates a Matrix that orients an object at `eye` position
// to look at a `target` position. The `up` vector defines the upward direction
// and helps stabilize the rotation.
func LookAtMatrix(eye, target, up *Vector3d) *Matrix {
	// 1. Calculate the forward vector (from target to eye) and normalize it.
	forward := NewVector3d(
		eye.Normal[0]-target.Normal[0],
		eye.Normal[1]-target.Normal[1],
		eye.Normal[2]-target.Normal[2],
	)
	forward.Normalize()

	// 2. Calculate the right vector by taking the cross product of the 'up'
	//    vector and the forward vector, then normalize it.
	rightNormal := Cross(up.Normal[:], forward.Normal[:])
	right := NewVector3dFromArray(rightNormal)
	right.Normalize()

	// 3. Recalculate the true 'up' vector by taking the cross product of the
	//    forward and right vectors. This ensures all three axes are orthogonal.
	newUpNormal := Cross(forward.Normal[:], right.Normal[:])
	newUp := NewVector3dFromArray(newUpNormal)
	// No need to normalize newUp as forward and right are already unit vectors.

	// 4. Construct the rotation matrix from the calculated vectors.
	m := IdentMatrix()
	m.ThisMatrix[0][0] = right.Normal[0]
	m.ThisMatrix[0][1] = right.Normal[1]
	m.ThisMatrix[0][2] = right.Normal[2]

	m.ThisMatrix[1][0] = newUp.Normal[0]
	m.ThisMatrix[1][1] = newUp.Normal[1]
	m.ThisMatrix[1][2] = newUp.Normal[2]

	m.ThisMatrix[2][0] = forward.Normal[0]
	m.ThisMatrix[2][1] = forward.Normal[1]
	m.ThisMatrix[2][2] = forward.Normal[2]

	// 5. Convert the rotation matrix to a quaternion.
	return m
}

// MatrixToQuaternion converts a 3x3 or 4x4 rotation matrix into a quaternion.
func MatrixToQuaternion(m *Matrix) *Quaternion {
	q := &Quaternion{}
	trace := m.ThisMatrix[0][0] + m.ThisMatrix[1][1] + m.ThisMatrix[2][2]

	if trace > 0 {
		s := 0.5 / math.Sqrt(trace+1.0)
		q.W = 0.25 / s
		q.X = (m.ThisMatrix[2][1] - m.ThisMatrix[1][2]) * s
		q.Y = (m.ThisMatrix[0][2] - m.ThisMatrix[2][0]) * s
		q.Z = (m.ThisMatrix[1][0] - m.ThisMatrix[0][1]) * s
	} else {
		if m.ThisMatrix[0][0] > m.ThisMatrix[1][1] && m.ThisMatrix[0][0] > m.ThisMatrix[2][2] {
			s := 2.0 * math.Sqrt(1.0+m.ThisMatrix[0][0]-m.ThisMatrix[1][1]-m.ThisMatrix[2][2])
			q.W = (m.ThisMatrix[2][1] - m.ThisMatrix[1][2]) / s
			q.X = 0.25 * s
			q.Y = (m.ThisMatrix[0][1] + m.ThisMatrix[1][0]) / s
			q.Z = (m.ThisMatrix[0][2] + m.ThisMatrix[2][0]) / s
		} else if m.ThisMatrix[1][1] > m.ThisMatrix[2][2] {
			s := 2.0 * math.Sqrt(1.0+m.ThisMatrix[1][1]-m.ThisMatrix[0][0]-m.ThisMatrix[2][2])
			q.W = (m.ThisMatrix[0][2] - m.ThisMatrix[2][0]) / s
			q.X = (m.ThisMatrix[0][1] + m.ThisMatrix[1][0]) / s
			q.Y = 0.25 * s
			q.Z = (m.ThisMatrix[1][2] + m.ThisMatrix[2][1]) / s
		} else {
			s := 2.0 * math.Sqrt(1.0+m.ThisMatrix[2][2]-m.ThisMatrix[0][0]-m.ThisMatrix[1][1])
			q.W = (m.ThisMatrix[1][0] - m.ThisMatrix[0][1]) / s
			q.X = (m.ThisMatrix[0][2] + m.ThisMatrix[2][0]) / s
			q.Y = (m.ThisMatrix[1][2] + m.ThisMatrix[2][1]) / s
			q.Z = 0.25 * s
		}
	}
	return q
}

// Cross calculates the cross product of two 3-element vectors.
func Cross(a, b []float64) []float64 {
	return []float64{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
		0, // W component is 0 for vectors
	}
}

// NewCameraWithQuaternion creates a camera at a specific position with an
// orientation defined by a quaternion.
func NewCameraWithQuaternion(pos *Point3d, orientation *Quaternion) *Camera {
	// The camera's view matrix is the inverse of its transformation matrix.
	// For a rotation, the inverse is the conjugate of the quaternion, which
	// for a unit quaternion is (-x, -y, -z, w).
	invOrientation := &Quaternion{
		W: orientation.W,
		X: -orientation.X,
		Y: -orientation.Y,
		Z: -orientation.Z,
	}

	// The final view matrix is a rotation followed by a translation.
	rotMatrix := invOrientation.ToRotationMatrix()
	transMatrix := TransMatrix(-pos.GetX(), -pos.GetY(), -pos.GetZ())
	viewMatrix := rotMatrix.MultiplyBy(transMatrix)

	return &Camera{
		camMatrixRev:   viewMatrix,
		cameraPosition: pos,
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
		dest.ThisMatrix[x][0] = m.ThisMatrix[0][0]*sx + m.ThisMatrix[1][0]*sy + m.ThisMatrix[2][0]*sz + m.ThisMatrix[3][0]
		dest.ThisMatrix[x][1] = m.ThisMatrix[0][1]*sx + m.ThisMatrix[1][1]*sy + m.ThisMatrix[2][1]*sz + m.ThisMatrix[3][1]
		dest.ThisMatrix[x][2] = m.ThisMatrix[0][2]*sx + m.ThisMatrix[1][2]*sy + m.ThisMatrix[2][2]*sz + m.ThisMatrix[3][2]
	}
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

type Mesh struct {
	Points *Matrix
}

func NewMesh() *Mesh {
	return &Mesh{Points: NewMatrix()}
}

func (m *Mesh) AddPoint(point []float64) ([]float64, int) {
	for i, p := range m.Points.ThisMatrix {
		if p[0] == point[0] && p[1] == point[1] && p[2] == point[2] {
			return p, i
		}
	}

	pointCopy := make([]float64, len(point))
	copy(pointCopy, point)
	m.Points.AddRow(pointCopy)
	newIndex := len(m.Points.ThisMatrix) - 1
	return pointCopy, newIndex
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

func (fm *FaceMesh) AddFace(f *Face) (*Face, []int) { // Return indices
	newPoints := make([][]float64, len(f.Points))
	indices := make([]int, len(f.Points))
	for i, p := range f.Points {
		newPoints[i], indices[i] = fm.AddPoint(p) // Capture the returned index
	}
	return NewFace(newPoints, f.Col, f.GetNormal()), indices
}

func (fm *FaceMesh) Copy() *FaceMesh {
	return &FaceMesh{Mesh: *fm.Mesh.Copy()}
}

type NormalMesh struct {
	Mesh
}

func NewNormalMesh() *NormalMesh {
	return &NormalMesh{Mesh: *NewMesh()}
}

func (nm *NormalMesh) AddNormal(pnts []float64) ([]float64, int) {
	// The previous fix was a patch. This is the correct implementation.
	return nm.AddPoint(pnts)
}

func (nm *NormalMesh) Copy() *NormalMesh {
	return &NormalMesh{Mesh: *nm.Mesh.Copy()}
}

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

type BspNode struct {
	normal           []float64
	Left             *BspNode
	Right            *BspNode
	colRed           uint8
	colGreen         uint8
	colBlue          uint8
	facePointIndices []int
	xp               []float32
	yp               []float32

	normalIndex int
}

func NewBspNode(facePoints [][]float64, faceNormal []float64, faceColor color.RGBA, pointIndices []int, normalIdx int) *BspNode {
	b := &BspNode{
		normal:           faceNormal, // We can still keep the original for other uses if needed
		colRed:           faceColor.R,
		colGreen:         faceColor.G,
		colBlue:          faceColor.B,
		facePointIndices: pointIndices,
		normalIndex:      normalIdx,
		xp:               make([]float32, len(facePoints)),
		yp:               make([]float32, len(facePoints)),
	}
	return b
}

func (b *BspNode) PaintWithoutColorChange(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix) {
	b.PaintWithColor(screen, x, y, transPoints, transNormals, false)
}

// In BspNode.Paint, change the signature and logic
// func (b *BspNode) PaintWithColor(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, changeColor bool) {
// 	if len(b.facePointIndices) == 0 {
// 		return
// 	}

// 	transformedNormal := transNormals.ThisMatrix[b.normalIndex]

// 	// Look up the first transformed point of the face.
// 	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

// 	// Perform the dot product using two vectors that are now in the SAME coordinate space.
// 	// this is also used in the lighting calculation to determine the angle of the face relative to the POSITION of camera (doesn't take into account the direction of the camera).
// 	where := transformedNormal[0]*firstTransformedPoint[0] +
// 		transformedNormal[1]*firstTransformedPoint[1] +
// 		transformedNormal[2]*firstTransformedPoint[2]

// 	if where <= 0 {
// 		if b.Left != nil {
// 			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 		if b.Right != nil {
// 			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 	} else {
// 		if b.Right != nil {
// 			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}

// 		// Project the points for this polygon
// 		for i, pointIndex := range b.facePointIndices {
// 			// Look up the transformed point using its index!
// 			pnt := transPoints.ThisMatrix[pointIndex]
// 			if pnt[2] < 0 { // Z-clipping
// 				if b.Left != nil {
// 					b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 				}
// 				return
// 			}
// 			b.xp[i] = float32((400*pnt[0])/pnt[2]) + float32(x)
// 			b.yp[i] = float32((400*pnt[1])/pnt[2]) + float32(y)
// 		}

// 		polyColor := color.RGBA{R: uint8(b.colRed), G: uint8(b.colGreen), B: uint8(b.colBlue), A: 255}

// 		if changeColor {
// 			// cosTheta := where / GetLength(firstTransformedPoint)

// 			// c := 240 - int(cosTheta*240)

// 			// min := 7
// 			// r1 := clamp(int(b.colRed)-c, min, 255)
// 			// g1 := clamp(int(b.colGreen)-c, min, 255)
// 			// b1 := clamp(int(b.colBlue)-c, min, 255)
// 			// polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}

// 			cosTheta := transformedNormal[2]

// 			// The original lighting formula works perfectly with this new cosTheta.
// 			c := 240 - int(cosTheta*240)

// 			min := 7
// 			r1 := clamp(int(b.colRed)-c, min, 255)
// 			g1 := clamp(int(b.colGreen)-c, min, 255)
// 			b1 := clamp(int(b.colBlue)-c, min, 255)
// 			polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}

// 		}
// 		fillConvexPolygon(screen, b.xp, b.yp, polyColor)

// 		if b.Left != nil {
// 			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 	}
// }

// In BspNode.PaintWithColor, update the lighting model to include a spotlight effect.
// func (b *BspNode) PaintWithColor(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, changeColor bool) {
// 	if len(b.facePointIndices) == 0 {
// 		return
// 	}

// 	transformedNormal := transNormals.ThisMatrix[b.normalIndex]
// 	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

// 	// Original back-face culling logic remains unchanged.
// 	where := transformedNormal[0]*firstTransformedPoint[0] +
// 		transformedNormal[1]*firstTransformedPoint[1] +
// 		transformedNormal[2]*firstTransformedPoint[2]

// 	if where <= 0 {
// 		if b.Left != nil {
// 			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 		if b.Right != nil {
// 			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 	} else {
// 		if b.Right != nil {
// 			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}

// 		// Project the points for this polygon
// 		for i, pointIndex := range b.facePointIndices {
// 			pnt := transPoints.ThisMatrix[pointIndex]
// 			if pnt[2] < 0 { // Z-clipping
// 				if b.Left != nil {
// 					b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 				}
// 				return
// 			}
// 			b.xp[i] = float32((400*pnt[0])/pnt[2]) + float32(x)
// 			b.yp[i] = float32((400*pnt[1])/pnt[2]) + float32(y)
// 		}

// 		polyColor := color.RGBA{R: uint8(b.colRed), G: uint8(b.colGreen), B: uint8(b.colBlue), A: 255}

// 		if changeColor {
// 			// ==================== CHANGE IS HERE ====================
// 			// This lighting model combines two factors for a spotlight effect.
// 			// Higher values create a sharper, more focused spotlight cone.
// 			const spotlightConePower = 2.0

// 			// --- 1. Diffuse Factor (based on face orientation) ---
// 			// This is the same as the previous change. It's 1.0 if the face is
// 			// pointing at the camera, and 0.0 if it's at a right angle.
// 			diffuseFactor := transformedNormal[2]
// 			if diffuseFactor < 0 {
// 				diffuseFactor = 0
// 			}

// 			// --- 2. Spotlight Factor (based on face position in view) ---
// 			// This factor makes faces in the center of the view brighter than
// 			// faces at the edges. It's the cosine of the angle between the
// 			// camera's view direction (0,0,1) and the vector to the face.
// 			var spotlightFactor float64
// 			lenVecToPoint := GetLength(firstTransformedPoint)

// 			if lenVecToPoint > 0 {
// 				// cos(angle) = dot((0,0,1), vecToPoint) / |vecToPoint|
// 				cosAngle := firstTransformedPoint[2] / lenVecToPoint
// 				if cosAngle < 0 {
// 					cosAngle = 0 // Don't light things behind the camera's origin
// 				}
// 				// We raise the cosine to a power to create the cone effect.
// 				spotlightFactor = math.Pow(cosAngle, spotlightConePower)
// 			} else {
// 				// The point is exactly at the camera's origin, give it full light.
// 				spotlightFactor = 1.0
// 			}

// 			// --- 3. Combine Factors ---
// 			// The final brightness is the product of both factors. A face must be
// 			// both facing the camera AND be in the center of the spotlight to be
// 			// at its brightest.
// 			finalBrightness := diffuseFactor * spotlightFactor

// 			// Use the final brightness with the original color modification formula.
// 			// Brightness of 1.0 -> subtract 0. Brightness of 0.0 -> subtract 240.
// 			c := 240 - int(finalBrightness*240)

// 			min := 7
// 			r1 := clamp(int(b.colRed)-c, min, 255)
// 			g1 := clamp(int(b.colGreen)-c, min, 255)
// 			b1 := clamp(int(b.colBlue)-c, min, 255)
// 			polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}
// 			// ================= END OF CHANGE =====================
// 		}
// 		fillConvexPolygon(screen, b.xp, b.yp, polyColor)

// 		if b.Left != nil {
// 			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
// 		}
// 	}
// }

// In BspNode.PaintWithColor, update the lighting model to include ambient light.
func (b *BspNode) PaintWithColor(screen *ebiten.Image, x, y int, transPoints *Matrix, transNormals *Matrix, changeColor bool) {
	if len(b.facePointIndices) == 0 {
		return
	}

	transformedNormal := transNormals.ThisMatrix[b.normalIndex]
	firstTransformedPoint := transPoints.ThisMatrix[b.facePointIndices[0]]

	// Original back-face culling logic remains unchanged.
	where := transformedNormal[0]*firstTransformedPoint[0] +
		transformedNormal[1]*firstTransformedPoint[1] +
		transformedNormal[2]*firstTransformedPoint[2]

	if where <= 0 {
		if b.Left != nil {
			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
		if b.Right != nil {
			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
	} else {
		if b.Right != nil {
			b.Right.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}

		// Project the points for this polygon
		for i, pointIndex := range b.facePointIndices {
			pnt := transPoints.ThisMatrix[pointIndex]
			if pnt[2] < 0 { // Z-clipping
				if b.Left != nil {
					b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
				}
				return
			}
			b.xp[i] = float32((400*pnt[0])/pnt[2]) + float32(x)
			b.yp[i] = float32((400*pnt[1])/pnt[2]) + float32(y)
		}

		polyColor := color.RGBA{R: uint8(b.colRed), G: uint8(b.colGreen), B: uint8(b.colBlue), A: 255}

		if changeColor {
			// Define Lighting Parameters
			// The minimum brightness for any surface (e.g., 20% light).
			const ambientLight = 0.65
			// Higher values create a sharper, more focused spotlight cone.
			const spotlightConePower = 10.0
			// The amount of light available for the spotlight effect is the total (1.0)
			// minus the ambient light we've already added.
			const spotlightLightAmount = 1.0 - ambientLight

			// Diffuse Factor (based on face orientation)
			diffuseFactor := transformedNormal[2]
			if diffuseFactor < 0 {
				diffuseFactor = 0
			}

			// potlight Factor (based on face position in view)
			var spotlightFactor float64
			lenVecToPoint := GetLength(firstTransformedPoint)

			if lenVecToPoint > 0 {
				cosAngle := firstTransformedPoint[2] / lenVecToPoint
				if cosAngle < 0 {
					cosAngle = 0
				}
				spotlightFactor = math.Pow(cosAngle, spotlightConePower)
			} else {
				spotlightFactor = 1.0
			}

			// The spotlight's contribution to brightness is scaled by the amount
			// of light available for it.
			spotlightBrightness := diffuseFactor * spotlightFactor * spotlightLightAmount

			// Add the ambient light to the calculated spotlight brightness.
			// This ensures a minimum light level for the entire scene.
			finalBrightness := ambientLight + spotlightBrightness

			// Use the final brightness with the original color modification formula.
			// A brightness of 1.0 means no color change (subtract 0).
			// A brightness of 0.0 means maximum darkness (subtract 240).
			c := 240 - int(finalBrightness*240)

			min := 7
			r1 := clamp(int(b.colRed)-c, min, 255)
			g1 := clamp(int(b.colGreen)-c, min, 255)
			b1 := clamp(int(b.colBlue)-c, min, 255)
			polyColor = color.RGBA{R: uint8(r1), G: uint8(g1), B: uint8(b1), A: 255}

		}
		fillConvexPolygon(screen, b.xp, b.yp, polyColor)

		if b.Left != nil {
			b.Left.PaintWithColor(screen, x, y, transPoints, transNormals, changeColor)
		}
	}
}

func clamp(value, min, max int) int {
	if value < min {
		return min
	}
	if value > max {
		return max
	}
	return value
}

func (o *Object3d) PaintSolid(screen *ebiten.Image, x, y int, lightingChange bool) {
	if o.root != nil {
		o.root.PaintWithColor(screen, x, y, o.transFaceMesh.Points, o.transNormalMesh.Points, lightingChange)
	}
}

type Object3d struct {
	faceMesh        *FaceMesh
	normalMesh      *NormalMesh
	transFaceMesh   *FaceMesh
	transNormalMesh *NormalMesh
	theFaces        *FaceStore
	root            *BspNode
	rotMatrix       *Matrix
}

func NewObject_3d() *Object3d {
	return &Object3d{
		transFaceMesh:   NewFaceMesh(),
		transNormalMesh: NewNormalMesh(),
		theFaces:        NewFaceStore(),
		rotMatrix:       IdentMatrix(),
	}
}

func (o *Object3d) Clone() *Object3d {
	clone := &Object3d{
		// shared
		faceMesh:   o.faceMesh,
		normalMesh: o.normalMesh,
		theFaces:   o.theFaces,
		root:       o.root,

		// instance-specific
		transFaceMesh:   o.transFaceMesh.Copy(),
		transNormalMesh: o.transNormalMesh.Copy(),
		rotMatrix:       IdentMatrix(),
	}
	return clone
}

func (o *Object3d) Finished() {
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

func (o *Object3d) ApplyMatrixBatch(m *Matrix) {
	o.rotMatrix = m.MultiplyBy(o.rotMatrix)
}

func (o *Object3d) ApplyMatrixTemp(aMatrix *Matrix) {
	rotMatrixTemp := aMatrix.MultiplyBy(o.rotMatrix)

	// Use the new, correct method to transform the normals (rotation only).
	rotMatrixTemp.TransformNormals(o.normalMesh.Points, o.transNormalMesh.Points)

	// The renormalization step you already had is still good practice to prevent
	// floating-point drift from affecting the length of the normal.
	for _, n := range o.transNormalMesh.Points.ThisMatrix {
		v := NewVector3dFromArray(n)
		v.Normalize()
		copy(n, v.Normal[:])
	}

	// Use the original method to transform the vertex positions (rotation and translation).
	rotMatrixTemp.TransformObj(o.faceMesh.Points, o.transFaceMesh.Points)
}

func (o *Object3d) createBspTree(faces *FaceStore, newFaces *FaceMesh, newNormMesh *NormalMesh) *BspNode {
	if faces.FaceCount() == 0 {
		return nil
	}

	parentFace := o.choosePlane(faces)
	originalNormal, normalIndex := newNormMesh.AddNormal(parentFace.GetNormal())
	parentFace.SetNormal(originalNormal)
	newFace, parentIndices := newFaces.AddFace(parentFace)
	parent := NewBspNode(newFace.Points, newFace.GetNormal(), newFace.Col, parentIndices, normalIndex)
	pPlane := NewPlane(newFace, newFace.GetNormal())

	// Create two new lists to hold the faces that fall on either side of the plane.
	fvLeft := NewFaceStore()
	fvRight := NewFaceStore()

	// Partition the *remaining* faces against the splitting plane.
	for a := 0; a < faces.FaceCount(); a++ {
		currentFace := faces.GetFace(a)
		if pPlane.FaceIntersect(currentFace) {
			// If the face is split by the plane...
			split := pPlane.SplitFace(currentFace)
			if split == nil {
				continue
			}
			for _, facePart := range split {
				if facePart != nil && len(facePart.Points) > 0 {
					if pPlane.Where(facePart) <= 0 {
						fvLeft.AddFace(NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal()))
					} else {
						fvRight.AddFace(NewFace(facePart.Points, currentFace.Col, currentFace.GetNormal()))
					}
				}
			}
		} else {
			// If the face is entirely on one side...
			w := pPlane.Where(currentFace)
			if w <= 0 {
				fvLeft.AddFace(currentFace)
			} else {
				fvRight.AddFace(currentFace)
			}
		}
	}

	// Build the left and right sub-trees from the new lists.
	if fvLeft.FaceCount() > 0 {
		parent.Left = o.createBspTree(fvLeft, newFaces, newNormMesh)
	}
	if fvRight.FaceCount() > 0 {
		parent.Right = o.createBspTree(fvRight, newFaces, newNormMesh)
	}

	return parent
}

func (o *Object3d) choosePlane(fs *FaceStore) *Face {
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

func NewCube() *Object3d {
	obj := NewObject_3d()
	s := 40.0 // size

	// The 8 vertices of the cube
	points := [][3]float64{
		{-s, -s, -s}, // 0
		{s, -s, -s},  // 1
		{s, s, -s},   // 2
		{-s, s, -s},  // 3
		{-s, -s, s},  // 4
		{s, -s, s},   // 5
		{s, s, s},    // 6
		{-s, s, s},   // 7
	}

	colors := []color.RGBA{
		{255, 0, 0, 255},   // Red
		{0, 255, 0, 255},   // Green
		{0, 0, 255, 255},   // Blue
		{255, 255, 0, 255}, // Yellow
		{0, 255, 255, 255}, // Cyan
		{255, 0, 255, 255}, // Magenta
	}

	// Define quads with a consistent Counter-Clockwise (CCW) winding order
	// This ensures the calculated normals point outwards from the cube.
	quads := [][]int{
		{0, 3, 2, 1}, // Front face  (Normal: 0, 0, -1)
		{1, 2, 6, 5}, // Right face  (Normal: 1, 0, 0)
		{5, 6, 7, 4}, // Back face   (Normal: 0, 0, 1)
		{4, 7, 3, 0}, // Left face   (Normal: -1, 0, 0)
		{3, 7, 6, 2}, // Top face    (Normal: 0, 1, 0)
		{4, 0, 1, 5}, // Bottom face (Normal: 0, -1, 0)
	}

	for i, q := range quads {
		face := NewFace(nil, colors[i], nil)
		// Vertices are added in the specified order to ensure correct normal
		face.AddPoint(points[q[0]][0], points[q[0]][1], points[q[0]][2])
		face.AddPoint(points[q[1]][0], points[q[1]][1], points[q[1]][2])
		face.AddPoint(points[q[2]][0], points[q[2]][1], points[q[2]][2])
		face.AddPoint(points[q[3]][0], points[q[3]][1], points[q[3]][2])

		// Use FACE_NORMAL because our winding order is now correct.
		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	obj.Finished()
	return obj
}

// NewRectangle creates a new Object_3d in the shape of a cuboid (a 3D rectangle).
// It's centered at the origin and has the specified dimensions and color.
func NewRectangle(width, height, length float64, clr color.RGBA) *Object3d {
	// create a new, empty object to populate.
	obj := NewObject_3d()

	// calculate half-dimensions for centering the rectangle at the origin.
	w2 := width / 2.0
	h2 := height / 2.0
	l2 := length / 2.0

	// Define the 8 vertices of the cuboid based on the dimensions.
	points := [][3]float64{
		{-w2, -h2, -l2}, // 0: front-bottom-left
		{w2, -h2, -l2},  // 1: front-bottom-right
		{w2, h2, -l2},   // 2: front-top-right
		{-w2, h2, -l2},  // 3: front-top-left
		{-w2, -h2, l2},  // 4: back-bottom-left
		{w2, -h2, l2},   // 5: back-bottom-right
		{w2, h2, l2},    // 6: back-top-right
		{-w2, h2, l2},   // 7: back-top-left
	}

	quads := [][]int{
		{0, 3, 2, 1}, // Front face  (Normal: 0, 0, -1)
		{1, 2, 6, 5}, // Right face  (Normal: 1, 0, 0)
		{5, 6, 7, 4}, // Back face   (Normal: 0, 0, 1)
		{4, 7, 3, 0}, // Left face   (Normal: -1, 0, 0)
		{3, 7, 6, 2}, // Top face    (Normal: 0, 1, 0)
		{4, 0, 1, 5}, // Bottom face (Normal: 0, -1, 0)
	}

	// Create each face and add it to the object.
	for _, q := range quads {
		// Use the single color passed into the function for every face.
		face := NewFace(nil, clr, nil)

		// Add the four vertices that make up this face.
		face.AddPoint(points[q[0]][0], points[q[0]][1], points[q[0]][2])
		face.AddPoint(points[q[1]][0], points[q[1]][1], points[q[1]][2])
		face.AddPoint(points[q[2]][0], points[q[2]][1], points[q[2]][2])
		face.AddPoint(points[q[3]][0], points[q[3]][1], points[q[3]][2])

		// Finalize the face's geometry. Use FACE_NORMAL because the winding
		// order is defined correctly to produce an outward normal.
		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	obj.Finished()

	return obj
}

// NewSubdividedRectangle creates a new Object_3d in the shape of a cuboid where each face
// is subdivided into a grid of smaller triangles.
// width, height, length: The dimensions of the cuboid.
// clr: The color for all faces of the cuboid.
// subdivisions: The number of divisions along each edge of a face. For example, a value
//
//	of 2 will split a face into a 2x2 grid of quads (8 triangles). A value of 1
//	will result in one quad per face (2 triangles).
func NewSubdividedRectangle(width, height, length float64, clr color.RGBA, subdivisions int) *Object3d {
	obj := NewObject_3d()
	w2, h2, l2 := width/2.0, height/2.0, length/2.0

	if subdivisions < 1 {
		subdivisions = 1 // Ensure at least one subdivision.
	}

	// generateFace is a helper function that constructs one of the six faces of the cuboid.
	// It takes an origin point and two vectors (u, v) that define the plane and dimensions
	// of the face. It then creates a grid of vertices and generates triangles.
	generateFace := func(origin, u, v [3]float64) {
		// Create a grid of vertices for the current face.
		vertices := make([][][3]float64, subdivisions+1)
		for i := range vertices {
			vertices[i] = make([][3]float64, subdivisions+1)
			for j := range vertices[i] {
				// Calculate the position of vertex (i, j) on the grid plane.
				ui := float64(i) / float64(subdivisions)
				vj := float64(j) / float64(subdivisions)
				vertices[i][j] = [3]float64{
					origin[0] + ui*u[0] + vj*v[0],
					origin[1] + ui*u[1] + vj*v[1],
					origin[2] + ui*u[2] + vj*v[2],
				}
			}
		}

		// Create two triangles for each quad in the subdivision grid.
		for i := 0; i < subdivisions; i++ {
			for j := 0; j < subdivisions; j++ {
				// Get the four corner vertices of the current quad.
				p1 := vertices[i][j]
				p2 := vertices[i+1][j]
				p3 := vertices[i+1][j+1]
				p4 := vertices[i][j+1]

				// Create the first triangle for the quad (p1, p2, p3).
				face1 := NewFace(nil, clr, nil)
				face1.AddPoint(p1[0], p1[1], p1[2])
				face1.AddPoint(p2[0], p2[1], p2[2])
				face1.AddPoint(p3[0], p3[1], p3[2])
				// The vertices are wound counter-clockwise to produce an outward-facing normal.
				face1.Finished(FACE_NORMAL)
				obj.theFaces.AddFace(face1)

				// Create the second triangle for the quad (p1, p3, p4).
				face2 := NewFace(nil, clr, nil)
				face2.AddPoint(p1[0], p1[1], p1[2])
				face2.AddPoint(p3[0], p3[1], p3[2])
				face2.AddPoint(p4[0], p4[1], p4[2])
				face2.Finished(FACE_NORMAL)
				obj.theFaces.AddFace(face2)
			}
		}
	}

	// Define the 6 faces of the cuboid by specifying their origin, u-vector, and v-vector.
	// The u and v vectors are chosen so their cross-product (which determines the normal)
	// points outwards from the center of the cuboid.

	// Back face (-Z direction)
	generateFace([3]float64{w2, -h2, -l2}, [3]float64{-width, 0, 0}, [3]float64{0, height, 0})

	// Front face (+Z direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{width, 0, 0}, [3]float64{0, height, 0})

	// Left face (-X direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{0, 0, -length}, [3]float64{0, height, 0})

	// Right face (+X direction)
	generateFace([3]float64{w2, -h2, -l2}, [3]float64{0, 0, length}, [3]float64{0, height, 0})

	// Bottom face (-Y direction)
	generateFace([3]float64{-w2, -h2, l2}, [3]float64{width, 0, 0}, [3]float64{0, 0, -length})

	// Top face (+Y direction)
	generateFace([3]float64{-w2, h2, -l2}, [3]float64{width, 0, 0}, [3]float64{0, 0, length})

	// Finalize the object by building its BSP tree.
	obj.Finished()
	return obj
}

// NewSphere creates a new Object_3d in the shape of an icosphere.
// An icosphere is a sphere made of a mesh of triangles, which is more
// uniform than a traditional UV sphere.
func NewSphere(radius float64, subdivisions int, clr color.RGBA) *Object3d {
	obj := NewObject_3d()

	// Define the 12 vertices of an Icosahedron.
	// An icosahedron is a 20-sided polyhedron that forms the base of our sphere.
	// The 't' value is the golden ratio, which helps define the vertex positions.
	t := (1.0 + math.Sqrt(5.0)) / 2.0

	vertices := [][3]float64{
		{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0},
		{0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
		{t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1},
	}

	// Define the 20 triangular faces of the Icosahedron using indices
	// into the vertex list. The order is important for correct normals.
	faces := [][]int{
		{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
		{1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
		{3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
		{4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1},
	}

	// Create initial list of faces for subdivision
	subdivisionFaces := make([][3][3]float64, len(faces))
	for i, faceIndices := range faces {
		subdivisionFaces[i] = [3][3]float64{
			vertices[faceIndices[0]],
			vertices[faceIndices[1]],
			vertices[faceIndices[2]],
		}
	}

	// Subdivide the faces recursively to make the sphere smoother.
	for i := 0; i < subdivisions; i++ {
		newFaces := make([][3][3]float64, 0)
		for _, face := range subdivisionFaces {
			v1 := face[0]
			v2 := face[1]
			v3 := face[2]

			// Calculate the midpoint of each edge of the triangle.
			a := [3]float64{(v1[0] + v2[0]) / 2, (v1[1] + v2[1]) / 2, (v1[2] + v2[2]) / 2}
			b := [3]float64{(v2[0] + v3[0]) / 2, (v2[1] + v3[1]) / 2, (v2[2] + v3[2]) / 2}
			c := [3]float64{(v3[0] + v1[0]) / 2, (v3[1] + v1[1]) / 2, (v3[2] + v1[2]) / 2}

			// Split the original triangle into 4 new ones.
			newFaces = append(newFaces, [3][3]float64{v1, a, c})
			newFaces = append(newFaces, [3][3]float64{v2, b, a})
			newFaces = append(newFaces, [3][3]float64{v3, c, b})
			newFaces = append(newFaces, [3][3]float64{a, b, c})
		}
		subdivisionFaces = newFaces
	}

	// Create the final faces for the Object_3d
	for _, faceVerts := range subdivisionFaces {
		face := NewFace(nil, clr, nil)

		for _, v := range faceVerts {
			// Normalize the vertex to push it onto the surface of a unit sphere.
			length := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
			normX := v[0] / length
			normY := v[1] / length
			normZ := v[2] / length

			// Scale by the desired radius and add the point.
			face.AddPoint(normX*radius, normY*radius, normZ*radius)
		}

		face.Finished(FACE_REVERSE)
		obj.theFaces.AddFace(face)
	}

	obj.Finished()
	return obj
}

type Camera struct {
	camMatrixRev   *Matrix
	cameraPosition *Point3d
}

func NewCamera(xp, yp, zp, xa, ya, za float64) *Camera {
	c := &Camera{}
	x := NewRotationMatrix(ROTX, -xa)
	y := NewRotationMatrix(ROTY, -ya)
	z := NewRotationMatrix(ROTZ, -za)
	c.camMatrixRev = z.MultiplyBy(y)
	c.camMatrixRev = c.camMatrixRev.MultiplyBy(x)
	c.cameraPosition = NewPoint3d(xp, yp, zp)

	return c
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

func (c *Camera) GetMatrix() *Matrix {
	return c.camMatrixRev
}

func (c *Camera) SetMatrix(m *Matrix) {
	c.camMatrixRev = m
}

func (c *Camera) AddAngle(x, y float64) {
	rotY := NewRotationMatrix(ROTY, -y)
	rotX := NewRotationMatrix(ROTX, -x)
	c.camMatrixRev = rotY.MultiplyBy(rotX).MultiplyBy(c.camMatrixRev)
}

// create a camera at position lookin at another point
func NewCameraLookAt(x, y, z, lookX, lookY, lookZ float64) *Camera {
	// Calculate the direction vector from the camera to the look-at point
	dx := (lookX - x)
	dy := (lookY - y)
	dz := (lookZ - z)

	// Normalize the direction vector
	length := math.Sqrt(dx*dx + dy*dy + dz*dz)
	if length == 0 {
		return NewCamera(x, y, z, 0, 0, 0)
	}
	dx /= length
	dy /= length
	dz /= length

	// Calculate angles based on the direction vector
	xa := math.Atan2(dy, math.Sqrt(dx*dx+dz*dz)) * (180 / math.Pi)
	ya := math.Atan2(dx, dz) * (180 / math.Pi)

	return NewCamera(x, y, z, xa, ya, 0)
}

type World_3d struct {
	objects               []*Object3d
	objXpos               []float64
	objYpos               []float64
	objZpos               []float64
	cameras               []*Camera
	camXpos               []float64
	camYpos               []float64
	camZpos               []float64
	currentCamera         int
	objectToDrawFirst     []*Object3d
	objectToDrawFirstXpos []float64
	objectToDrawFirstYpos []float64
	objectToDrawFirstZpos []float64
}

func NewWorld_3d() *World_3d {
	return &World_3d{
		currentCamera: -1,
	}
}

func (w *World_3d) AddObject(obj *Object3d, x, y, z float64) {
	w.objects = append(w.objects, obj)
	w.objXpos = append(w.objXpos, x)
	w.objYpos = append(w.objYpos, y)
	w.objZpos = append(w.objZpos, z)
}

// AddObjectDrawFirst
func (w *World_3d) AddObjectDrawFirst(obj *Object3d, x, y, z float64) {
	w.objectToDrawFirst = append(w.objectToDrawFirst, obj)
	w.objectToDrawFirstXpos = append(w.objectToDrawFirstXpos, x)
	w.objectToDrawFirstYpos = append(w.objectToDrawFirstYpos, y)
	w.objectToDrawFirstZpos = append(w.objectToDrawFirstZpos, z)
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

	// sort objects by distance to camera
	var sortedIndices []int
	for i := range w.objects {
		sortedIndices = append(sortedIndices, i)
	}
	sort.Slice(sortedIndices, func(i, j int) bool {
		distanceI := math.Sqrt(math.Pow(w.objXpos[sortedIndices[i]]-camX, 2) +
			math.Pow(w.objYpos[sortedIndices[i]]-camY, 2) +
			math.Pow(w.objZpos[sortedIndices[i]]-camZ, 2))

		distanceJ := math.Sqrt(math.Pow(w.objXpos[sortedIndices[j]]-camX, 2) +
			math.Pow(w.objYpos[sortedIndices[j]]-camY, 2) +
			math.Pow(w.objZpos[sortedIndices[j]]-camZ, 2))

		return distanceI > distanceJ
	})

	// Draw objects that should be drawn first
	for i, obj := range w.objectToDrawFirst {
		m := TransMatrix(w.objectToDrawFirstXpos[i]-camX, w.objectToDrawFirstYpos[i]-camY, w.objectToDrawFirstZpos[i]-camZ)
		obj.ApplyMatrixTemp(cam.camMatrixRev.MultiplyBy(m))
		obj.PaintSolid(screen, xsize/2, ysize/2, true)
	}

	for _, i := range sortedIndices {
		obj := w.objects[i]
		m := TransMatrix(w.objXpos[i]-camX, w.objYpos[i]-camY, w.objZpos[i]-camZ)
		obj.ApplyMatrixTemp(cam.camMatrixRev.MultiplyBy(m))
		obj.PaintSolid(screen, xsize/2, ysize/2, true)
	}
}

var (
	whiteImage = ebiten.NewImage(3, 3)
)

func init() {
	whiteImage.Fill(color.White)
}

type Game struct {
	world        *World_3d
	i, p         float64
	cube         *Object3d
	lastX, lastY int
	dragged      bool
}

func NewGame() *Game {
	g := &Game{}
	log.Println("Initializing World...")
	g.world = NewWorld_3d()
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

// NewUVSphere creates a sphere based on latitude/longitude rings (sectors and stacks).
// This structure allows for a perfectly straight horizontal stripe.
func NewUVSphere(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Object3d {
	obj := NewObject_3d()

	// We loop through stacks (latitude) and sectors (longitude).
	vertices := make([][3]float64, 0)
	for i := 0; i <= stacks; i++ {
		stackAngle := math.Pi/2 - float64(i)*math.Pi/float64(stacks) // phi
		xy := radius * math.Cos(stackAngle)
		z := radius * math.Sin(stackAngle)

		for j := 0; j <= sectors; j++ {
			sectorAngle := float64(j) * 2 * math.Pi / float64(sectors) // theta
			x := xy * math.Cos(sectorAngle)
			y := xy * math.Sin(sectorAngle)
			vertices = append(vertices, [3]float64{x, y, z})
		}
	}

	// Determine the start and end stacks for the stripe.
	// The stripe is centered around the equator (the middle stack).
	middleStack := stacks / 2
	stripeStart := middleStack - (stripeStacks / 2)
	stripeEnd := middleStack + (stripeStacks / 2)

	for i := 0; i < stacks; i++ {
		// Determine the color for this entire ring of faces.
		var faceColor color.RGBA
		if i >= stripeStart && i < stripeEnd {
			faceColor = stripeClr
		} else {
			faceColor = bodyClr
		}

		k1 := i * (sectors + 1)
		k2 := k1 + sectors + 1

		for j := 0; j < sectors; j++ {
			k1j := k1 + j
			k1j1 := k1j + 1
			k2j := k2 + j
			k2j1 := k2j + 1

			// For each quad, we create two triangles.
			// Special handling for the poles.
			if i != 0 {
				// First triangle of the quad
				f1 := NewFace(nil, faceColor, nil)
				f1.AddPoint(vertices[k1j][0], vertices[k1j][1], vertices[k1j][2])
				f1.AddPoint(vertices[k2j][0], vertices[k2j][1], vertices[k2j][2])
				f1.AddPoint(vertices[k1j1][0], vertices[k1j1][1], vertices[k1j1][2])
				f1.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(f1)
			}

			if i != (stacks - 1) {
				// Second triangle of the quad
				f2 := NewFace(nil, faceColor, nil)
				f2.AddPoint(vertices[k1j1][0], vertices[k1j1][1], vertices[k1j1][2])
				f2.AddPoint(vertices[k2j][0], vertices[k2j][1], vertices[k2j][2])
				f2.AddPoint(vertices[k2j1][0], vertices[k2j1][1], vertices[k2j1][2])
				f2.Finished(FACE_REVERSE)
				obj.theFaces.AddFace(f2)
			}
		}
	}

	// 3. Finalize the object by building its BSP tree.
	obj.Finished()
	return obj
}

func NewUVSphere2(radius float64, sectors, stacks int, bodyClr, stripeClr color.RGBA, stripeStacks int) *Object3d {
	obj := NewObject_3d()

	vertices := make([][3]float64, 0)
	for i := 0; i <= stacks; i++ {
		stackAngle := math.Pi/2 - float64(i)*math.Pi/float64(stacks)
		xy := radius * math.Cos(stackAngle)
		z := radius * math.Sin(stackAngle)

		for j := 0; j <= sectors; j++ {
			sectorAngle := float64(j) * 2 * math.Pi / float64(sectors)
			x := xy * math.Cos(sectorAngle)
			y := xy * math.Sin(sectorAngle)
			vertices = append(vertices, [3]float64{x, y, z})
		}
	}

	middleStack := stacks / 2
	stripeStart := middleStack - (stripeStacks / 2)
	stripeEnd := middleStack + (stripeStacks / 2)

	// A helper function to create a face and set its normal correctly
	createFace := func(p1, p2, p3 [3]float64, clr color.RGBA) *Face {
		face := NewFace(nil, clr, nil)
		face.AddPoint(p1[0], p1[1], p1[2])
		face.AddPoint(p2[0], p2[1], p2[2])
		face.AddPoint(p3[0], p3[1], p3[2])

		// Find the center of the triangle face.
		centerX := (p1[0] + p2[0] + p3[0]) / 3.0
		centerY := (p1[1] + p2[1] + p3[1]) / 3.0
		centerZ := (p1[2] + p2[2] + p3[2]) / 3.0

		// The normal is the vector from the sphere's center (0,0,0) to the face's center, normalized.
		normalVec := []float64{centerX, centerY, centerZ, 0.0}
		length := GetLength(normalVec)
		if length > 0 {
			normalVec[0] /= length
			normalVec[1] /= length
			normalVec[2] /= length
		}

		normalVec[0] *= -1
		normalVec[1] *= -1
		normalVec[2] *= -1

		// Explicitly set this perfect normal on the face.
		face.SetNormal(normalVec)
		// We no longer need the FACE_NORMAL/FACE_REVERSE flag.
		face.Finished(FACE_REVERSE)
		return face
	}

	for i := 0; i < stacks; i++ {
		var faceColor color.RGBA
		if i >= stripeStart && i < stripeEnd {
			faceColor = stripeClr
		} else {
			faceColor = bodyClr
		}

		k1 := i * (sectors + 1)
		k2 := k1 + sectors + 1

		for j := 0; j < sectors; j++ {
			// Get the four vertices for the quad on the sphere surface
			v1 := vertices[k1+j]
			v2 := vertices[k2+j]
			v3 := vertices[k1+j+1]
			v4 := vertices[k2+j+1]

			if i != 0 {
				// First triangle of the quad: v1, v2, v3
				obj.theFaces.AddFace(createFace(v1, v2, v3, faceColor))
			}

			if i != (stacks - 1) {
				// Second triangle of the quad: v3, v2, v4
				obj.theFaces.AddFace(createFace(v3, v2, v4, faceColor))
			}
		}
	}

	obj.Finished()
	return obj
}

func (g *Game) Update() error {
	moveSpeed := 15.0

	// Automatic object rotation
	g.i += 0.02
	g.p += 0.05
	for index, obj := range g.world.objects {
		s := NewForwardMatrix(math.Sin(g.i+float64(index))/8/5, math.Cos(g.p)/8/5, math.Cos(g.i+1.14)/8/5)
		obj.ApplyMatrixBatch(s)
	}

	// Get the current camera
	cam := g.world.cameras[g.world.currentCamera]
	if cam != nil {
		// The camera's view matrix (`CamMatrixRev`) is the inverse of its world
		// transformation matrix. For a pure rotation matrix, the inverse is the
		// transpose. Therefore, the columns of the view matrix represent the
		// camera's axes (right, up, forward) in world space.

		// // Right vector is the first column of the view matrix.
		// rightVecX := cam.camMatrixRev.ThisMatrix[0][0]
		// rightVecY := cam.camMatrixRev.ThisMatrix[1][0]
		// rightVecZ := cam.camMatrixRev.ThisMatrix[2][0]

		// // Forward vector is the third column of the view matrix.
		// // Note: In a right-handed coordinate system, the camera looks down its
		// // negative Z-axis, but the "forward" direction for movement is typically
		// // along the positive Z-axis of the camera's coordinate space.
		// forwardVecX := cam.camMatrixRev.ThisMatrix[0][2]
		// forwardVecY := cam.camMatrixRev.ThisMatrix[1][2]
		// forwardVecZ := cam.camMatrixRev.ThisMatrix[2][2]

		// // Handle keyboard input for movement
		// if ebiten.IsKeyPressed(ebiten.KeyW) { // Move forward
		// 	g.world.camXpos[g.world.currentCamera] += forwardVecX * moveSpeed
		// 	g.world.camYpos[g.world.currentCamera] += forwardVecY * moveSpeed
		// 	g.world.camZpos[g.world.currentCamera] += forwardVecZ * moveSpeed
		// }
		// if ebiten.IsKeyPressed(ebiten.KeyS) { // Move backward
		// 	g.world.camXpos[g.world.currentCamera] -= forwardVecX * moveSpeed
		// 	g.world.camYpos[g.world.currentCamera] -= forwardVecY * moveSpeed
		// 	g.world.camZpos[g.world.currentCamera] -= forwardVecZ * moveSpeed
		// }
		// if ebiten.IsKeyPressed(ebiten.KeyA) { // Strafe left
		// 	g.world.camXpos[g.world.currentCamera] -= rightVecX * moveSpeed
		// 	g.world.camYpos[g.world.currentCamera] -= rightVecY * moveSpeed
		// 	g.world.camZpos[g.world.currentCamera] -= rightVecZ * moveSpeed
		// }
		// if ebiten.IsKeyPressed(ebiten.KeyD) { // Strafe right
		// 	g.world.camXpos[g.world.currentCamera] += rightVecX * moveSpeed
		// 	g.world.camYpos[g.world.currentCamera] += rightVecY * moveSpeed
		// 	g.world.camZpos[g.world.currentCamera] += rightVecZ * moveSpeed
		// }

		// The camera's view matrix
		viewMatrix := cam.GetMatrix().ThisMatrix

		// Right vector is the FIRST ROW of the matrix
		rightVecX := viewMatrix[0][0]
		rightVecY := viewMatrix[0][1]
		rightVecZ := viewMatrix[0][2]

		// Forward vector is the THIRD ROW of the matrix
		// (In a right-handed view matrix, the Z-axis row points away from what you're looking at)
		forwardVecX := viewMatrix[2][0]
		forwardVecY := viewMatrix[2][1]
		forwardVecZ := viewMatrix[2][2]

		// Handle keyboard input for movement
		if ebiten.IsKeyPressed(ebiten.KeyW) { // Move forward
			// Pressing W should move you IN THE OPPOSITE direction of the forward vector
			cam.SetCameraPosition(
				cam.GetPosition().GetX()-forwardVecX*moveSpeed,
				cam.GetPosition().GetY()-forwardVecY*moveSpeed,
				cam.GetPosition().GetZ()-forwardVecZ*moveSpeed,
			)
		}
		if ebiten.IsKeyPressed(ebiten.KeyS) { // Move backward
			cam.SetCameraPosition(
				cam.GetPosition().GetX()+forwardVecX*moveSpeed,
				cam.GetPosition().GetY()+forwardVecY*moveSpeed,
				cam.GetPosition().GetZ()+forwardVecZ*moveSpeed,
			)
		}
		if ebiten.IsKeyPressed(ebiten.KeyA) { // Strafe left
			cam.SetCameraPosition(
				cam.GetPosition().GetX()-rightVecX*moveSpeed,
				cam.GetPosition().GetY()-rightVecY*moveSpeed,
				cam.GetPosition().GetZ()-rightVecZ*moveSpeed,
			)
		}
		if ebiten.IsKeyPressed(ebiten.KeyD) { // Strafe right
			cam.SetCameraPosition(
				cam.GetPosition().GetX()+rightVecX*moveSpeed,
				cam.GetPosition().GetY()+rightVecY*moveSpeed,
				cam.GetPosition().GetZ()+rightVecZ*moveSpeed,
			)
		}
	}
	// --- END of new code for camera movement ---

	// Mouse look controls (existing code)
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

// func main() {
// 	rand.Seed(time.Now().UnixNano())
// 	game := NewGame()
// 	ebiten.SetWindowSize(screenWidth, screenHeight)
// 	ebiten.SetWindowTitle("sie3d")
// 	if err := ebiten.RunGame(game); err != nil {
// 		log.Fatal(err)
// 	}
// }

func fillConvexPolygon(screen *ebiten.Image, xp, yp []float32, clr color.RGBA) {
	if len(xp) < 3 {
		return
	}

	indices := make([]uint16, 0, (len(xp)-2)*3)
	for i := 2; i < len(xp); i++ {
		indices = append(indices, 0, uint16(i-1), uint16(i))
	}

	vertices := make([]ebiten.Vertex, len(xp))
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

func LoadObjectFromDXFFile(fileName string, reverse int) (*Object3d, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, fmt.Errorf("could not open DXF file %s: %w", fileName, err)
	}
	defer file.Close()

	obj, err := NewObjectFromDXF(file, reverse)
	if err != nil {
		return nil, fmt.Errorf("error parsing DXF file %s: %w", fileName, err)
	}

	return obj, nil
}

// NewObjectFromDXF creates a new Object_3d by reading a simplified DXF file from
// the given reader. It returns a fully constructed object with its BSP tree
// already built, or an error if the file cannot be parsed.
func NewObjectFromDXF(reader io.Reader, reverse int) (*Object3d, error) {
	// 1. Create a new, empty object to populate.
	obj := NewObject_3d()

	scanner := bufio.NewScanner(reader)

	// Helper function to read the next line and parse it as a float64
	readFloatLine := func() (float64, error) {
		if !scanner.Scan() {
			if err := scanner.Err(); err != nil {
				return 0, err
			}
			return 0, io.EOF // Clean end of file
		}
		val, err := strconv.ParseFloat(strings.TrimSpace(scanner.Text()), 64)
		if err != nil {
			return 0, fmt.Errorf("could not parse float value '%s': %w", scanner.Text(), err)
		}
		return val, nil
	}

	for scanner.Scan() {
		if !strings.HasPrefix(scanner.Text(), "3DFACE") {
			continue
		}

		for i := 0; i < 3; i++ {
			if !scanner.Scan() {
				// On error, return a nil object and the error
				return nil, fmt.Errorf("unexpected end of file while parsing 3DFACE header")
			}
		}

		faceColor := color.RGBA{
			R: uint8(rand.Intn(256)),
			G: uint8(rand.Intn(256)),
			B: uint8(rand.Intn(256)),
			A: 255,
		}
		aFace := NewFace(nil, faceColor, nil)

		for c := 0; c < 4; c++ {
			x, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading X coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			y, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Y coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			z, err := readFloatLine()
			if err != nil {
				return nil, fmt.Errorf("error reading Z coordinate for vertex %d: %w", c, err)
			}
			scanner.Scan() // Skip line

			aFace.AddPoint(x, y, z)
		}

		aFace.Finished(reverse)

		// 2. Add the face to the new object we created at the start.
		obj.theFaces.AddFace(aFace)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading from DXF source: %w", err)
	}

	// 3. Finalize the new object by building its BSP tree.
	obj.Finished()

	// 4. On success, return the fully populated object and a nil error.
	return obj, nil
}
