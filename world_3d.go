package gosie3d

import (
	"math"
	"sort"

	"github.com/hajimehoshi/ebiten/v2"
)

const UP_DIR = -1.0

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
	// draw last
	objectToDrawLast     []*Object3d
	objectToDrawLastXpos []float64
	objectToDrawLastYpos []float64
	objectToDrawLastZpos []float64
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

// AddObjectDrawLast
func (w *World_3d) AddObjectDrawLast(obj *Object3d, x, y, z float64) {
	w.objectToDrawLast = append(w.objectToDrawLast, obj)
	w.objectToDrawLastXpos = append(w.objectToDrawLastXpos, x)
	w.objectToDrawLastYpos = append(w.objectToDrawLastYpos, y)
	w.objectToDrawLastZpos = append(w.objectToDrawLastZpos, z)
}

func (w *World_3d) AddCamera(c *Camera, x, y, z float64) {
	w.cameras = append(w.cameras, c)
	w.camXpos = append(w.camXpos, x)
	w.camYpos = append(w.camYpos, y)
	w.camZpos = append(w.camZpos, z)
	w.currentCamera = len(w.cameras) - 1
}

func paint(screen *ebiten.Image, xsize, ysize int, obj *Object3d, x, y, z float64, cam *Camera) {

	objToWorld := TransMatrix(x, y, z)

	objToCam := cam.camMatrixRev.MultiplyBy(objToWorld)
	obj.ApplyMatrixTemp(objToCam)
	obj.PaintObject(screen, xsize/2, ysize/2, true)
}

func distBetweenObjectAndCamera(obj *Object3d, cam *Camera) float64 {
	objX, objY, objZ := obj.GetPosition().X, obj.GetPosition().Y, obj.GetPosition().Z
	camX, camY, camZ := cam.GetPosition().X, cam.GetPosition().Y, cam.GetPosition().Z

	return math.Sqrt(math.Pow(objX-camX, 2) + math.Pow(objY-camY, 2) + math.Pow(objZ-camZ, 2))
}

func draw(screen *ebiten.Image, xsize, ysize int, objects []*Object3d, cam *Camera) {
	// draw background objects
	for _, obj := range objects {
		objToWorld := TransMatrix(
			obj.GetPosition().X,
			obj.GetPosition().Y,
			obj.GetPosition().Z,
		)
		objToCam := cam.camMatrixRev.MultiplyBy(objToWorld)
		obj.ApplyMatrixTemp(objToCam)
		obj.PaintObject(screen, xsize/2, ysize/2, true)
	}

}

func sortObjects(backgroundObjects []*Object3d, cam *Camera) {
	sort.Slice(backgroundObjects, func(i, j int) bool {
		distanceI := distBetweenObjectAndCamera(backgroundObjects[i], cam)
		distanceJ := distBetweenObjectAndCamera(backgroundObjects[j], cam)
		return distanceI > distanceJ
	})
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
	sortFunc := func(i, j int) bool {
		distanceI := math.Sqrt(math.Pow(w.objXpos[sortedIndices[i]]-camX, 2) +
			math.Pow(w.objYpos[sortedIndices[i]]-camY, 2) +
			math.Pow(w.objZpos[sortedIndices[i]]-camZ, 2))

		distanceJ := math.Sqrt(math.Pow(w.objXpos[sortedIndices[j]]-camX, 2) +
			math.Pow(w.objYpos[sortedIndices[j]]-camY, 2) +
			math.Pow(w.objZpos[sortedIndices[j]]-camZ, 2))

		return distanceI > distanceJ
	}
	sort.Slice(sortedIndices, sortFunc)

	// // Draw objects that should be drawn first (that dont have a direction vector)
	for i, obj := range w.objectToDrawFirst {
		if obj.objectDirection != nil {
			continue
		}

		paint(screen, xsize, ysize, obj, w.objectToDrawFirstXpos[i], w.objectToDrawFirstYpos[i], w.objectToDrawFirstZpos[i], cam)
	}

	// get objects which are poing at and from the camera. objects pointing towards the camera are drawn first
	backgroundObjects := make([]*Object3d, 0)
	foregroundObjects := make([]*Object3d, 0)
	for i, obj := range w.objectToDrawFirst {
		if obj.objectDirection == nil {
			continue
		}

		objToWorld := TransMatrix(w.objectToDrawFirstXpos[i], w.objectToDrawFirstYpos[i], w.objectToDrawFirstZpos[i])
		objToCam := cam.camMatrixRev.MultiplyBy(objToWorld)

		// Transform the direction vector (the direction the cusion is pointing) to camera space
		vecMatrix := NewMatrix()
		vecMatrix.AddRow([]float64{obj.objectDirection.X, obj.objectDirection.Y, obj.objectDirection.Z, 1.0})
		destMatrix := IdentMatrix()
		objToCam.TransformNormals(vecMatrix, destMatrix)
		direction := NewVector3(destMatrix.ThisMatrix[0][0], destMatrix.ThisMatrix[0][1], destMatrix.ThisMatrix[0][2])
		direction.Normalize()

		obj.ApplyMatrixTemp(objToCam)
		pointX := obj.transFaceMesh.Points.ThisMatrix[0][0]
		pointY := obj.transFaceMesh.Points.ThisMatrix[0][1]
		pointZ := obj.transFaceMesh.Points.ThisMatrix[0][2]

		plane := NewPlaneFromPoint(NewPoint3d(pointX, pointY, pointZ), direction)
		where := plane.PointOnPlane(0, 0, 0)

		if where > 0 {
			backgroundObjects = append(backgroundObjects, obj)
		} else {
			foregroundObjects = append(foregroundObjects, obj)
		}
	}

	// sort background objects by distance to camera and draw them
	sortObjects(backgroundObjects, cam)
	draw(screen, xsize, ysize, backgroundObjects, cam)

	for _, i := range sortedIndices {
		obj := w.objects[i]

		// object space to world space trans
		objToWorld := TransMatrix(
			w.objXpos[i],
			w.objYpos[i],
			w.objZpos[i],
		)

		// objToCam := objToWorld.MultiplyBy(cam.camMatrixRev)
		objToCam := cam.camMatrixRev.MultiplyBy(objToWorld)

		// then, world space to camera space
		obj.ApplyMatrixTemp(objToCam)
		obj.PaintObject(screen, xsize/2, ysize/2, true)
	}

	// draw foreground objects
	sortObjects(foregroundObjects, cam)
	draw(screen, xsize, ysize, foregroundObjects, cam)

	// Draw objects that should be drawn last
	for i, obj := range w.objectToDrawLast {

		paint(screen, xsize, ysize, obj, w.objectToDrawLastXpos[i], w.objectToDrawLastYpos[i], w.objectToDrawLastZpos[i], cam)
	}
}
