package gosie3d

import (
	"math"
	"sort"

	"github.com/hajimehoshi/ebiten/v2"
)

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

	// // Draw objects that should be drawn first
	// for _, obj := range w.objectToDrawFirst {
	// 	// m := TransMatrix(w.objectToDrawFirstXpos[i]-camX, w.objectToDrawFirstYpos[i]-camY, w.objectToDrawFirstZpos[i]-camZ)
	// 	// obj.ApplyMatrixTemp(cam.camMatrixRev.MultiplyBy(m))
	// 	obj.ApplyMatrixTemp(cam.camMatrixRev)
	// 	obj.PaintSolid(screen, xsize/2, ysize/2, true)
	// }

	for _, i := range sortedIndices {
		obj := w.objects[i]
		// m := TransMatrix(w.objXpos[i]-camX, w.objYpos[i]-camY, w.objZpos[i]-camZ)
		// obj.ApplyMatrixTemp(cam.camMatrixRev.MultiplyBy(m))

		// object space to world space trans
		objToWorld := TransMatrix(
			w.objXpos[i],
			w.objYpos[i],
			w.objZpos[i],
		)

		// objToCam := objToWorld.MultiplyBy(cam.camMatrixRev)
		objToCam := cam.camMatrixRev.MultiplyBy(objToWorld)

		// then, world space to camera space
		// obj.ApplyMatrixTemp(cam.camMatrixRev)
		obj.ApplyMatrixTemp(objToCam)
		obj.PaintSolid(screen, xsize/2, ysize/2, true)
	}
}
