package gosie3d

func NewNormalMesh() *NormalMesh {
	return &NormalMesh{Mesh: *NewMesh()}
}

type NormalMesh struct {
	Mesh
}

func (nm *NormalMesh) AddNormal(pnts []float64) ([]float64, int) {
	// The previous fix was a patch. This is the correct implementation.
	return nm.AddPoint(pnts)
}

func (nm *NormalMesh) Copy() *NormalMesh {
	return &NormalMesh{Mesh: *nm.Mesh.Copy()}
}
