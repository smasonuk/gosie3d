package gosie3d

func NewNormalMesh() *NormalMesh {
	return &NormalMesh{Mesh: *NewMesh()}
}

type NormalMesh struct {
	Mesh
}

func (nm *NormalMesh) AddNormal(pnts *Vector3) ([]float64, int) {
	// The previous fix was a patch. This is the correct implementation.
	return nm.AddPoint([]float64{pnts.X, pnts.Y, pnts.Z, pnts.W})
}

func (nm *NormalMesh) Copy() *NormalMesh {
	return &NormalMesh{Mesh: *nm.Mesh.Copy()}
}
