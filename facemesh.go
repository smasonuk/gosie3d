package gosie3d

func NewFaceMesh() *FaceMesh {
	return &FaceMesh{Mesh: *NewMesh()}
}

type FaceMesh struct {
	Mesh
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
