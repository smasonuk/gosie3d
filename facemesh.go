package gosie3d

type FaceMesh struct {
	Mesh
	// faces []*Face
}

func NewFaceMesh() *FaceMesh {
	return &FaceMesh{Mesh: *NewMesh()}
}

// AddFace adds a face to the FaceMesh and returns the new face and its indices.
func (fm *FaceMesh) AddFace(f *Face) (*Face, []int) { // Return indices
	newPoints := make([][]float64, len(f.Points))
	indices := make([]int, len(f.Points))
	for i, p := range f.Points {
		newPoints[i], indices[i] = fm.AddPoint(p)
	}
	newface := NewFace(newPoints, f.Col, f.GetNormal())

	// fm.faces = append(fm.faces, newface)

	return newface, indices
}

func (fm *FaceMesh) Copy() *FaceMesh {
	m := &FaceMesh{
		Mesh: *fm.Mesh.Copy(),
		// faces: make([]*Face, len(fm.faces)),
	}

	// for i, f := range fm.faces {
	// 	m.faces[i] = f.Copy()
	// }

	return m
}
