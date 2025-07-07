package gosie3d

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
