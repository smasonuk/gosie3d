package gosie3d

type Mesh struct {
	Points     *Matrix
	pointIndex map[[3]float64]int
}

func NewMesh() *Mesh {
	return &Mesh{
		Points:     NewMatrix(),
		pointIndex: make(map[[3]float64]int),
	}

}

// func (m *Mesh) AddPoint(point []float64) ([]float64, int) {
// 	// Check if the point already exists in the mesh and return it if found
// 	for i, p := range m.Points.ThisMatrix {
// 		if p[0] == point[0] && p[1] == point[1] && p[2] == point[2] {
// 			return p, i
// 		}
// 	}

// 	pointCopy := make([]float64, len(point))
// 	copy(pointCopy, point)
// 	m.Points.AddRow(pointCopy)
// 	newIndex := len(m.Points.ThisMatrix) - 1
// 	return pointCopy, newIndex
// }

// func (m *Mesh) Copy() *Mesh {
// 	return &Mesh{
// 		Points: m.Points.Copy(),
// 	}
// }

// AddPoint now uses the map for an average O(1) lookup.
func (m *Mesh) AddPoint(point []float64) ([]float64, int) {
	// Create a comparable key from the point slice.
	// Assumes point slice will always have 3 elements.
	pointKey := [3]float64{point[0], point[1], point[2]}

	// Check if the point already exists using the map.
	if index, found := m.pointIndex[pointKey]; found {
		return m.Points.ThisMatrix[index], index
	}

	// If not found, add the new point.
	pointCopy := make([]float64, len(point))
	copy(pointCopy, point)
	m.Points.AddRow(pointCopy)
	newIndex := len(m.Points.ThisMatrix) - 1

	// Add the new point's key and index to the map.
	m.pointIndex[pointKey] = newIndex

	return pointCopy, newIndex
}

// Copy must also duplicate the pointIndex map.
func (m *Mesh) Copy() *Mesh {
	newPointIndex := make(map[[3]float64]int, len(m.pointIndex))
	for key, value := range m.pointIndex {
		newPointIndex[key] = value
	}

	return &Mesh{
		Points:     m.Points.Copy(),
		pointIndex: newPointIndex,
	}
}
