package gosie3d

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
