package gosie3d

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
