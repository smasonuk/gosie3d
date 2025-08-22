package gosie3d

import "image/color"

// PolygonBatcher is a mock for testing purposes
type PolygonBatcher struct{}

func (b *PolygonBatcher) AddPolygon(xp, yp []float32, clr color.RGBA) {}

func (b *PolygonBatcher) AddPolygonAndOutline(xp, yp []float32, fillClr, strokeClr color.RGBA, strokeWidth float32) {}
