##Copyright 2010-2014 Thomas Paviot (tpaviot@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function


from OCC.TCollection import TCollection_ExtendedString

from OCC.TDocStd import Handle_TDocStd_Document
from OCC.XCAFApp import XCAFApp_Application
from OCC.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                         XCAFDoc_DocumentTool_ColorTool,
                         XCAFDoc_DocumentTool_LayerTool,
                         XCAFDoc_DocumentTool_MaterialTool,
                         XCAFDoc_ColorGen,
                         XCAFDoc_ColorSurf,
                         XCAFDoc_ColorCurv)
from OCC.STEPCAFControl import STEPCAFControl_Reader
from OCC.IFSelect import IFSelect_RetDone
from OCC.TDF import TDF_LabelSequence, TDF_Label, TDF_Tool
from OCC.TDataStd import Handle_TDataStd_Name, TDataStd_Name_GetID
from OCC.TCollection import TCollection_ExtendedString, TCollection_AsciiString
from OCC import Quantity
#from OCC import Standard

from OCC.Display.SimpleGui import init_display

filename = './models/as1_pe_203.stp'
#filename = './as1-oc-214.stp'
#filename = './as1_pe_203.stp'
_shapes = []

# create an handle to a document
h_doc = Handle_TDocStd_Document()

# Create the application
app = XCAFApp_Application.GetApplication().GetObject()
app.NewDocument(TCollection_ExtendedString("MDTV-CAF"), h_doc)

# Get root assembly
doc = h_doc.GetObject()
h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(doc.Main())
h_color_tool = XCAFDoc_DocumentTool_ColorTool(doc.Main())
h_layer_tool = XCAFDoc_DocumentTool_LayerTool(doc.Main())
h_mat_tool = XCAFDoc_DocumentTool_MaterialTool(doc.Main())

step_reader = STEPCAFControl_Reader()
step_reader.SetColorMode(True)
step_reader.SetLayerMode(True)
step_reader.SetNameMode(True)
step_reader.SetMatMode(True)

status = step_reader.ReadFile(filename)
if status == IFSelect_RetDone:
    step_reader.Transfer(doc.GetHandle())


color_tool = h_color_tool.GetObject()

def getColors():
	print()
	print("========================================================================")
	print("                          COLORS")

	color_labels = TDF_LabelSequence()
	color_tool.GetColors(color_labels)

	print()
	print("Number of colors:", color_labels.Length())
	print()

	for i in range(color_labels.Length()):
		l_color = color_labels.Value(i+1)
		print(l_color.DumpToString())
		c = Quantity.Quantity_Color()
		r = color_tool.GetColor(l_color, c)
		#print('color table label, color: ', l_color, c, r)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Table color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())

def getLayers():
	print()
	print("========================================================================")
	print("                          LAYERS")

	layer_tool = h_layer_tool.GetObject()
	layer_labels = TDF_LabelSequence()
	layer_tool.GetLayerLabels(layer_labels)

	print()
	print("Number of layers:", layer_labels.Length())
	print()


	for i in range(layer_labels.Length()):
		l_layer = layer_labels.Value(i+1)
		print(l_layer.DumpToString())

def getMaterials():
	print()
	print("========================================================================")
	print("                          MATERIALS")

	mat_tool = h_mat_tool.GetObject()
	mat_labels = TDF_LabelSequence()
	mat_tool.GetMaterialLabels(mat_labels)

	print()
	print("Number of materials:", mat_labels.Length())
	print()


	for i in range(mat_labels.Length()):
		l_mat = mat_labels.Value(i+1)
		print(l_mat.DumpToString())


def getLabelName(lab):
	entry = TCollection_AsciiString()
	TDF_Tool.Entry(lab, entry)
	N = Handle_TDataStd_Name()
	lab.FindAttribute(TDataStd_Name_GetID(), N)
	n = N.GetObject()
	if n:
#		return unicode(n.Get().PrintToString(),errors='ignore')#.decode("latin-1")
		return n.Get().PrintToString()
	return "No Name"

cnt = 0
def handleLabel(l_shape):
	global cnt
	cnt += 1
	print("\n[%d]  handling LABEL %s\n" % (cnt, getLabelName(l_shape)))
	print()
	print(l_shape.DumpToString())
	print()
	print("Is Assembly    :", shape_tool.IsAssembly(l_shape))
	print("Is Free        :", shape_tool.IsFree(l_shape))
	print("Is Shape       :", shape_tool.IsShape(l_shape))
	print("Is Compound    :", shape_tool.IsCompound(l_shape))
	print("Is Component   :", shape_tool.IsComponent(l_shape))
	print("Is SimpleShape :", shape_tool.IsSimpleShape(l_shape))
	print("Is Reference   :", shape_tool.IsReference(l_shape))
#	print("Is SubShape    :", shape_tool.IsSubShape(l_shape))

	if shape_tool.IsShape(l_shape):
		if shape_tool.IsAssembly(l_shape):
			handleAssembly(l_shape)
		elif shape_tool.IsCompound(l_shape):
			handleCompound(l_shape)
		elif shape_tool.IsFree(l_shape):
			handleFree(l_shape)

	if shape_tool.IsComponent(l_shape):
		handleComponent(l_shape)

	if shape_tool.IsSimpleShape(l_shape):
		handleSimpleShape(l_shape)
#	else:
#		print("UNHANDLED")
#		exit(1)

	if shape_tool.IsReference(l_shape):
		handleReference(l_shape)

def handleAssembly(l_shape):
	print("\nhandling ASSEMBLY\n")
	l_comps = TDF_LabelSequence()
	l_subss = TDF_LabelSequence()

	r = shape_tool.GetSubShapes(l_shape, l_subss)
	print("    Nb subshapes     :", l_subss.Length())
	print()
	for j in range(l_subss.Length()):
		print("\nhandling ASSEMBLY SUBSHAPE\n")
		l_subs = l_subss.Value(j+1)
		handleLabel(l_subs)

	print("    Nb comp of shape :", shape_tool.NbComponents(l_shape))
	r = shape_tool.GetComponents(l_shape, l_comps)
	print("    Nb comp loaded   :", l_comps.Length())
	print()
	for j in range(l_comps.Length()):
		l_comp = l_comps.Value(j+1)
		handleLabel(l_comp)

def handleCompound(l_shape):
	print("\nhandling COMPOUND\n")
	l_comps = TDF_LabelSequence()
	l_subss = TDF_LabelSequence()

	r = shape_tool.GetSubShapes(l_shape, l_subss)
	print("    Nb subshapes     :", l_subss.Length())
	print()
	for j in range(l_subss.Length()):
		print("\nhandling COMPOUND SUBSHAPE\n")
		l_subs = l_subss.Value(j+1)
#		print(l_subs.DumpToString())
		handleLabel(l_subs)

	print("    Nb comp of shape :", shape_tool.NbComponents(l_shape))
	r = shape_tool.GetComponents(l_shape, l_comps)
	print("    Nb comp loaded   :", l_comps.Length())
	print()
	for j in range(l_comps.Length()):
		l_comp = l_comps.Value(j+1)
		handleLabel(l_comp)

def handleFree(l_shape):
	print("\nhandling FREE\n")
	l_comps = TDF_LabelSequence()
	l_subss = TDF_LabelSequence()

	r = shape_tool.GetSubShapes(l_shape, l_subss)
	print("    Nb subshapes     :", l_subss.Length())
	print()
	for j in range(l_subss.Length()):
		print("\nhandling FREE SUBSHAPE\n")
		l_subs = l_subss.Value(j+1)
		handleLabel(l_subs)

	print("    Nb comp of shape :", shape_tool.NbComponents(l_shape))
	r = shape_tool.GetComponents(l_shape, l_comps)
	print("    Nb comp loaded   :", l_comps.Length())
	print()
	for j in range(l_comps.Length()):
		l_comp = l_comps.Value(j+1)
		handleLabel(l_comp)

def handleSimpleShape(l_shape):
	print("handling SHAPE")
	cg = color_tool.IsSet(l_shape, XCAFDoc_ColorGen)
	cs = color_tool.IsSet(l_shape, XCAFDoc_ColorSurf)
	cc = color_tool.IsSet(l_shape, XCAFDoc_ColorCurv)
	print("color_tool.IsSet: Gen=", cg, "Surf=", cs, "Curv=", cc)

	shape = shape_tool.GetShape(l_shape)

	if cg:
		c = Quantity.Quantity_Color()
		r = color_tool.GetColor(l_shape, XCAFDoc_ColorGen, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Gen color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
		c = Quantity.Quantity_Color()
		r = color_tool.GetInstanceColor(shape, XCAFDoc_ColorGen, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Gen instance color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())


	if cs:
		c = Quantity.Quantity_Color()
		r = color_tool.GetColor(l_shape, XCAFDoc_ColorSurf, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Surf color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
		c = Quantity.Quantity_Color()
		r = color_tool.GetInstanceColor(shape, XCAFDoc_ColorGen, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Surf instance color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())

	if cc:
		c = Quantity.Quantity_Color()
		r = color_tool.GetColor(l_shape, XCAFDoc_ColorCurv, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Curv color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
		c = Quantity.Quantity_Color()
		r = color_tool.GetInstanceColor(shape, XCAFDoc_ColorGen, c)
		if r:
			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('Curv instance color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())

#	_shapes.append(shape)
	c = Quantity.Quantity_Color()
	if (color_tool.GetInstanceColor(shape, 0, c) or
		color_tool.GetInstanceColor(shape, 1, c) or
		color_tool.GetInstanceColor(shape, 2, c)):
		for i in (0, 1, 2):
			color_tool.SetInstanceColor(shape, i, c)

		n = c.Name(c.Red(), c.Green(), c.Blue())
		print('1 set instance color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
#		_shapes.append(shape)
		display.DisplayColoredShape(shape, c)
		return c

	if (color_tool.GetColor(l_shape, 0, c) or
		color_tool.GetColor(l_shape, 1, c) or
		color_tool.GetColor(l_shape, 2, c)):
		#shape = shape_tool.GetShape(label)
		for i in (0, 1, 2):
			color_tool.SetInstanceColor(shape, i, c)

		n = c.Name(c.Red(), c.Green(), c.Blue())
		print('2 set instance color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
#		_shapes.append(shape)
		display.DisplayColoredShape(shape, c)
		return c

def handleComponent(l_shape):
	print("handling COMPONENT")
#	exit(1)
	handleSimpleShape(l_shape)

def handleReference(l_shape):
	print("handling REFERENCE")
	l_ref = TDF_Label()
	r = shape_tool.GetReferredShape(l_shape, l_ref)
	print("    refered shape  :", r, l_ref)
	handleLabel(l_ref)

print()
print("========================================================================")
print("                          SHAPES")

shape_tool = h_shape_tool.GetObject()
shape_tool.SetAutoNaming(True)

def getShapes():
	labels = TDF_LabelSequence()
	h_shape_tool.GetObject().GetFreeShapes(labels)

	print()
	print("Number of shapes at root :", labels.Length())
	print()
	cnt = 0
	for i in range(labels.Length()):
		l_shape = labels.Value(i+1)
		handleLabel(l_shape)

def getShapes2():
	l_shapes = TDF_LabelSequence()
	shape_tool.GetShapes(l_shapes)
	for i in range(l_shapes.Length()):
		l_shape = l_shapes.Value(i+1)
		shape = shape_tool.GetShape(l_shape)
		label = getLabelName(l_shape)
		print("label: ", label)

def run(event=None):
	display.EraseAll()
	getShapes()
	print()
	print("Handled %d labels" % cnt)
	print()
	display.FitAll()

def exit(event=None):
    sys.exit()

#getColors()

#getShapes2()

#print()
#print("Handled %d labels" % cnt)
#print()
#print("Have %d shapes" % len(_shapes))

#print()
#print()
#print("========================================================================")
#print("                          DRAW")

#for i in range(labels.Length()):
#	l_shape = labels.Value(i+1)
#	a_shape = shape_tool.GetShape(l_shape)
#	_shapes.append(a_shape)

#
# Display
#
#display, start_display, add_menu, add_function_to_menu = init_display()
#display.DisplayShape(_shapes, update=True)
#start_display()

if __name__ == '__main__':
	display, start_display, add_menu, add_function_to_menu = init_display()
	add_menu('STEP import')
	add_function_to_menu('STEP import', run)
	start_display()

