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

from sys import exit

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
from OCC.gp import gp_Ax1, gp_Pnt, gp_Dir, gp_Trsf, gp_Quaternion, gp_Vec, gp_GTrsf
from OCC.TopLoc import TopLoc_Datum3D, TopLoc_Location
from OCC.gp import gp_Vec, gp_Pnt, gp_Trsf, gp_OX, gp_OY, gp_OZ
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Display.SimpleGui import init_display

# from core_geometry_utils import translate_shp, rotate_shp_3_axis

filename = './models/as1-oc-214.stp'
#filename = './models/as1_pe_203.stp'
filename = './models/cubez.stp'
#filename = './as1-oc-214.stp'
#filename = './as1_pe_203.stp'
#filename = './models/SPKcorrector.stp'

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


def getMatrix1(tran):
    m = tran.VectorialPart()
    gp = m.Row(1)
    x1 = gp.X()
    x2 = gp.Y()
    x3 = gp.Z()
    x4 = tran.Transforms()[0]
    gp = m.Row(2)
    y1 = gp.X()
    y2 = gp.Y()
    y3 = gp.Z()
    y4 = tran.Transforms()[1]
    gp = m.Row(3)
    z1 = gp.X()
    z2 = gp.Y()
    z3 = gp.Z()
    z4 = tran.Transforms()[2]
    return [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4]

def getMatrix(location):
	"""
	 == == == == == = ==
	x1 x2 x3 x4  x = x'
	y1 y2 y3 y4  y = y'
	z1 z2 z3 z4  z = z'
	0  0  0  1   1 = 1
	== == == == == = ==
	:param location: location
	:type location: :class:`.OCC.TopLoc.TopLoc_Location`
	"""
	m = location.VectorialPart()
	gp = m.Row(1)
	x1 = gp.X()
	x2 = gp.Y()
	x3 = gp.Z()
	x4 = location.Transforms()[0]
	gp = m.Row(2)
	y1 = gp.X()
	y2 = gp.Y()
	y3 = gp.Z()
	y4 = location.Transforms()[1]
	gp = m.Row(3)
	z1 = gp.X()
	z2 = gp.Y()
	z3 = gp.Z()
	z4 = location.Transforms()[2]
	the_trsf = gp_Trsf()
	the_trsf.SetValues(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4)

#	print("    transfroms:", location.Transforms())
	print("    matrix:", x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4)
	return the_trsf
	#the_gtrsf = gp_GTrsf()
	#the_gtrsf.SetTrsf(the_trsf)
	#return the_gtrsf

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
	users = TDF_LabelSequence()
	users_cnt = shape_tool.GetUsers(l_shape, users)
	print("Nr Users       :", users_cnt)

	loc = shape_tool.GetLocation(l_shape)
	print("    location     :", loc)
	tran = loc.Transformation()
	print("    transformation     :", tran)
	m = getMatrix(tran)
	print("    location matrix     :", m)

	handled = False

	if shape_tool.IsReference(l_shape):
		handleReference(l_shape)
#		handled = True

	if shape_tool.IsSimpleShape(l_shape):
		handleSimpleShape(l_shape)
		handled = True

	if not handled and shape_tool.IsShape(l_shape):
		if shape_tool.IsAssembly(l_shape):
			handleAssembly(l_shape)
			handled = True
		elif shape_tool.IsCompound(l_shape):
			handleCompound(l_shape)
			handled = True
		elif shape_tool.IsFree(l_shape):
			handleFree(l_shape)
			handled = True

	if not handled and shape_tool.IsComponent(l_shape):
		handleComponent(l_shape)
		handled = True

#	else:
#		print("UNHANDLED")
#		exit(1)

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
	print("handling SIMPLE SHAPE")
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
	colorSet = False
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
		colorSet = True
#		return c

	if not colorSet:
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
#			return c

	l_subss = TDF_LabelSequence()
	r = shape_tool.GetSubShapes(l_shape, l_subss)
	print("    Nb subshapes     :", l_subss.Length())
	print()
#	subs = False
	for j in range(l_subss.Length()):
		print("\nhandling SHAPE SUBSHAPE\n")
		l_subs = l_subss.Value(j+1)
		handleLabel(l_subs)
#		subs = True
#	if subs:
#		return True

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
		l_shape = (i+1)
		handleLabel(l_shape)

def getShapes2():
	l_shapes = TDF_LabelSequence()
	shape_tool.GetShapes(l_shapes)

	print()
	print("Number of shapes :", l_shapes.Length())
	print()
	for i in range(l_shapes.Length()):
		l_shape = l_shapes.Value(i+1)
#		shape = shape_tool.GetShape(l_shape)
#		label = getLabelName(l_shape)
#		print("label: ", label)
		handleLabel(l_shape)


def translate_shp(shp, vec, copy=False):
	trns = gp_Trsf()
	trns.SetTranslation(vec)
	brep_trns = BRepBuilderAPI_Transform(shp, trns, copy)
	brep_trns.Build()
	return brep_trns.Shape()


def rotate_shp_3_axis(shape, rx, ry, rz, unity="deg"):
	""" Rotate a shape around (O,x), (O,y) and (O,z).

	@param rx_degree : rotation around (O,x)
	@param ry_degree : rotation around (O,y)
	@param rz_degree : rotation around (O,z)

	@return : the rotated shape.
	"""
	if unity == "deg":  # convert angle to radians
		rx = radians(rx)
		ry = radians(ry)
		rz = radians(rz)
	alpha = gp_Trsf()
	alpha.SetRotation(gp_OX(), rx)
	beta = gp_Trsf()
	beta.SetRotation(gp_OY(), ry)
	gamma = gp_Trsf()
	gamma.SetRotation(gp_OZ(), rz)
	brep_trns = BRepBuilderAPI_Transform(shape, alpha*beta*gamma, False)
	brep_trns.Build()
	shp = brep_trns.Shape()
	return shp

lvl = 0
locs = []
def getSubShapes3(lab, loc):
	global cnt, lvl
	cnt += 1
	print("\n[%d] level %d, handling LABEL %s\n" % (cnt, lvl, getLabelName(lab)))
	print()
	print(lab.DumpToString())
	print()
	print("Is Assembly    :", shape_tool.IsAssembly(lab))
	print("Is Free        :", shape_tool.IsFree(lab))
	print("Is Shape       :", shape_tool.IsShape(lab))
	print("Is Compound    :", shape_tool.IsCompound(lab))
	print("Is Component   :", shape_tool.IsComponent(lab))
	print("Is SimpleShape :", shape_tool.IsSimpleShape(lab))
	print("Is Reference   :", shape_tool.IsReference(lab))

	users = TDF_LabelSequence()
	users_cnt = shape_tool.GetUsers(lab, users)
	print("Nr Users       :", users_cnt)

	l_subss = TDF_LabelSequence()
	r = shape_tool.GetSubShapes(lab, l_subss)
	print("Nb subshapes   :", l_subss.Length())
	l_comps = TDF_LabelSequence()
	shape_tool.GetComponents(lab, l_comps)
	print("Nb components  :", l_comps.Length())
	print()

	if shape_tool.IsAssembly(lab):
		l_c = TDF_LabelSequence()
		shape_tool.GetComponents(lab, l_c)
		for i in range(l_c.Length()):
			label = l_c.Value(i+1)
			if shape_tool.IsReference(label):
				print("\n########  reference label :", label)
				label_reference=TDF_Label()
				shape_tool.GetReferredShape(label, label_reference)
				loc = shape_tool.GetLocation(label)
				print("    loc          :", loc)
				trans = loc.Transformation()
				print("    tran form    :", trans.Form())
				rot = trans.GetRotation()
				print("    rotation     :", rot)
				print("    X            :", rot.X())
				print("    Y            :", rot.Y())
				print("    Z            :", rot.Z())
				print("    W            :", rot.W())
				tran = trans.TranslationPart()
				print("    translation  :", tran)
				print("    X            :", tran.X())
				print("    Y            :", tran.Y())
				print("    Z            :", tran.Z())

				locs.append(loc)
				print(">>>>")
				lvl += 1
				getSubShapes3(label_reference, loc)
				lvl -= 1
				print("<<<<")
				locs.pop()

	elif shape_tool.IsSimpleShape(lab):
		print("\n########  simpleshape label :", lab)
		shape = shape_tool.GetShape(lab)
		print("    all ass locs   :", locs)

		loc = TopLoc_Location()
		for i in range(len(locs)):
			print("    take loc       :", locs[i])
			loc = loc.Multiplied(locs[i])

		trans = loc.Transformation()
		print("    FINAL loc    :")
		print("    tran form    :", trans.Form())
		rot = trans.GetRotation()
		print("    rotation     :", rot)
		print("    X            :", rot.X())
		print("    Y            :", rot.Y())
		print("    Z            :", rot.Z())
		print("    W            :", rot.W())
		tran = trans.TranslationPart()
		print("    translation  :", tran)
		print("    X            :", tran.X())
		print("    Y            :", tran.Y())
		print("    Z            :", tran.Z())
		shape = BRepBuilderAPI_Transform(shape, loc.Transformation()).Shape()

		c = Quantity.Quantity_Color()
		if (color_tool.GetColor(lab, 0, c) or
			color_tool.GetColor(lab, 1, c) or
			color_tool.GetColor(lab, 2, c)):
			for i in (0, 1, 2):
				color_tool.SetInstanceColor(shape, i, c)

		n = c.Name(c.Red(), c.Green(), c.Blue())
		print('    color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
		cs = display.DisplayColoredShape(shape, c)

		for i in range(l_subss.Length()):
			lab = l_subss.Value(i+1)
			print("\n########  simpleshape label :", lab)
			shape = shape_tool.GetShape(lab)

			c = Quantity.Quantity_Color()
			if (color_tool.GetColor(lab, 0, c) or
				color_tool.GetColor(lab, 1, c) or
				color_tool.GetColor(lab, 2, c)):
				for i in (0, 1, 2):
					color_tool.SetInstanceColor(shape, i, c)

			n = c.Name(c.Red(), c.Green(), c.Blue())
			print('    color Name & RGB: ', c, n, c.Red(), c.Green(), c.Blue())
			cs = display.DisplayColoredShape(shape, c)

def getShapes3():
	labels = TDF_LabelSequence()
	h_shape_tool.GetObject().GetFreeShapes(labels)
	global cnt
	cnt += 1

	print()
	print("Number of shapes at root :", labels.Length())
	print()
	root = labels.Value(1)

	getSubShapes3(root, None)

def run(event=None):
	display.EraseAll()
#	getShapes()
#	getShapes2()
	getShapes3()
	print()
	print("Handled %d labels" % cnt)
	print()
	display.FitAll()

def exit1(event=None):
    exit()

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
	add_function_to_menu('STEP import', exit1)
	start_display()
