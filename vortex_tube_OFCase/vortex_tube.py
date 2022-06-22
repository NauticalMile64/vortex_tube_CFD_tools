# This script generates the vortex tube meshes for the paper entitled:
# 	"The The Impact of Boundary Treatment and Turbulence Model on CFD Simulations of the Ranque-Hilsch Vortex Tube"

# Running this script overwrites the blockMeshDict file in the system folder
# The blockmeshbuilder module is available from https://github.com/NauticalMile64/blockmeshbuilder

import blockmeshbuilder as bmh
import numpy as np

# Settings
blocksOnly = False  # Each block written to file contains a single cell
GradingFlag = True  # Turn Grading on/off
Wedge = False  # [Not Working] Only write the blocks for 0 -> pi/2
Twist = False  # [Not tested with release version] Twist the mesh to accommodate a perpendicular inlet shroud
IsInlets = True  # If the inlets (vortex generator) should be modelled, or the inlets should be flush with the tube
RibbonInlet = False  # If the inlet should be modelled as a ribbon (3D vortex tube, but geometry is axis-symmetric)
ExtCone = False  # If the cone block is to be meshed by an external program
IsPlenums = False  # Should the plenums be included
Axisym = False  # Should the vortex tube be modelled using axis-symmetry (model as wedge with periodic boundaries)
SmallYPlus = True  # Yplus ~ 1 based on wall speed of 328 m/s
DENSITY_SCALE = 0.8  # 1.4 is default	# Changes the overall density of the VT mesh
LENGTH_SCALE = 1.0  # Can use this to adjust the length of the VT while keeping the mesh density roughly the same

# Dimensions (mm)
rad = 12.52 / 2
workingLength = 313.06 * LENGTH_SCALE

inletLength = 4.5
inletWidth = 1.5 if not (Axisym or RibbonInlet) else 0.68635
inletHeight = 1.5

shroudRad = 19.4469 / 2

coldTubeRad = 2.5
coldTubeOuterRad = 4.5
coldTubeLen = 22.0
coldPlenumLen = 25
coldPlenumOffset = 11.69
coldPlenumTubeRad = 17.5

hotAnnulusOffset = 1.0
coneOffset = 2.985584
hotPlenumLen = 8.18
hotPlenumRad = 9.725

cone_angle = np.radians(30)
cone_slope = np.tan(cone_angle)
nose_drad = 2.5
rear_drad = 7.407477
cone_dist = 1

cone_offset_point = workingLength + cone_dist
cone_rear_z = cone_offset_point + (rear_drad - rad) / cone_slope
cone_tip_z = cone_offset_point + (nose_drad - rad) / cone_slope


def get_cone_z(r):
	return cone_tip_z + (r - nose_drad) / cone_slope


def get_cone_r(z):
	return nose_drad + (z - cone_tip_z) * cone_slope


dz = 0.25
print(f'cone tip z = {cone_tip_z + dz}')
print(f'cone tip r = {get_cone_r(cone_tip_z + dz)}')

print(f'cone rear z = {cone_rear_z}')
print(f'cone rear r = {rear_drad}')

plenumExitTubeCoord = 76.7
petc = plenumExitTubeCoord

wl = workingLength
iw = inletWidth
hao = hotAnnulusOffset
ctl = coldTubeLen
cpl = coldPlenumLen
hpl = hotPlenumLen
cpo = coldPlenumOffset

wlo = wl + hao

hpr = hotPlenumRad
cof = coneOffset
ctr = coldTubeRad
ctor = coldTubeOuterRad

hotExitTubeRad = 2.5
hotExitTubeZ = wl + 5
hetr = hotExitTubeRad
hetz = hotExitTubeZ

coldExitTubeRad = 2.5
coldExitTubeZ = -(cpo + 20)
coldExitTubeY = 10.0
cetr = coldExitTubeRad
cetz = coldExitTubeZ
cety = coldExitTubeY

# Radial Divisions
rs = np.array([0.2, 0.9 * ctr / rad, ctr / rad, 1.0, hpr / rad]) * rad
ndr = np.array([3, 1, 30 if SmallYPlus else 13, 4])

if Axisym:
	rs[0] = 0.0
	ndr[0] = 9

# Angular Divisions
if Axisym:
	WA = np.radians(1.)
	ts = np.array([-WA / 2, WA / 2]) + np.pi / 2
	ndt = np.array([1])
else:
	Nt = 9
	ND_t = 7

	ts = np.linspace(0, 2 * np.pi, Nt, endpoint=True)
	if not Twist: ts -= np.pi / 8
	ndt = np.full_like(ts, ND_t, dtype=int)

# Axial Divisions
nIW = 7 if Axisym else 10

zs = np.array([-(cpo + cpl), cetz - 0.75 * cetr, cetz + 0.75 * cetr, -ctl, -cpo, 0.0,
			   iw, get_cone_z(nose_drad), wl,
			   wl + 1, get_cone_z(rear_drad),
			   hetz + hetr * 0.8, wl + hpl])

ndz = np.array([7, 6, 11, 10, 35, nIW, int(round(290 * LENGTH_SCALE)), 27, 10, 9, 15, 8])

# Calculate indicies of important locations
inInd = np.where(zs == 0.0)[0][0]
print('Inlet Plane Index -> inInd = {}'.format(inInd))

pInd = np.where(zs == get_cone_z(rear_drad))[0][0]
print('Phantom Index -> pInd = {}'.format(pInd))

csInd = np.where(zs == get_cone_z(nose_drad))[0][0]
print('Cone Start Index -> csInd = {}'.format(csInd))

ceInd = np.where(zs == get_cone_z(rear_drad))[0][0]
print('Cone End Index -> ceInd = {}'.format(ceInd))

ctInd = np.where(rs == coldTubeRad)[0][0]
print('Cold Tube Radius Index -> ctInd = {}'.format(ctInd))

wrInd = np.where(rs == rad)[0][0]
print('Working Region Radius Index -> wrInd = {}'.format(wrInd))

plInd = np.where(rs == hotPlenumRad)[0][0]
print('Hot plenum radius Index -> plInd = {}'.format(plInd))

cplInd = plInd
print('Cold plenum radius Index -> cplInd = {}'.format(cplInd))

crInd = ctInd
print('Cone radius Index -> crInd = {}'.format(crInd))

if Axisym:
	tube = bmh.TubeBlockStruct(rs, ts, zs, ndr, ndt, ndz, zone='ts')
else:
	vortex_tube = bmh.CylBlockStructContainer(rs, ts, zs, ndr, ndt, ndz, zone='ts', inner_arc_curve=0.0,
											  is_core_aligned=True)

	core = vortex_tube.core_struct
	c_bmask = core.block_mask
	c_vts = core.vertices
	c_divs = core.num_divisions
	c_grad = core.grading

	tube = vortex_tube.tube_struct

t_bmask = tube.block_mask
t_vmask = tube.vertex_mask
t_emask = tube.edge_mask
t_fmask = tube.face_mask

t_vts = tube.vertices
t_divs = tube.num_divisions
t_grad = tube.grading
t_bvts = tube.baked_vertices
t_edges = tube.edges

# O-grid
if not Axisym:
	XS, YS = c_vts[:, :, :pInd, 0], c_vts[:, :, :pInd, 1]
	RS = np.sqrt(XS ** 2 + YS ** 2)
	TS = np.arctan2(YS, XS)
	RS[[0, 1, 1, 2], [1, 0, 2, 1]] *= 1.2
	if not Twist: TS[:] -= np.pi / 8
	XS[:] = RS * np.cos(TS)
	YS[:] = RS * np.sin(TS)

# #################
# Geometry Tweaking
# #################

t_bmask[wrInd:, :, inInd:csInd + 1] = True
t_bmask[ctInd:, :, inInd - 1] = True

# #############
# Hot Exit Tube
# #############

if not Axisym:
	htp1 = bmh.Point(np.array([-1e3, 0.0, hetz]))
	htp2 = bmh.Point(np.array([1e3, 0.0, hetz]))
	hotTubeCyl = bmh.Cylinder(htp1, htp2, hetr, 'hot_exit_tube')

	for pt in t_bvts[-1, 3:5, -3:-1].flatten():
		pt.proj_geom(hotTubeCyl)
	for edge in t_edges[-1, 3:5, -3, 2]:
		edge.proj_geom(hotTubeCyl)
	for edge in t_edges[-1, 3, -3:-1, 1]:
		edge.proj_geom(hotTubeCyl)

	hsc = 0.8
	xht = np.array([hpr + 2, hpr * 2.75, petc])
	yht = np.array([-hetr * hsc, hetr * hsc])
	zht = np.array([hetz - hetr, hetz + hetr])

	nxht = np.array([30, 40])
	nyht = np.full_like(yht, ND_t, dtype=int)
	nzht = np.full_like(zht, ndz[-2])

	hotTube = bmh.CartBlockStruct(xht, yht, zht, nxht, nyht, nzht, zone='hot_tube')

	ht_faces = hotTube.faces
	ht_edges = hotTube.edges
	ht_bvts = hotTube.baked_vertices

	for vt in ht_bvts[1:].flatten():
		vt.proj_geom(hotTubeCyl)
	for edge in ht_edges[1:, ..., 1:].flatten():
		if edge:
			edge.proj_geom(hotTubeCyl)
	for edge in ht_edges[..., 0].flatten():
		if edge:
			edge.proj_geom(hotTubeCyl)
	for face in ht_faces[..., 1:].flatten():
		if face:
			face.proj_geom(hotTubeCyl)

	ht_bvts[0] = t_bvts[-1, 3:5, -3:-1]

	t_fmask[wrInd + 1:, [2, 5], -3:, 0] = False

	t_vts[-1, 3:5, -2, 2] += 0.35
	t_vts[-1, 3:5, -3, 2] -= 0.35

	t_grad[:, 1, -4:, 1] = bmh.SimpleGradingElement(2.)
	t_grad[:, -2, -4:, 1] = bmh.SimpleGradingElement(1 / 2.)

	if Twist:
		t_vts[:, :-1, pInd:, 1] -= np.pi / 32

# ######
# Inlets
# ######

if not Axisym or RibbonInlet:
	t_bmask[wrInd + 1, :, inInd] = False

	if IsInlets:
		t_bmask[wrInd, ::2, inInd] = False

	vin_top = bmh.Point([t_vts[wrInd, 0, inInd, 0], inletLength, 0])
	vin_bot = vin_top - bmh.Point([inletHeight, 0, 0])

	rs_top = t_vts[wrInd + 1, ::2, inInd:(inInd + 2), 0]
	ts_top = t_vts[wrInd + 1, ::2, inInd:(inInd + 2), 1]

	rs_bot = t_vts[wrInd + 1, 1::2, inInd:(inInd + 2), 0]
	ts_bot = t_vts[wrInd + 1, 1::2, inInd:(inInd + 2), 1]

	t_bot = np.arctan2(vin_bot[1], vin_bot[0])
	r_bot = vin_bot[0] / np.cos(t_bot)

	t_top = np.arctan2(vin_top[1], vin_top[0])
	r_top = vin_top[0] / np.cos(t_top)

	ts_top += t_top
	rs_top[:] = r_top

	ts_bot += t_top
	ts_bot[:] = ts_top[:-1] + (t_bot - t_top)
	rs_bot[:] = r_bot

	t_vts[wrInd, 1::2, inInd:(inInd + 2), 1] = t_vts[wrInd, ::2, inInd:(inInd + 2), 1][:-1] + np.arccos(
		vin_bot[0] / t_vts[wrInd, 1::2, inInd:(inInd + 2), 0])

	t_emask[wrInd + 1, :, inInd:(inInd + 2), 1] = True

	t_divs[wrInd, :, inInd:(inInd + 2), 0] = 6
	t_grad[wrInd, :, inInd:(inInd + 2), 0] = bmh.SimpleGradingElement(1 / 3)

	if IsInlets:
		t_grad[wrInd, ::2, inInd:(inInd + 2), 1] = bmh.SimpleGradingElement(1 / 4)

		if SmallYPlus:
			t_divs[wrInd, :, inInd:inInd + 2, 0] *= 4
			t_grad[wrInd, :, inInd:inInd + 2, 0] = bmh.SimpleGradingElement(4.)

# ####
# Cone
# ####

# Mask the necessary blocks
t_bmask[:crInd, ..., csInd:] = True
t_bmask[crInd - 1, ..., ceInd:] = False

if Axisym:
	t_bmask[0, ..., csInd:] = True
else:
	c_bmask[..., csInd:] = True

# If an external cone mesh is generated, mask the cone and change the divisions
if ExtCone:
	t_bmask[crInd:, :, csInd] = True
	t_divs[crInd:, :, csInd + 1:, 0] = round(6 * SCALE)

# Position the vertices
if not ExtCone and not Axisym: t_vts[wrInd, :, csInd, 2] -= 0.35

t_vts[:crInd + 1, :, csInd + 1, 2] += 0.4 * cone_dist / 2

t_vts[wrInd:, :, ceInd, 2] -= 0.2
t_vts[crInd, :, csInd + 1:ceInd + 1, 0] = get_cone_r(t_vts[crInd, :, csInd + 1:ceInd + 1, 2])
t_vts[crInd + 1, :, ceInd, 0] = t_vts[crInd, :, ceInd, 0] + 0.8
t_vts[crInd:crInd + 3, :, ceInd + 1:, 0] = t_vts[crInd:crInd + 3, :, ceInd, 0, np.newaxis]

t_vts[crInd:, :, ceInd - 1] = (t_vts[crInd:, :, ceInd - 2] + t_vts[crInd:, :, ceInd]) / 2
t_vts[plInd, :, ceInd - 1, 2] -= 0.25
t_vts[crInd, :, ceInd - 1, 0] = get_cone_r(t_vts[crInd, :, ceInd - 1, 2])

# Vertices past cone
t_vts[crInd, :, ceInd + 1:, 0] -= 0.5 if ExtCone else 1.25
t_vts[crInd + 1, :, ceInd + 1:, 0] -= 0.5 if ExtCone else 0.25

# Cone start
t_grad[crInd:, :, csInd, 2] = bmh.SimpleGradingElement(1. / 4)

# Cone downstream of working section
t_grad[crInd:, :, csInd + 1, 2] = bmh.SimpleGradingElement(2.5)

t_grad[wrInd, :, csInd + 1, 0] = bmh.SimpleGradingElement(10)
t_grad[wrInd, :, csInd + 2, 0] = bmh.SimpleGradingElement(5)
t_grad[wrInd, :, ceInd, 0] = bmh.SimpleGradingElement(2.5)
t_divs[wrInd, :, csInd + 1:, 0] = 15

# Create splines along the rear portion of the cone
for i in range(ts.size - 1):

	# Cone axial
	spl_vts = t_vts[crInd + 1, i, csInd + 1:csInd + 3]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt = bmh.Point(spl_vts[0] + 0.4 * spl_vec, bmh.cyl_conv_pair)
	spl_pt[0] -= 0.125
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd + 1, i, csInd + 1:csInd + 3], [spl_pt])
	t_edges[crInd + 1, i, csInd + 1, 2] = spl_edge

	spl_vts = t_vts[crInd + 1, i, csInd + 2:csInd + 4]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt = bmh.Point(spl_vts[0] + 0.6 * spl_vec, bmh.cyl_conv_pair)
	spl_pt[0] += 0.125
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd + 1, i, csInd + 2:csInd + 4], [spl_pt])
	t_edges[crInd + 1, i, csInd + 2, 2] = spl_edge

	# Cone radial
	if not ExtCone and IsPlenums:
		spl_vts = t_vts[crInd:crInd + 2, i, csInd]
		spl_vec = spl_vts[1] - spl_vts[0]
		spl_pt = bmh.Point(spl_vts[0] + 0.4 * spl_vec, bmh.cyl_conv_pair)
		spl_pt[0] -= 0.1
		spl_pt[2] -= 0.1
		spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd:crInd + 2, i, csInd], [spl_pt])
		t_edges[crInd, i, csInd, 0] = spl_edge

		spl_vts = t_vts[crInd:crInd + 2, i, csInd + 1]
		spl_vec = spl_vts[1] - spl_vts[0]
		spl_pt = bmh.Point(spl_vts[0] + 0.6 * spl_vec, bmh.cyl_conv_pair)
		spl_pt[2] -= 0.125 * cone_dist / 2
		spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd:crInd + 2, i, csInd + 1], [spl_pt])
		t_edges[crInd, i, csInd + 1, 0] = spl_edge

	spl_vts = t_vts[crInd:crInd + 2, i, csInd + 2]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt = bmh.Point(spl_vts[0] + 0.4 * spl_vec, bmh.cyl_conv_pair)
	spl_pt[2] += 0.125
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd:crInd + 2, i, csInd + 2], [spl_pt])
	# t_edges[crInd,i,csInd+2,0] = spl_edge

	spl_vts = t_vts[crInd + 1:crInd + 3, i, csInd + 2]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt = bmh.Point(spl_vts[0] + 0.2 * spl_vec, bmh.cyl_conv_pair)
	spl_pt[2] -= 0.125
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd + 1:crInd + 3, i, csInd + 2], [spl_pt])
	t_edges[crInd + 1, i, csInd + 2, 0] = spl_edge

	# Behind cone axial
	spl_vts = t_vts[crInd, i, ceInd:ceInd + 2]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt1 = bmh.Point(spl_vts[0] + 0.25 * spl_vec, bmh.cyl_conv_pair)
	spl_pt1[0] += 0.15
	spl_pt2 = bmh.Point(spl_vts[0] + 0.75 * spl_vec, bmh.cyl_conv_pair)
	spl_pt2[0] -= 0.15
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd, i, ceInd:ceInd + 2], [spl_pt1, spl_pt2])
	t_edges[crInd, i, ceInd, 2] = spl_edge

	spl_vts = t_vts[crInd + 1, i, ceInd:ceInd + 2]
	spl_vec = spl_vts[1] - spl_vts[0]
	spl_pt1 = bmh.Point(spl_vts[0] + 0.25 * spl_vec, bmh.cyl_conv_pair)
	spl_pt1[0] += 0.1
	spl_pt2 = bmh.Point(spl_vts[0] + 0.75 * spl_vec, bmh.cyl_conv_pair)
	spl_pt2[0] -= 0.1
	spl_edge = bmh.SplineCurvedEdge(t_bvts[crInd + 1, i, ceInd:ceInd + 2], [spl_pt1, spl_pt2])
	t_edges[crInd + 1, i, ceInd, 2] = spl_edge

# Past the cone
t_divs[crInd - 1, :, ceInd:, 0] = 18
t_grad[crInd - 1, :, ceInd, 0] = bmh.SimpleGradingElement(1 / 10)
t_grad[crInd - 1, :, ceInd + 1:, 0] = bmh.SimpleGradingElement(1 / 3)
t_grad[crInd, :, ceInd, 2] = bmh.SimpleGradingElement(2)

# #########
# Cold Exit
# #########

t_vts[cplInd, :, :inInd, 0] = coldPlenumTubeRad
t_vts[ctInd + 1, :, :inInd, 0] = coldTubeOuterRad

t_bmask[ctInd, :, inInd - 2] = True

t_divs[ctInd, :, :inInd, 0] = 6
t_divs[cplInd - 1, :, :inInd, 0] = 16

len_pcts = np.array([0.2, 0.6, 0.2])
dens = np.array([3., 1., 1., 3.5])
grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))
t_grad[cplInd - 1, :, :inInd, 0] = grd_elm

# Cold Exit Tube
if not Axisym:
	ctp1 = bmh.Point(np.array([-1e3, cety, cetz]))
	ctp2 = bmh.Point(np.array([1e3, cety, cetz]))
	coldTubeCyl = bmh.Cylinder(ctp1, ctp2, cetr, name='cold_exit_tube')

	t_vts[cplInd - 1:, :4, :4, 1] -= np.pi / 16

	if Twist:
		t_vts[2:, :-1, :inInd, 1] -= np.pi / 16
		t_vts[3:, :-1, :inInd, 1] -= np.pi / 32
		t_vts[4:, :-1, :inInd, 1] -= np.pi / 32

	t_vts[cplInd - 1:, 4, 1:3, 1] += np.pi / 32
	t_vts[cplInd - 1:, 5, 1:3, 1] -= np.pi / 8

	t_vts[cplInd - 1:, 4, [0, 3], 1] += np.pi / 32
	t_vts[cplInd - 1:, 5, [0, 3], 1] -= np.pi / 32
	t_vts[cplInd - 1:, 4:6, [0, 3], 1] -= np.pi / 32

	t_grad[cplInd - 1:, 3, :3, 1] = bmh.SimpleGradingElement(1 / 2.)
	t_grad[cplInd - 1:, 5, :3, 1] = bmh.SimpleGradingElement(2.)

	for vt in t_bvts[-1, 4:6, 1:3].flatten():
		vt.proj_geom(coldTubeCyl)
	for edge in t_edges[-1, 4:6, 1, 2]:
		edge.proj_geom(coldTubeCyl)
	for edge in t_edges[-1, 4, 1:3, 1]:
		edge.proj_geom(coldTubeCyl)

	xct = np.array([coldPlenumLen - 4, coldPlenumLen, petc])
	yct = np.array([cety - cetr, cety + cetr])
	zct = np.array([cetz - cetr, cetz + cetr])

	nxct = np.array([19, 40])
	nyct = np.full_like(yct, ND_t)
	nzct = np.full_like(zct, ndz[1])

	coldTube = bmh.CartBlockStruct(xct, yct, zct, nxct, nyct, nzct, zone='cold_tube')

	ct_faces = coldTube.faces
	ct_edges = coldTube.edges
	ct_bvts = coldTube.baked_vertices
	ct_grad = coldTube.grading

	for vt in ct_bvts[1:].flatten():
		vt.proj_geom(coldTubeCyl)
	for edge in ct_edges[1:, ..., 1:].flatten():
		if edge:
			edge.proj_geom(coldTubeCyl)
	for edge in ct_edges[..., 0].flatten():
		if edge:
			edge.proj_geom(coldTubeCyl)
	for face in ct_faces[..., 1:].flatten():
		if face:
			face.proj_geom(coldTubeCyl)

	ct_bvts[0] = t_bvts[-1, 4:6, 1:3]

	t_fmask[-1, :, :3, 0] = False

	ct_grad[0, ..., 0] = bmh.SimpleGradingElement(1 / 5)

# #######
# Grading
# #######

# Inlet
len_pcts = np.array([0.3, 0.4, 0.3])
dens = np.array([1.25, 1., 1., 1.25])

grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))

# Cold tube
t_divs[:ctInd, :, inInd - 2, 2] = 16
if not Axisym:
	c_divs[:, :, inInd - 2, 2] = 16

t_grad[:ctInd + 1, :, inInd - 2, 2] = bmh.SimpleGradingElement(2)
if not Axisym:
	c_grad[:, :, inInd - 2, 2] = bmh.SimpleGradingElement(2)

if IsPlenums:
	grd_elm = bmh.SimpleGradingElement(1 / 3)
else:
	n_elms = 30
	t_divs[:ctInd, :, inInd - 1, 2] = n_elms
	if not Axisym: c_divs[:ctInd, :, inInd - 1, 2] = n_elms

	len_pcts = np.array([0.2, 0.6, 0.2])
	dens = np.array([3.5, 1., 1., 2.5])
	grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))

t_grad[:wrInd, :, inInd - 1, 2] = grd_elm
if not Axisym: c_grad[:, :, inInd - 1, 2] = grd_elm

grd_elm = bmh.SimpleGradingElement(0.5)
t_grad[:ctInd + 1, :, inInd - 3, 2] = grd_elm
if not Axisym: c_grad[:, :, inInd - 3, 2] = grd_elm

# Cold Plenum
len_pcts = np.array([0.7, 0.3])
dens = np.array([1., 1., 2.5])
grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))
t_grad[ctInd + 1:, :, inInd - 2, 2] = grd_elm

# Main Tube
st_pct = 0.025
len_pcts = np.array([st_pct, 0.2 - st_pct, 0.6, 0.2])
dens = np.array([10., 4., 1., 1., 3.])
grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))
t_grad[..., inInd + 1, 2] = grd_elm
if not Axisym: c_grad[..., inInd + 1, 2] = grd_elm

# Working section inflation layers
len_pcts = np.array([0.8, 0.15, 0.05] if SmallYPlus else [0.8, 0.2])
dens = np.array([2., 2., 5., 120.] if SmallYPlus else [1., 1., 3.])
grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))
t_grad[wrInd - 1, :, inInd:csInd + 1, 0] = grd_elm

len_pcts = np.array([0.8, 0.15, 0.05] if SmallYPlus else [0.8, 0.2])
dens = np.array([2., 2., 5., 30.] if SmallYPlus else [1., 1., 3.])
grd_elm = bmh.MultiGradingElement(*bmh.get_grading_info(len_pcts, dens))
t_grad[wrInd - 1, :, csInd + 1, 0] = grd_elm

if ExtCone:
	t_grad[wrInd - 1, :, csInd, 0] = bmh.SimpleGradingElement(1 / 2.8)

# Hot Exit tube
if not Axisym: hotTube.grading[0, ..., 0] = bmh.SimpleGradingElement(8)

# Cold Exit tube
if not Axisym: coldTube.grading[0, ..., 0] = bmh.SimpleGradingElement(2)

# ###############
# Exclude Plenums
# ###############

if not IsPlenums:
	t_bmask[:, :, :inInd - 1] = True
	if not Axisym: c_bmask[:, :, :inInd - 1] = True

	t_bmask[:, :, csInd + 1:] = True
	if not Axisym: c_bmask[:, :, csInd + 1:] = True

# ##########
# Boundaries
# ##########

# Inlet
if IsInlets:
	inletInd = wrInd + 1
else:
	inletInd = wrInd

if Axisym or RibbonInlet:
	tube.boundary_tags[wrInd, :, inInd, 0].flatten()[:-1] = bmh.BoundaryTag('inlet')

else:
	tube.boundary_tags[inletInd, ::2, inInd, 0][:-1] = np.array([bmh.BoundaryTag(f'inlet{i}-inner') for i in range(4)])

n_arr = np.array(None)

if ExtCone:
	tube.boundary_tags[crInd, :, csInd, 2] = bmh.BoundaryTag(f'cone_start')
	tube.boundary_tags[crInd, :, csInd + 1, 2] = bmh.BoundaryTag(f'cone_end')

# Cold Exit
ce_bnd_tag = bmh.BoundaryTag('cold_outlet')

if IsPlenums:
	coldTube.boundary_tags[-1, ..., 0] = ce_bnd_tag

else:
	tube.boundary_tags[:, :, inInd - 1, 2] = ce_bnd_tag

	if not Axisym:
		core.boundary_tags[:, :, inInd - 1, 2] = ce_bnd_tag

# Hot Exit
he_bnd_tag = bmh.BoundaryTag('hot_outlet')

if IsPlenums:
	hotTube.boundary_tags[-1, ..., 0] = he_bnd_tag
else:
	tube.boundary_tags[:, :, csInd + 1, 2] = he_bnd_tag

# Testing mask
n = 6
# t_bmask[..., n:] = True
# c_bmask[..., n:] = True

# Periodic Faces
if Axisym:
	p1_tags = tube.boundary_tags[:-1, 0, :-1, 1]
	p1_tags[np.logical_not(t_bmask[:-1, 0, :-1])] = bmh.BoundaryTag('p1')

	p2_tags = tube.boundary_tags[:-1, 1, :-1, 1]
	p2_tags[np.logical_not(t_bmask[:-1, 0, :-1])] = bmh.BoundaryTag('p2')

if blocksOnly:
	tube.num_divisions[:] = 1

	if not Axisym:
		core.num_divisions[:] = 1

	if IsPlenums:
		hotTube.num_divisions[:] = 1
		coldTube.num_divisions[:] = 1

if Wedge:
	t_bmask[:, :(Nt - 1) // 2] = True

bmd = bmh.BlockMeshDict(metric='mm')

if IsPlenums:
	hotTube.write(bmd)
	coldTube.write(bmd)

if Axisym:
	tube.write(bmd)
else:
	vortex_tube.write(bmd)

bmd.write_file(run_blockMesh=False, run_renumberMesh=False, density_scale=DENSITY_SCALE, block_structure_only=blocksOnly)
