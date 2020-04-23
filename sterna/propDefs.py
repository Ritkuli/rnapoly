import bpy, math, bmesh
from bpy.props import *

#scale: 10^-9 m per blender unit
SIZE = 0.04
RADIUS = 0.87
TWIST = 32.73 / 180 * math.pi
AXIS = 139.9 / 180 * math.pi
INCLINATION = -0.745
RISE = 0.281
SCALE = 1.0
BASE_SCALE = 1.0
PADDING = ((RADIUS * TWIST) ** 2 + RISE ** 2) ** 0.5


def update_scale(self, context):
    if context.mode != "EDIT_MESH": return
    ob = context.object
    prev = ob.prevScaleParam
    cur = ob.scaleParam
    if ob.isSternaHelix:
        bm = bmesh.from_edit_mesh(ob.data)
        for v in bm.verts:
            v.select = 1
        bpy.ops.transform.resize(value=3*(prev/cur,))
    ob.prevScaleParam = cur



bpy.types.Object.rotationParam = FloatProperty(
    name="Rotation",
    min = -math.pi, max = math.pi,
    precision = 3,
    default = 0
)

bpy.types.Object.sizeParam = FloatProperty(
    name="Size",
    min = 0.001, max = 0.2,
    precision = 4,
    default = SIZE
)

bpy.types.Object.paddingParam = FloatProperty(
    name="Base Distance",
    min = 0.01, max = 10,
    precision = 4,
    default = PADDING
)

bpy.types.Object.minPadding = IntProperty(
    name="Min padding",
    default = 1,
    min = 0
)

bpy.types.Object.maxPadding = IntProperty(
    name="Max padding",
    default = 4,
    min = 0
)
"""-------------------------"""
bpy.types.Object.scaleParam = FloatProperty(
    name="Scale",
    min = 1, max = 100.0,
    precision = 4,
    step = 10,
    default = SCALE,
    update=update_scale
)
bpy.types.Object.prevScaleParam = FloatProperty(
    name="Previous scale",
    precision = 4,
    default = SCALE
)


bpy.types.Object.radiusParam = FloatProperty(
    name="Radius",
    min = 0.01, max = 1.0,
    precision = 4,
    default = RADIUS
)

bpy.types.Object.twistParam = FloatProperty(
    name="Twist",
    min = -3.14, max = 3.14,
    precision = 4,
    default = TWIST
)

bpy.types.Object.axisParam = FloatProperty(
    name="Axis",
    min = -3.14, max = 3.14,
    precision = 4,
    default = AXIS
)

bpy.types.Object.riseParam = FloatProperty(
    name="Rise",
    min = 0.01, max = 1.0,
    precision = 4,
    default = RISE
)

bpy.types.Object.inclinationParam = FloatProperty(
    name="Inclination",
    min = -1.0, max = 1.0,
    precision = 4,
    default = INCLINATION
)
"""-------------------------"""

bpy.types.Object.rstParam = BoolProperty(
    name="Use random spanning tree",
    default = False
)


bpy.types.Object.isSternaHelix = BoolProperty(
    name="isSternaHelix",
    default = False
)

bpy.types.Object.fivePrime = IntProperty(
    name="Five Prime",
    default = 0,
)


bpy.types.Object.genMesh = BoolProperty(
    name="Generate a mesh",
    default = False
)

bpy.types.Object.springRelaxOrder = FloatProperty(
    name="Spring order",
    min = 1.0, max = 100.0,
    default = 5
)


bpy.types.Object.setNickAtVertex = BoolProperty(
    name="Set nick at a vertex.",
    default = False
)


bpy.types.Object.ulParam = BoolProperty(
    name="Try to force integer number of turns",
    default = False
)


bpy.types.Object.springRelax = BoolProperty(
    name="Spring relax",
    default = False
)


bpy.types.Object.useAdaptiveOffset = BoolProperty(
    name="Use adaptive offset",
    default = False
)


bpy.types.Object.numTurnsParam = IntProperty(
    name="Number of full turns",
    default = 3,
    min = 0
)

bpy.types.Object.springRelaxSteps = IntProperty(
    name="Relaxation steps",
    min = 1, max = 10000,
    default = 100
)

bpy.types.Object.offsetParam = FloatProperty(
    name="Corner offset multiplier",
    min = -100.0, max = 100.0,
    precision = 4,
    default = 1
)

bpy.types.Object.klOffsetParam = IntProperty(
    name="KL offset",
    default = 0,
)

bpy.types.Object.klArmOffsetParam = IntProperty(
    name="KL Arm offset",
    default = 0,
)

bpy.types.Object.stemOffsetParam = IntProperty(
    name="Stem offset",
    default = 0,
)
