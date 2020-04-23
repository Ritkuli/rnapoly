import bpy, bmesh, sys, mathutils, traceback, math
from bpy.props import FloatVectorProperty, FloatProperty, IntProperty, BoolProperty
from . import sternaStrand
from .sternaStrand import get_bases, create_sterna_helix
from bpy_extras import view3d_utils

edit_props = ["mousePositionParam", "rotationParam", "klOffsetParam", "klArmOffsetParam"]


class SternaEditPanel(bpy.types.Panel):
    """A panel with the buttons and parameters required to generate an RNA nano-
    structure
    """
    bl_label = "Sterna"
    bl_idname = "STERNA_PT_edit_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Sterna'


    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj and obj.type in {'MESH', 'LATTICE'}) and obj.isSternaHelix

    def draw(self, context):
        layout = self.layout

        obj = context.object

        row = layout.row()
        row.label(text="Sterna", icon='WORLD_DATA')
        row = layout.row()
        if not validate_sterna(context):
            row.label(text="Invalid Sterna Helix")
            row = layout.row()
        row = layout.row()
        if obj.mode == "EDIT":
            row.prop(obj, 'scaleParam')
            row.label(text = "Bases: " + str(len(bmesh.from_edit_mesh(obj.data).verts)))
        else:
            row.label(text = "Scale: " + str(round(obj.scaleParam, 4)))
            row.label(text = "Bases: " + str(len(obj.data.vertices)))
        row = layout.row()
        row.operator("object.mark_five_prime")
        row.operator("object.select_five_prime")
        row = layout.row()
        row.operator("object.mark_base_pair")
        row.operator("object.select_base_pairs")
        row = layout.row()
        row.operator("object.mark_pseudoknot")
        row.operator("object.select_pseudoknots")
        row = layout.row()
        row.operator("mesh.hide")
        row.operator("mesh.reveal")
        row = layout.row()
        """ These are bugged in 2.81 and probably useless anyway
        row.label(text="Add")
        row = layout.row()
        row.operator("object.add_helix")
        row.operator("object.add_double_helix")
        row = layout.row()
        row.operator("object.add_loop")
        row.operator("object.add_kissing_loop")
        row = layout.row()
        """
        #row.label(text="Modify")
        #row = layout.row()
        row.operator("object.select_strand")


class SelectPseudoknots(bpy.types.Operator):
    """ An operator selecting pseudoknots.
    """
    """Tooltip"""
    bl_idname = "object.select_pseudoknots"
    bl_label = "Select Pseudoknots"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        bm = 1#blender_edge_groups.dic.get(obj.name)
        return obj.isSternaHelix == 1

    def execute(self, context):
        bpy.ops.object.mode_set(mode = 'OBJECT')
        for e in context.object.data.edges:
            if e.use_edge_sharp:
                e.select = True
        bpy.ops.object.mode_set(mode = 'EDIT')
        return {'FINISHED'}


class SelectBasePairs(bpy.types.Operator):
    """ An operator selecting base pairs.
    """
    """Tooltip"""
    bl_idname = "object.select_base_pairs"
    bl_label = "Select Base Pairs"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix == 1

    def execute(self, context):
        bpy.ops.object.mode_set(mode = 'OBJECT')
        for e in context.object.data.edges:
            if e.use_seam:
                e.select = True
        bpy.ops.object.mode_set(mode = 'EDIT')
        return {'FINISHED'}


class MarkFivePrime(bpy.types.Operator):
    """ An operator selecting base pairs.
    """
    """Tooltip"""
    bl_idname = "object.mark_five_prime"
    bl_label = "Mark Five Prime"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix and context.mode == "EDIT_MESH"

    def execute(self, context):
        bpy.ops.object.mode_set(mode = 'OBJECT')
        selection = []
        for v in context.object.data.vertices:
            if v.select:
                selection.append(v.index)
        if len(selection) == 1:
            context.object.fivePrime = selection[0]
        bpy.ops.object.mode_set(mode = 'EDIT')
        return {'FINISHED'}


class SelectFivePrime(bpy.types.Operator):
    """ An operator selecting base pairs.
    """
    """Tooltip"""
    bl_idname = "object.select_five_prime"
    bl_label = "Select Five Prime"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix == 1

    def execute(self, context):
        bpy.ops.object.mode_set(mode = 'EDIT')
        key = context.object.fivePrime
        ob = context.object
        bm = bmesh.from_edit_mesh(ob.data)
        bm.verts.ensure_lookup_table()
        v = bm.verts[key]
        v.select = True
        bm.select_history.add(v)
        # update view
        bmesh.update_edit_mesh(ob.data)
        return {'FINISHED'}


class MarkBasePair(bpy.types.Operator):
    """ An operator marking an edge as a base pair.
    """
    """Tooltip"""
    bl_idname = "object.mark_base_pair"
    bl_label = "Mark Base Pair"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix and context.mode == "EDIT_MESH"

    def execute(self, context):
        bpy.ops.mesh.mark_sharp(clear=True)
        bpy.ops.mesh.mark_seam()
        return {'FINISHED'}


class MarkPseudoknot(bpy.types.Operator):
    """ An operator marking an edge as a pseudoknot base pair.
    """
    """Tooltip"""
    bl_idname = "object.mark_pseudoknot"
    bl_label = "Mark Pseudoknot"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix and context.mode == "EDIT_MESH"

    def execute(self, context):
        bpy.ops.mesh.mark_seam(clear=True)
        bpy.ops.mesh.mark_sharp()
        return {'FINISHED'}


class SelectStrand(bpy.types.Operator):
    """ An operator for selecting the vertices of a strand
    """
    """Tooltip"""
    bl_idname = "object.select_strand"
    bl_label = "Select Strand"
    bl_options = {'REGISTER', 'UNDO'}

    selectionCount: IntProperty(
        name="Select #",
    )

    selectBothStrands: BoolProperty(
        name="Also select paired strand.",
    )

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.isSternaHelix and context.mode == "EDIT_MESH"


    def execute(self, context):
        obj = context.object
        bm = bmesh.from_edit_mesh(obj.data)
        start_points = []
        for v in bm.verts:
            if v.select:
                start_points.append(v)
        visited_edges = set()
        visited_vertices = set(start_points)
        bm.select_history.clear()
        for v in start_points:
            c_v = v
            for i in range(abs(self.selectionCount) + 1):
                c_v.select = True
                bm.select_history.add(c_v)
                if self.selectBothStrands:
                    pair_edges = [e for e in c_v.link_edges if e.seam]
                    if len(pair_edges) <1:
                        pass
                    else:
                        v1, v2 = pair_edges[0].verts
                        if c_v == v1:
                            v2.select = True
                        else:
                            v1.select = True
                edges = [e for e in c_v.link_edges if e.smooth and not e.seam and e not in visited_edges]
                if len(edges) < 1:
                    break
                elif self.selectionCount > 0:
                    edge = edges[0]
                else:
                    edge = edges[-1]
                visited_edges.add(edge)
                v1,v2 = edge.verts
                if v1 not in visited_vertices:
                    c_v = v1
                elif v2 not in visited_vertices:
                    c_v = v2
                else:
                    break
                visited_vertices.add(c_v)
        bmesh.update_edit_mesh(obj.data)
        return {'FINISHED'}


class AddHelix(bpy.types.Operator):
    """ An operator generating a sterna helix from a mesh.
    """
    """Tooltip"""
    bl_idname = "object.add_helix"
    bl_label = "Add Helix"
    bl_options = {'REGISTER', 'UNDO'}

    mousePositionParam: FloatVectorProperty(
        name="mousePositionParam",
        size=3,
    )

    rotationParam: bpy.types.Object.rotationParam

    @classmethod
    def poll(cls, context):
        obj = context.object
        return context.mode == "EDIT_MESH" and obj.isSternaHelix == 1

    def modal(self, context, event):
        return sterna_edit_modal(self, context, event)

    def invoke(self, context, event):
        return sterna_edit_invoke(self, context, event)

    def execute(self, context):
        return sterna_edit_execute(self, context, sternaStrand.SternaGenModes.SINGLE_HELIX)


class AddDoubleHelix(bpy.types.Operator):
    """ An operator generating a sterna helix from a mesh.
    """
    """Tooltip"""
    bl_idname = "object.add_double_helix"
    bl_label = "Add Double Helix"
    bl_options = {'REGISTER', 'UNDO'}

    mousePositionParam: FloatVectorProperty(
        name="mousePositionParam",
        size=3,
    )

    rotationParam: bpy.types.Object.rotationParam

    @classmethod
    def poll(cls, context):
        obj = context.object
        return context.mode == "EDIT_MESH" and obj.isSternaHelix == 1

    def modal(self, context, event):
        return sterna_edit_modal(self, context, event)

    def invoke(self, context, event):
        return sterna_edit_invoke(self, context, event)

    def execute(self, context):
        return sterna_edit_execute(self, context, sternaStrand.SternaGenModes.DOUBLE_HELIX)


class AddKissingLoop(bpy.types.Operator):
    """ An operator generating a sterna helix from a mesh.
    """
    """Tooltip"""
    bl_idname = "object.add_kissing_loop"
    bl_label = "Add Kissing Loop"
    bl_options = {'REGISTER', 'UNDO'}

    mousePositionParam: FloatVectorProperty(
        name="mousePositionParam",
        size=3,
    )

    rotationParam: bpy.types.Object.rotationParam
    #klOffsetParam = bpy.types.Object.klOffsetParam
    klArmOffsetParam: bpy.types.Object.klArmOffsetParam

    @classmethod
    def poll(cls, context):
        obj = context.object
        return context.mode == "EDIT_MESH" and obj.isSternaHelix == 1

    def modal(self, context, event):
        return sterna_edit_modal(self, context, event)

    def invoke(self, context, event):
        return sterna_edit_invoke(self, context, event)

    def execute(self, context):
        return sterna_edit_execute(self, context, sternaStrand.SternaGenModes.KISSING_LOOP)

"""
from . import propDefs
aks = {p: getattr(AddKissingLoop, p) for p in ["execute", "poll", "modal", "invoke", "bl_idname", "bl_label", "bl_options", "offset"]}
aks.update(propDefs.edit_props)
AddKissingLoop = type("aks", (bpy.types.Operator,), aks)
"""


class AddLoop(bpy.types.Operator):
    """ An operator generating a sterna helix from a mesh.
    """
    """Tooltip"""
    bl_idname = "object.add_loop"
    bl_label = "Add Loop"
    bl_options = {'REGISTER', 'UNDO'}

    mousePositionParam: FloatVectorProperty(
        name="mousePositionParam",
        size=3,
    )

    rotationParam: bpy.types.Object.rotationParam

    @classmethod
    def poll(cls, context):
        obj = context.object
        return context.mode == "EDIT_MESH" and obj.isSternaHelix == 1

    def modal(self, context, event):
        return sterna_edit_modal(self, context, event)

    def invoke(self, context, event):
        return sterna_edit_invoke(self, context, event)

    def execute(self, context):
        return sterna_edit_execute(self, context, sternaStrand.SternaGenModes.LOOP)


def sterna_edit_modal(self, context, event):
    #print(event.type)
    if event.type in {'MOUSEMOVE', 'INBETWEEN_MOUSEMOVE'}:
        cursor = context.scene.cursor_location
        diff = mathutils.Vector((event.mouse_x - event.mouse_prev_x, event.mouse_y - event.mouse_prev_y, 0))
        accuracy = 2 if self.shift else 1
        self.mouse_offset += diff * 10 ** (-accuracy + 1)
        mouse_3d = view3d_utils.region_2d_to_location_3d(self.region, self.region_3d, self.first_mouse + self.mouse_offset, cursor)
        if self.ctrl:
            mouse_3d = mathutils.Vector((round(c, int(math.log(context.object.scaleParam, 10) + accuracy - 1)) for c in mouse_3d))
        if self.lock_x:
            mouse_3d.y = cursor.y
            mouse_3d.z = cursor.z
        if self.lock_y:
            mouse_3d.x = cursor.x
            mouse_3d.z = cursor.z
        if self.lock_z:
            mouse_3d.x = cursor.x
            mouse_3d.y = cursor.y
        self.mousePositionParam = mouse_3d
        bpy.ops.mesh.delete()
        self.execute(context)
        context.area.header_text_set("Offset %.4f %.4f %.4f" % tuple(self.mousePositionParam))

    elif event.type == 'LEFTMOUSE':
        context.area.header_text_set()
        return {'FINISHED'}

    elif event.type in {'LEFT_ALT', 'RIGHT_ALT'}:
        if event.value == 'PRESS':
            self.alt = True
        if event.value == 'RELEASE':
            self.alt = False
        return {'RUNNING_MODAL'}

    elif event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}:
        if event.value == 'PRESS':
            self.ctrl = True
        if event.value == 'RELEASE':
            self.ctrl = False
        return {'RUNNING_MODAL'}

    elif event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}:
        if event.value == 'PRESS':
            self.shift = True
        if event.value == 'RELEASE':
            self.shift = False

    elif event.type in {'RIGHTMOUSE', 'ESC'}:
        context.area.header_text_set()
        return {'CANCELLED'}

    elif event.type in 'xX':
        if event.value == 'PRESS':
            self.lock_x = not self.lock_x
            self.lock_y = self.lock_z = False

    elif event.type in 'yY':
        if event.value == 'PRESS':
            self.lock_y = not self.lock_y
            self.lock_x = self.lock_z = False

    elif event.type in 'zZ':
        if event.value == 'PRESS':
            self.lock_z = not self.lock_z
            self.lock_x = self.lock_y = False

    else:
        self.first_mouse = mathutils.Vector((event.mouse_x, event.mouse_y, 0)) - self.region_offset
        self.mouse_offset = mathutils.Vector()
        return {'PASS_THROUGH'}

    return {'RUNNING_MODAL'}


def sterna_edit_invoke(self, context, event):
    for area in context.screen.areas:
        if area.type == "VIEW_3D":
            self.view_distance = area.spaces[0].region_3d.view_distance
            self.region_3d = area.spaces[0].region_3d
            for region in area.regions:
                if region.type == "WINDOW":
                    self.region = region
                    self.region_offset = mathutils.Vector((region.x, region.y, 0))
            break
    self._initial_mouse = mathutils.Vector((event.mouse_x, event.mouse_y, 0.0))
    context.window_manager.modal_handler_add(self)
    self.selected = set()
    for v in bmesh.from_edit_mesh(context.object.data).verts:
        if v.select:
            self.selected.add(v)
        v.select = False
    bpy.context.scene.objects.active = bpy.context.scene.objects.active

    self.shift = self.ctrl = self.alt = False
    self.lock_x = self.lock_y = self.lock_z = False
    self.first_mouse = mathutils.Vector((event.mouse_x, event.mouse_y, 0)) - self.region_offset
    self.mouse_offset = mathutils.Vector()
    return {'RUNNING_MODAL'}


def sterna_edit_execute(self, context, cmd):
    mod_params = {p: getattr(self, p) for p in edit_props if hasattr(self, p)}
    #mod_params.update({p: getattr(self, p) for p in propDefs.edit_props})
    try:
        bases, scale = get_bases(context, mode = cmd, mod_params=mod_params)
        bpy.ops.object.mode_set(mode="OBJECT")
        active_object = context.active_object
        new_object = create_sterna_helix(context, bases, scale)
        join_objects(context, active_object, new_object)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print("Unknown exception:", e.with_traceback(exc_traceback))
        traceback.print_exc()
        bpy.ops.object.mode_set(mode="EDIT")
        return {'CANCELLED'}
    return {'FINISHED'}


def validate_sterna(context):
    """ Determines the validity of a sterna helix.

    Args:
        context -- Blender Context

    Returns:
        bool
    """
    if len(context.object.data.polygons) > 0:
        return False

    return True


def join_objects(context, active_object, new_object):
    """ Joins Blender objects into one object

    Args:
        context -- Blender Context
        active_object -- Blender Object
        new_object -- Blender Object
    """
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.object.mode_set(mode="OBJECT")
    context.scene.objects.active = active_object
    new_object.select = True
    bpy.ops.object.join()
    bpy.ops.object.mode_set(mode="EDIT")


def register():
    bpy.utils.register_class(AddHelix)
    bpy.utils.register_class(AddDoubleHelix)
    bpy.utils.register_class(AddKissingLoop)
    bpy.utils.register_class(AddLoop)
    bpy.utils.register_class(SelectPseudoknots)
    bpy.utils.register_class(SelectBasePairs)
    bpy.utils.register_class(SelectFivePrime)
    bpy.utils.register_class(MarkFivePrime)
    bpy.utils.register_class(MarkBasePair)
    bpy.utils.register_class(MarkPseudoknot)
    bpy.utils.register_class(SelectStrand)
    bpy.utils.register_class(SternaEditPanel)


def unregister():
    bpy.utils.unregister_class(AddHelix)
    bpy.utils.unregister_class(AddDoubleHelix)
    bpy.utils.unregister_class(AddKissingLoop)
    bpy.utils.unregister_class(AddLoop)
    bpy.utils.unregister_class(SelectPseudoknots)
    bpy.utils.unregister_class(SelectBasePairs)
    bpy.utils.unregister_class(MarkFivePrime)
    bpy.utils.unregister_class(SelectFivePrime)
    bpy.utils.unregister_class(MarkBasePair)
    bpy.utils.unregister_class(MarkPseudoknot)
    bpy.utils.unregister_class(SelectStrand)
    bpy.utils.unregister_class(SternaEditPanel)
