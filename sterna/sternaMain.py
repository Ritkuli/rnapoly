import bpy, bmesh, mathutils, sys, traceback
from . import blenderEdgeGroups, sternaStrand, propDefs
from .sternaStrand import get_bases, Strand, create_sterna_helix, create_mesh

class BASE_PRIMITIVE():

    verts = [
        mathutils.Vector((1, 0, 0)),
        mathutils.Vector((-1, 0, 0)),
        mathutils.Vector((0, 1, 0)),
        mathutils.Vector((0, -1, 0)),
        mathutils.Vector((0, 0, 1)),
        mathutils.Vector((0, 0, -1))
    ]

    faces = [
        (0, 2, 4),
        (0, 2, 5),
        (0, 3, 4),
        (0, 3, 5),
        (1, 2, 4),
        (1, 2, 5),
        (1, 3, 4),
        (1, 3, 5),
    ]

"""
def init():
    global blenderEdgeGroups
    try:
        blenderEdgeGroups
    except:
        blenderEdgeGroups = sys.modules[modulesNames['blenderEdgeGroups']]
    global sterna_data
    sterna_data = {}
"""

class AddSterna(bpy.types.Operator):
    """Add an empty Sterna object"""
    bl_idname = "mesh.add_sterna"
    bl_label = "Sterna Object"
    bl_options = {'REGISTER', 'UNDO'}

    #asdf = bpy.types.Object.numTurnsParam

    @classmethod
    def poll(cls, context):
        obj = context.object
        return context.mode == "OBJECT"


    def execute(self, context):
        mesh = bpy.data.meshes.new("test2")
        obj = bpy.data.objects.new("test", mesh)
        obj.isSternaHelix = True
        context.scene.objects.link(obj)
        obj.select = True
        bpy.context.scene.objects.active = obj
        return {'FINISHED'}

class AddSternaMenu(bpy.types.Menu):
    bl_label = "Sterna"
    bl_idname = "add_sterna_menu"

    def draw(self, context):
        self.layout.operator(AddSterna.bl_idname, icon='MESH_CUBE')
        self.layout.operator(AddSterna.bl_idname, icon='MESH_CUBE')

def draw_sterna_menu(self, context):
    layout = self.layout
    layout.operator(AddSterna.bl_idname, icon="PMARKER_SEL")



class SternaPanel(bpy.types.Panel):
    """A panel with the buttons and parameters required to generate an RNA nano-
    structure
    """
    bl_category = 'MyTab'
    bl_label = "Sterna"
    bl_idname = "sterna_main_panel"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "modifier"

    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj and obj.type in {'MESH', 'LATTICE'}) and not obj.isSternaHelix

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'

        obj = context.object
        me = context.mesh

        row = layout.row()
        row.label(text="Sterna", icon='WORLD_DATA')
        row = layout.row()
        row.label(text="Parameters")
        row = layout.row()
        row.label(text="Resolution")
        row = layout.row()
        row.prop(obj, 'scaleParam')
        row = layout.row()
        row.label(text="Visual parameters")
        row = layout.row()
        row.prop(obj, 'genMesh')
        if obj.genMesh:
            row.prop(obj, 'sizeParam')
        row = layout.row()
        row.label(text="Physical parameters")
        row = layout.row()
        row.prop(obj, 'radiusParam')
        row.prop(obj, 'twistParam')
        row = layout.row()
        row.prop(obj, 'axisParam')
        row.prop(obj, 'riseParam')
        row = layout.row()
        row.prop(obj, 'inclinationParam')
        row = layout.row()
        row.label(text="Generator settings")
        row = layout.row()
        row.prop(obj, 'offsetParam')
        row = layout.row()
        row.prop(obj, 'paddingParam')
        row = layout.row()
        row.prop(obj, 'minPadding')
        row.prop(obj, 'maxPadding')
        row = layout.row()
        row.prop(obj, 'klOffsetParam')
        row.prop(obj, 'stemOffsetParam')
        row = layout.row()
        row.prop(obj, 'klArmOffsetParam')
        row = layout.row()
        row.prop(obj, 'rstParam')
        if not obj.rstParam:
            row = layout.row()
            col = row.column()
            col.template_list("UI_UL_list", "edge_groups", obj, "edge_groups", obj, "active_edge_group_index")

            col = row.column(align=True)
            col.operator("object.edge_group_add", icon='ZOOMIN', text="")
            col.operator("object.edge_group_delete", icon='ZOOMOUT', text="")
            if obj.edge_groups:
                col.separator()
                col.operator("object.edge_group_move", icon='TRIA_UP', text="").direction = 'UP'
                col.operator("object.edge_group_move", icon='TRIA_DOWN', text="").direction = 'DOWN'

            if obj.edge_groups and obj.mode == 'EDIT':
                row = layout.row()

                sub = row.row(align=True)
                sub.operator("object.edge_group_assign", text="Assign")
                sub.operator("object.edge_group_remove", text="Remove")

                sub = row.row(align=True)
                sub.operator("object.edge_group_select", text="Select")
                sub.operator("object.edge_group_deselect", text="Deselect")

        row = layout.row()
        row.prop(obj, 'useAdaptiveOffset')
        row = layout.row()
        row.prop(obj, 'ulParam')
        row = layout.row()
        if obj.ulParam:
            row.prop(obj, 'numTurnsParam')
            row = layout.row()
        row = layout.row()
        row.label(text="Post processing")
        row = layout.row()
        row.prop(obj, 'setNickAtVertex')
        row = layout.row()
        row.prop(obj, 'springRelax')
        row = layout.row()
        row.prop(obj, 'springRelaxOrder')
        row.prop(obj, 'springRelaxSteps')
        row = layout.row()
        row.operator("object.generate_sterna")

class GenerateSterna(bpy.types.Operator):
    """ An operator generating a sterna helix from a mesh.
    """
    """Tooltip"""
    bl_idname = "object.generate_sterna"
    bl_label = "Generate Structure"
    bl_options = {'REGISTER', 'UNDO'}

    params = None

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.mode == "OBJECT" and obj.isSternaHelix != 1

    def execute(self, context):
        prev_mode = context.object.mode
        i_params = {p: getattr(context.object, p) for p in propDefs.edit_props}
        if self.params != i_params:
            self.params = i_params
            for p in propDefs.edit_props:
                setattr(self, p, i_params[p])
        mod_params = {p: getattr(self, p) for p in propDefs.edit_props}
        try:
            print("Mod Parameters:", mod_params)
            print("UI Parameters:", i_params)
            print(context.object)
            print(context.object.rstParam)
            bpy.ops.object.mode_set(mode="EDIT")
            bases, scale = get_bases(context, mode = sternaStrand.SternaGenModes.FULL, mod_params = mod_params)
        except Strand.StrandGenerationException as e:
            print(e)
            return {'FINISHED'}
        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("Unknown exception:", e.with_traceback(exc_traceback))
            traceback.print_exc()
            return {'FINISHED'}
        finally:
            bpy.ops.object.mode_set(mode=prev_mode)
        if context.object.genMesh: create_mesh(bases, context)
        ob = create_sterna_helix(context, bases, scale)
        #context.scene.objects.active = ob
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.object.mode_set(mode=prev_mode)
        return {'FINISHED'}

gs = {p: getattr(GenerateSterna, p) for p in ["execute", "poll", "bl_idname", "bl_label", "bl_options", "params"]}
gs.update(propDefs.edit_props)
GenerateSterna = type("gs", (bpy.types.Operator,), gs)


def find_start(bm):
    """ Finds the start vertex of the Sterna Helix by finding the shorterst distance to the end.

    Args:
        bm -- bmesh, Sterna object

    Returns:
        bmesh.Vertex
    """
    v1, i1 = traverse_to_end(bm, True)
    v2, i2 = traverse_to_end(bm, False)
    print(i1)
    print(i2)
    if i1 == i2:
        return bm.verts[0]
    elif i1 <= i2:
        return v1
    else:
        return v2

def traverse_to_end(bm, forwards = True):
    """ Traverses a Sterna Helix and returns the index and vertex of the final vertex.

    Args:
        bm -- bmesh, Sterna object

    KWArgs:
        forwards -- bool, direction of the traversal

    Returns:
        tuple<bmesh.Vertex, int>
    """
    visited = set()
    c_v = bm.verts[0]
    i = 0
    while True:
        i += 1
        v = c_v
        edges = v.link_edges if forwards else reversed(v.link_edges)
        for e in edges:
            e_verts = set([vx.index for vx in e.verts])
            index_other = sum(e_verts) - v.index
            other = bm.verts[index_other]
            if not e.seam and e.smooth:
                if not index_other in visited:
                    c_v = other
        visited.add(v.index)
        if v == c_v: break
    return c_v, i

def export_sterna_helix(object):
    """ Serializes the BMesh.

    Args:
        object -- bmesh, Sterna Object

    Returns:
        tuple<
            coordinates -- str,
            secondary_structure -- str,
            pseudoknot_numbering -- str,
            >

    """
    assert object.isSternaHelix

    bpy.ops.object.mode_set(mode="EDIT")
    bm = bmesh.from_edit_mesh(object.data)
    dotbracket = []
    pseudoknot_numbering = []
    coords = []

    visited = set()
    p_d = {}
    p_num = 0
    in_pknot = False
    prev_other = None
    bm.verts.ensure_lookup_table()

    c_v = bm.verts[object.fivePrime]
    while True:
        print(c_v)
        v = c_v
        coords.append(v.co * object.scaleParam)
        edges = v.link_edges
        current = "."
        for e in edges:
            e_verts = set([vx.index for vx in e.verts])
            index_other = sum(e_verts) - v.index
            other = bm.verts[index_other]
            if e.seam:
                if not index_other in visited:
                    current = "("
                else:
                    current = ")"
            elif not e.smooth:
                if not index_other in visited:
                    current = "["
                    if not in_pknot or bm.edges.get((other, prev_other)) == None:
                        p_num = len(set(pseudoknot_numbering))
                    prev_other = other # Need to check whether two pseudoknots start immediately consecutively
                else:
                    current = "]"
                    p_num = p_d[index_other]
            else:
                if not index_other in visited:
                    c_v = other

        dotbracket.append(current)
        p_d[v.index] = p_num
        visited.add(v.index)
        if current in "[]":
            if not in_pknot:
                pseudoknot_numbering.append(p_num)
            in_pknot = True
        else:
            in_pknot = False
        if c_v == v: break

    dotbracket = "".join(dotbracket)
    pseudoknot_numbering = " ".join([str(x) for x in pseudoknot_numbering])
    coords = ", ".join(["{} {} {}".format(*c) for c in coords])

    bpy.ops.object.mode_set(mode="OBJECT")

    return coords, dotbracket, pseudoknot_numbering


def import_sterna_helix(data):
    """ Converts a dictionary to an array of bases

    Args:
        data -- dict<str, str>, must contain secondary_structure, pseudoknot_numbering and positions

    Returns:
        list<Base> -- bases
    """
    secondary_structure = data["secondary_structure"]
    pseudoknot_numbering = data["pseudoknot_numbering"].split(" ")
    positions = data["positions"].strip("[] ").split(",")

    bases = []
    pairs = {}
    stack = []
    p_stack = {}
    p_num = None

    for i, sym in enumerate(secondary_structure):
        pos = positions[i].strip().split(" ")
        b = sternaStrand.Base(mathutils.Vector([float(x) for x in pos]))
        if sym in "[]":
            if p_num == None:
                p_num = pseudoknot_numbering.pop(0)
            b.p_num = p_num
            if p_num not in pseudoknot_numbering:
                p_stack[p_num][-1].pair = b
                b.pair = p_stack[p_num].pop()
            else:
                p_stack.setdefault(p_num, []).append(b)
        else:
            p_num = None
            if sym == "(":
                stack.append(b)
            elif sym == ")":
                stack[-1].pair = b
                b.pair = stack.pop()
        bases.append(b)
    return bases

def create_base_primitive(origin, size):
    """ Creates a base primitive at origin with a radius of size.

    Args:
        origin -- Vector
        size -- float
    """
    # Create mesh and object
    me = bpy.data.meshes.new('Base Mesh')
    ob = bpy.data.objects.new("Base", me)
    ob.location = origin

    # Link object to scene and make active
    scn = bpy.context.scene
    scn.objects.link(ob)
    scn.objects.active = ob
    ob.select = True

    verts = [size * v for v in BASE_PRIMITIVE.verts]

    # Create mesh from given verts, faces.
    me.from_pydata(verts, [], BASE_PRIMITIVE.faces)
    # Update mesh with new data
    me.update()
    """
    bpy.ops.mesh.primitive_uv_sphere_add(location=origin, size=size, segments=6, ring_count=6)
    """



def register():
    bpy.utils.register_class(GenerateSterna)
    bpy.utils.register_class(SternaPanel)
    bpy.utils.register_class(AddSterna)
    bpy.types.INFO_MT_mesh_add.prepend(draw_sterna_menu)


def unregister():
    bpy.utils.unregister_class(GenerateSterna)
    bpy.utils.unregister_class(SternaPanel)
    bpy.utils.unregister_class(AddSterna)
    bpy.types.INFO_MT_mesh_add.remove(draw_sterna_menu)


if __name__ == "__main__":
    register()
