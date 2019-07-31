import bpy, bmesh, random, mathutils, math, statistics
from enum import Enum
from operator import itemgetter
from . import blenderEdgeGroups

class SternaGenModes(Enum):
    FULL = 0
    SINGLE_HELIX = 1
    DOUBLE_HELIX = 2
    LOOP = 3
    KISSING_LOOP = 4


def get_bases(context, mode = SternaGenModes.FULL, mod_params = None):
    """ Returns a list of bases in a strand.

    Args:
        context -- Blender Context
        mode -- SternaGenModes enum
        mod_params -- dict<param, val>
    Returns:
        list<Base> -- bases
    """
    strand = Strand(context, mod_params = mod_params)
    #print("Mode:", mode)

    if mode == SternaGenModes.FULL:
         strand.full_traverse(context)
    else:
        t_mesh = bmesh.new()
        a = context.scene.cursor_location
        b = mathutils.Vector(mod_params["mousePositionParam"])
        if mode == SternaGenModes.LOOP: b *= 2 # Loops are half-length, compensate
        e = t_mesh.verts.new(a), t_mesh.verts.new(b)

        if mode == SternaGenModes.SINGLE_HELIX or mode == SternaGenModes.DOUBLE_HELIX:
            x = Strand.BeamNetwork.Beam(e, strand.params, "e")
            strand.get_edge(False, strand.get_orientation_vectors(x, False))
            if mode == SternaGenModes.DOUBLE_HELIX:
                strand.add_strand_separator()
                strand.get_edge(True, strand.get_orientation_vectors(x, True))
        elif mode == SternaGenModes.LOOP or mode == SternaGenModes.KISSING_LOOP:
            x = Strand.BeamNetwork.Beam(e, strand.params, "p")
            strand.get_pseudoknot(False, 0, strand.get_orientation_vectors(x, False))
            if mode == SternaGenModes.KISSING_LOOP:
                strand.add_strand_separator()
                strand.get_pseudoknot(True, 0, strand.get_orientation_vectors(x, True))

    bases = strand.get_bases()
    return bases


def create_mesh(bases, context):
    """ Creates a 3D mesh of a Sterna Helix

    Args:
        bases -- list<Base>
        context -- Blender Context
    """
    obj = context.object
    scale = obj.scaleParam;
    size = obj.sizeParam;

    prev_objects = set(bpy.context.scene.objects)
    for o in prev_objects:
        o.select = False

    for b in bases:
        create_base_primitive(b.co, scale * size)
    objects = set(bpy.context.scene.objects) - prev_objects
    for o in objects:
        o.select = True
        bpy.context.scene.objects.active = o
    bpy.ops.object.join()

def create_sterna_helix(context, bases, scale):
    """ Creates a Sterna Helix

    Args:
        context -- Blender Context
        bases -- list<Base>
        scale -- float, the relation between blender units and nanometers

    Returns:
        Blender.Object
    """
    me = bpy.data.meshes.new('Backbone')
    ob = bpy.data.objects.new("Sterna Helix", me)
    ob.location = context.object.location if context.object else mathutils.Vector()

    scn = bpy.context.scene
    scn.objects.link(ob)
    ob.isSternaHelix = True
    ob.scaleParam = scale
    ob.prevScaleParam = scale # also set the previous scale, since default=1
    prev_act = scn.objects.active
    scn.objects.active = ob
    ob.select = True

    bases_proper = []
    nicks = set()
    for i, b in enumerate(bases):
        if b:
            bases_proper.append(b)
        else:
            nicks.add(i - len(nicks))

    verts = [b.co for b in bases_proper]
    edges = []
    backbone = [(i, i + 1) for i in range(len(bases_proper) - 1) if i + 1 not in nicks]
    base_pairs = []
    pseudoknots = []
    visited = set()
    indices = []
    pseudo_indices = {}
    in_pknot = False
    for i, b in enumerate(bases_proper):
        if b.p_num != None:
            if b.p_num in pseudo_indices and not in_pknot:
                pseudoknots.append((i, pseudo_indices[b.p_num].pop()))
            else:
                t = pseudo_indices.setdefault(b.p_num, [])
                t.append(i)
                in_pknot = True
        else:
            in_pknot = False
            if b.pair:
                if b.pair in visited:
                    base_pairs.append((i, indices.pop()))
                else:
                    visited.add(b)
                    indices.append(i)
    edges.extend(backbone)
    edges.extend(base_pairs)
    edges.extend(pseudoknots)
    me.from_pydata(verts, edges, [])


    bpy.ops.object.mode_set(mode="EDIT")
    bm = bmesh.from_edit_mesh(me)
    # Mark base pairs
    for e in bm.edges:
        if e.index >= len(backbone):
            e.select = True
        else:
            e.select = False
    bpy.ops.mesh.mark_seam()
    # Mark Pseudoknots
    for e in bm.edges:
        if e.index >= len(backbone) + len(base_pairs):
            e.select = True
        else:
            e.select = False
    bpy.ops.mesh.mark_seam(clear = True)
    bpy.ops.mesh.mark_sharp()

    for v in bm.verts:
        v.select = True

    scn.objects.active = prev_act
    return ob


class Base():
    """ A base object containing the position, pseudoknot number and its base pair.
    """

    def __init__(self, origin, p_num = None, pair = None):
        self.co = origin
        self.p_num = p_num
        self.pair = pair


    def bond(self, other):
        """ Pairs this base with another
        """
        self.pair = other
        other.pair = self

    def __repr__(self):
        return(str(self.co))


class Strand():
    """ A representation of an RNA strand.
    """
    def __init__(self, context, mod_params = None):
        """
        Args:
            context -- Blender Context
            mod_params -- dict<param, val>, extra parameters not available through context
        """

        self.bases = []
        self.base_pairs = {}
        self.visited = set() # a set of visited edges
        self.p_visited = {} # A dictionary of visited pseudo knots mapping to their numbering
        self.cur = None # The coordinates of the last base inserted
        self.stack = [] # A stack of unpaired bases

        self.set_parameters(context, mod_params = mod_params)


    def set_parameters(self, context, mod_params = None):
        """ Setups the generation Parameters

        Args:
            context -- Blender Context
            mod_params -- dict<param, val>, extra parameters not available through context
        """
        obj = context.object
        true_scale = 1.0 / obj.scaleParam
        self.params = type("Parameters", (), {
            "scale": true_scale,
            "radiusParam": obj.radiusParam,
            "twistParam": obj.twistParam,
            "axisParam": obj.axisParam,
            "riseParam": obj.riseParam,
            "inclinationParam": obj.inclinationParam,
            "offsetParam": obj.offsetParam,
            "klOffsetParam": obj.klOffsetParam,
            "klArmOffsetParam": obj.klArmOffsetParam,
            "stemOffsetParam": obj.stemOffsetParam,
            "ulParam": obj.ulParam,
            "numTurnsParam": obj.numTurnsParam,
            "paddingParam": obj.paddingParam,
            "minPadding": obj.minPadding,
            "maxPadding": obj.maxPadding,
            "springRelax": obj.springRelax,
            "springRelaxSteps": obj.springRelaxSteps,
            "springRelaxOrder": obj.springRelaxOrder,
            "useAdaptiveOffset": obj.useAdaptiveOffset,
            "setNickAtVertex": obj.setNickAtVertex
        })
        if mod_params:
            for param in mod_params:
                setattr(self.params, param, mod_params[param])

    def setup_beam_network(self, context):
        """ Setups the beam network used in relaxing the helix orientations

        Args:
            context -- Blender Context
        """

        obj = context.object
        eo = context.edit_object
        rst = obj.rstParam;
        eg_index = eo.edge_groups[eo.active_edge_group_index].id if not rst else None
        mesh = bmesh.from_edit_mesh(eo.data)

        self.structure = get_structure(mesh, rst, eg_index)
        self.beam_network = self.BeamNetwork(self.structure, self.params)

        if self.params.springRelax:
            self.beam_network.relax(self.params.springRelaxSteps, self.params.springRelaxOrder)



    #TODO: Get rid of this
    def get_parameters(self):
        return self.params


    def get_orientation_vectors(self, beam, cur_visited):
        """ Calculates the start point, end point and the normal vector for the given edge.

        Args:
            beam -- Beam
            cur_visited -- bool, visited during traversal, i.e., 5'-3' or 3'-5'

        Returns:
            bmesh.Vector -- start point
            bmesh.Vector -- end point
            bmesh.Vector -- normal vector with a magnitude
            bmesh.Vector -- d_t; (end - start) = d_t * num
            int          -- number of nucleotides
        """
        a, b, n = beam.get_orientation(cur_visited)
        num = beam.get_length_bases()
        t = (b - a).normalized() * self.get_parameters().riseParam * self.get_parameters().scale

        return a, b, n, t, num

    def get_beam(self, e):
        """ Returns a beam from the beam network associated with an edge.

        Args:
            edge -- tuple<Vertex, Vertex>

        Returns:
            beam -- Beam
            beam visited -- bool
        """
        v1, v2 = e[1]
        cur_visited = (v2, v1) in self.visited
        beam = self.beam_network.get_beam(e)
        return beam, cur_visited


    def get_padding(self, orientation_vectors):
        """ Creates padding bases between two helices.
        Args:
            orientation_vectors -- tuple<Vector, Vector, Vector, Vector, int>

        Returns:
            num -- number of bases added
        """
        a, b, n = orientation_vectors[:3]
        end = self.cur
        start = a + n

        base_distance = self.get_parameters().paddingParam * self.get_parameters().scale
        num = math.ceil((start - end).length / base_distance)
        num = max(min(num, self.get_parameters().maxPadding), self.get_parameters().minPadding)
        t = 1.0 / (num + 1) * (start - end)
        for i in range(1, num + 1):
            coord = end + i * t
            self.bases.append(Base(coord))
            #self.cur = coord
        return num

    def get_pseudoknot(self, visited, p_num, orientation_vectors):
        """ Creates pseudoknot bases.

        Args:
            visited -- bool, edge has previously been visited
            p_num -- int, the number of the pseudoknot
            orientation_vectors -- tuple<Vector, Vector, Vector, Vector, int>

        Returns:
            bmesh.Vector -- coordinates of the 3'-nucleotide
        """
        a0, b0, n, t, num = orientation_vectors
        a = mathutils.Vector(a0)
        inclinationParam = self.get_parameters().inclinationParam * self.get_parameters().scale * t.normalized()
        coord = None
        if visited:
            t_n = mathutils.Vector(n)
            half = math.ceil(num / 2 - 4.5) + self.get_parameters().klArmOffsetParam #mathutils.Quaternion(t, self.get_parameter("twistParam"))
        else:
            t_n = mathutils.Vector(n)
            half = math.floor(num / 2 - 4.5) - self.get_parameters().klArmOffsetParam
        if num < 13:
            raise Strand.StrandGenerationException("Edge too short for a pseudoknot.")
        rot = mathutils.Quaternion(t, self.get_parameters().twistParam)

        # First Half
        for i in range(half):
            coord = a0 + t_n + i * t
            self.bases.append(Base(coord))
            t_n.rotate(rot)
            self.stack.append(self.bases[-1])

        ## switch direction
        a = a0 + inclinationParam + (half + 6) * t # half + 1 paddingParam + 6 pseudoknot
        t_n = mathutils.Vector(n)
        t_n.rotate(mathutils.Quaternion(t, self.get_parameters().axisParam + (half + 6) * self.get_parameters().twistParam)) # half + 1 paddingParam + 6 pseudoknot + axisParam
        t = -t

        # AA padding
        next_coord = a + t_n
        self.bases.append(Base((2/3) * coord + (1/3) * next_coord))
        self.bases.append(Base((1/3) * coord + (2/3) * next_coord))

        #Pseudo complex
        rot = mathutils.Quaternion(t, self.get_parameters().twistParam)
        for i in range(6):
            coord = a + t_n
            self.bases.append(Base(coord, p_num = p_num))
            t_n.rotate(rot)
            a += t

        # A padding
        coord = a + t_n
        self.bases.append(Base(coord))
        t_n.rotate(rot)

        # Second half
        for i in range(half):
            coord = a + (i + 1) * t + t_n
            self.bases.append(Base(coord))
            self.bases[-1].bond(self.stack.pop())
            t_n.rotate(rot)
        self.cur = coord
        return coord

    def get_edge(self, cur_visited, orientation_vectors):
        """ Creates the bases of an edge.

        Args:
            cur_visited -- bool, edge previously visited
            orientation_vectors -- tuple<Vector, Vector, Vector, Vector, int>

        Returns:
            bmesh.Vector -- coordinates of the 3'-nucleotide
        """
        a, b, n, t, num = orientation_vectors
        n_t = None
        coord = None
        for i in range(num):
            rot = mathutils.Quaternion(t, i * self.get_parameters().twistParam)
            n_t = mathutils.Vector(n)
            n_t.rotate(rot)
            coord = a + i * t + n_t
            self.bases.append(Base(coord))
            if cur_visited:
                self.bases[-1].bond(self.stack.pop())
            else:
                self.stack.append(self.bases[-1])
        self.cur = coord
        return coord

    def add_strand_separator(self):
        """ Marks the following nucleotides as belonging to a new strand.
        """
        self.bases.append(None)

    def full_traverse(self, context):
        """ Does the routing and returns the bases.

        Args:
            context -- Blender Context

        Returns:
            array<Base> -- list of bases
            scale -- float, scale factor used
        """
        self.setup_beam_network(context)

        #find longest edge for the potential nick
        longest = sorted(self.structure, key=lambda e: self.beam_network.get_beam(e).get_length_bases() if e[0] == "e" else -1)[-1]
        longest_pos = 0
        for e in self.structure:
            type = e[0]
            v1, v2 = e[1]
            #nick
            if e == longest: longest_pos = int(len(self.bases) + self.beam_network.get_beam(e).get_length_bases() / 2)
            beam, cur_visited = self.get_beam(e)
            orientation_vectors = list(self.get_orientation_vectors(beam, cur_visited))
            #Padding
            if self.cur:
                self.get_padding(orientation_vectors)
            #Edge
            if type == "e":
                visited = (v2, v1) in self.visited
                c = self.get_edge(visited, orientation_vectors)
                self.visited.add((v1, v2))
            #Pseudoknot
            elif type == "p":
                visited = (v2, v1) in self.p_visited
                if visited:
                    p_num = self.p_visited[(v2, v1)]
                else:
                    p_num = self.p_visited.setdefault((v1, v2), len(self.p_visited))
                c = self.get_pseudoknot(visited, p_num, orientation_vectors)
                self.visited.add((v1, v2))
                #print(len(self.visited))
        #self.normalize_scale()
        if not self.get_parameters().setNickAtVertex:
            self.place_nick_at(longest_pos)

    def get_bases(self):
        """ Returns the list of nucleotides generated with this objects

        Returns:
            list<Base> -- list of bases
        """
        return self.bases, 1.0 / self.get_parameters().scale

    def place_nick_at(self, idx):
        """ Sets the nick at index.

        Args:
            idx -- index of the base, where the nick is to be placed
        """
        #set the linkers at the first vertex
        e = self.structure[0]
        beam = self.get_beam(e)[0]
        orientation_vectors = list(self.get_orientation_vectors(beam, False))
        num = self.get_padding(orientation_vectors)
        #shift bases to the nick
        self.bases = self.bases[idx + num:] + self.bases[:idx + num]
        return

    def normalize_scale(self):
        """ Sets scale to 1 and moves all bases accordingly.
        """
        for b in self.bases:
            b.co = (self.scale / self.true_scale) * b.co
        self.scale = self.true_scale


    class StrandGenerationException(Exception):
        def __init__(self, message):
            super().__init__(message)

    class BeamNetwork():
        """ A beam network containing double helices as beams and methods
        to rescale and orientate the beams.
        """
        def __init__(self, structure, params):
            """
            Args:
                structure -- list<char, tuple<Vertex, Vertex>>
                params -- parameters object
            """
            if params.ulParam:
                for i in range(10):
                    params = self.scale_to_unit_length(structure, params)
            self.beams, self.mapping = self.create_beams(structure, params)
            self.params = params

        def create_beams(self, structure, params):
            """ Creates beams according to the structure.

            Args:
                structure -- list<char, tuple<Vertex, Vertex>>
                params -- parameters object

            Returns:
                list<Beam> -- beams array
                dict<tuple<Vertex, Vertex>, Beam> -- mapping of edges to the associated beam
            """
            created = {}
            mapping = {}
            last = None
            out_slot = 0
            for e in structure:
                if e[1] in created:
                    b = created[e[1]]
                    b.connect(last, 2, out_slot)
                    if e[0] == "e":
                        out_slot = 3
                    else:
                        out_slot = 1
                else:
                    b = self.Beam(e[1], params, type=e[0])
                    created[e[1][::-1]] = b
                    b.connect(last, 0, out_slot)
                    if e[0] == "e":
                        out_slot = 1
                    else:
                        out_slot = 3

                last = b
                mapping[e] = b
            mapping[structure[0]].connect(mapping[structure[-1]], 0, out_slot)
            beams = list(created.values())
            #print(len(beams), beams)
            return beams, mapping


        def relax(self, steps = 50, order = 1.0):
            """ Relaxes the orientations of the beam network
            using spring relaxation.

            Args:
                steps -- iterations, int
                order -- the order of the spring strength

            Returns:
                float -- total torque in the system
            """
            for i in range(steps):
                force_multiplier = 0.5 ** int(i / (steps / 5))
                tot = 0
                for b in self.beams[::(-1)**i]:
                    tot += b.relax(force_multiplier, order)
                #print("Total torque:", tot)
            return tot

        def get_parameters(self):
            """ Returns the parameters associated with this beam network.

            Returns:
                parameters object
            """
            return self.params

        def scale_to_unit_length(self, structure, params):
            """ Sets self.scale in such a way as to maximize the number of edges that are unit length.

            Args:
                structure -- list<char, tuple<Vertex, Vertex>>
                params -- parameters object

            Returns:
                object -- new parameters object
            """
            beams, mapping = self.create_beams(structure, params)

            tolerance = 10000#round((2 * math.pi / self.twistParam))
            mods = {}
            actual = {}
            for b in beams:
                turns = b.get_length() / params.riseParam / params.scale * params.twistParam / (2 * math.pi)
                #print("Turns: ", turns)
                f = round((turns / params.numTurnsParam) * tolerance)
                actual[f] = turns / params.numTurnsParam
                if mods.get(f):
                    mods[f] += 1
                else:
                    mods[f] = 1
            #print(mods)
            #print(actual)
            try:
                mod = statistics.mode(mods)
            except statistics.StatisticsError as e:
                mod = next(mods.__iter__())
            #print(mod)
            new_scale = params.scale * actual[mod]
            if new_scale == 0: new_scale = params.scale / 10
            proportion = new_scale / params.scale
            #print("Scale proportions:", proportion)
            #print("New scale:", new_scale)
            params.scale = new_scale
            #if proportion > 10 or proportion < 0.1:
            #    raise Exception("Could not find unit length scale.")
            return params


        def get_beam(self, edge):
            """ Returns the beam associated with the edge.

            Returns:
                Beam
            """
            return self.mapping[edge]


        class Beam():
            """ A beam object represents a double helix. It contains four input /
            output sockets that can connect to other beams.
            """
            def __init__(self, edge, params, type = "e"):
                """
                Args:
                    edge -- tuple<Vertex, Vertex>
                    params -- parameters object
                    type -- char, "e" for edge or "p" for pseudoknot
                """
                self.helical_axisParam = (edge[1].co - edge[0].co).normalized()
                axisParam = params.axisParam
                twistParam = params.twistParam
                scale = params.scale
                radiusParam = params.radiusParam
                half_inclinationParam = params.inclinationParam * scale * self.helical_axisParam * 0.5
                riseParam = params.riseParam
                start = edge[0].co + self.__get_offset(edge, edge[0], params) * self.helical_axisParam - half_inclinationParam
                end = edge[1].co - self.__get_offset(edge, edge[1], params) * self.helical_axisParam + half_inclinationParam
                length = (end - start).length
                if type == "e":
                    num_offset = params.stemOffsetParam
                elif type == "p":
                    num_offset = params.klOffsetParam
                self.num = int(length / scale / params.riseParam) + 1 + num_offset
                true_offset = - 0.5 * num_offset * riseParam * scale * self.helical_axisParam

                self.rotation = params.rotationParam if hasattr(params, "rotationParam") else 0
                rot = mathutils.Quaternion(self.helical_axisParam, math.pi / 2)
                self.north = self.helical_axisParam.orthogonal().normalized() * scale * radiusParam
                self.east = mathutils.Vector(self.north)
                self.east.rotate(rot)
                self.start = start + true_offset
                self.end = start + (self.num - 1) * riseParam * scale * self.helical_axisParam + true_offset
                self.edge = edge

                self.forces = []
                #Helix 1 A, Helix 1 B, Helix 2 A, Helix 2 B
                self.slots = [
                    0,
                    (self.num - 1) * twistParam,
                    (self.num - 1) * twistParam + axisParam,
                    axisParam
                ]
                self.positions = [
                    self.start - half_inclinationParam,
                    self.end - half_inclinationParam,
                    self.end + half_inclinationParam,
                    self.start + half_inclinationParam
                ]


            def is_sharp(self, edge, v0):
                """ Returns whether the edge comes in at a sharp angle to the vertex

                Args:
                    edge -- bmesh.Edge
                    v0 -- Vertex
                Returns:
                    bool
                """
                sign = None
                if len(v0.link_edges) < 2:
                    return True
                inc = 2 * v0.co - edge[0].co - edge[1].co
                for e in v0.link_edges:
                    if set(e.verts) == set(edge): continue
                    dir = 2 * v0.co - e.verts[0].co - e.verts[1].co
                    if inc.angle(dir) < math.pi / 4:
                        return True
                return False


            def __get_offset(self, edge, v0, params):
                """ Returns the offset required to avoid overlap between beams

                Args:
                    edge -- bmesh.Edge
                    v0 -- bmesh.Vertex
                    params -- parameters object

                Returns:
                    float -- offset in nanometers
                """
                #return self.get_cylinder_offset(edge, v0, params)
                if params.useAdaptiveOffset:
                    if self.is_sharp(edge, v0):
                        offset = self.__get_sphere_offset(v0, params)
                    else:
                        offset = self.__get_cylinder_offset(edge, v0, params)
                else:
                    offset = self.__get_sphere_offset(v0, params)
                return offset



            def __get_cylinder_offset(self, edge, v0, params):
                """ Calculates the offset required to avoid overlapping helices at vertices.

                Args:
                    edge -- bmesh.Edge, incoming edge
                    v0 -- bmesh.Vertex, end point
                    params -- parameters object

                Returns:
                    float -- offset along the edge
                """
                radiusParam = params.radiusParam * params.scale
                vectors = []
                m = float("inf")
                inc = 2 * v0.co - edge[0].co - edge[1].co
                if len(v0.link_edges) <= 1: return 0
                for e in v0.link_edges:
                    if set(e.verts) == set(edge): continue
                    other = 2 * v0.co - e.verts[0].co - e.verts[1].co
                    angle = inc.angle(other)
                    if angle < m:
                        m = angle
                h = radiusParam / math.sin(m / 2)
                return (h**2.0 - radiusParam**2.0)**0.5 * params.offsetParam


            def __get_sphere_offset(self, v0, params):
                """ Calculates the offset required to avoid overlapping helices at vertices.

                Args:
                    edge -- bmesh.Edge, incoming edge
                    v0 -- bmesh.Vertex, end point
                    params -- parameters object

                Returns:
                    float -- offset along the edge
                """
                radiusParam = params.radiusParam * params.scale
                vectors = []
                m = float("inf")
                for e in v0.link_edges:
                    #print(e)
                    if e.verts[0] == v0: vectors.append(e.verts[1].co - v0.co)
                    else: vectors.append(e.verts[0].co - v0.co)
                if len(vectors) <= 1: return 0
                for v1 in vectors:
                    for v2 in vectors:
                        if v2 == v1: continue
                        t = mathutils.Vector.angle(v1, v2)
                        #print(t)
                        if m > t: m = t
                h = radiusParam / math.sin(m / 2)
                return (h**2.0 - radiusParam**2.0)**0.5 * params.offsetParam



            def get_orientation(self, visited):
                """ Returns the orientation of one of the helices associated with this beam.
                Args:
                    visited -- boolean, true if visited before

                Returns:
                    Vector -- start point
                    Vector -- end point
                    Vector -- normal vector of the first base
                """
                #print(visited)
                if visited:
                    rot = mathutils.Quaternion(self.helical_axisParam, self.get_angle(2))
                    a = self.positions[2]
                    b = self.positions[3]
                else:
                    rot = mathutils.Quaternion(self.helical_axisParam, self.get_angle(0))
                    a = self.positions[0]
                    b = self.positions[1]
                n = mathutils.Vector(self.east)
                n.rotate(rot)
                return a, b, n

            def get_length(self):
                """ Returns the absolute length of this beam. Note:
                    this is longer than a single helix.

                Returns:
                    float
                """
                return (self.positions[1] - self.positions[0]).length

            def get_length_bases(self):
                """ Returns the length of one helix in the number of bases

                Returns:
                    int
                """
                return self.num


            def connect(self, other, slot, other_slot):
                """ Connects a socket of this beam to another socket in another
                beam and vice versa.

                Args:
                    other -- Beam, other beam
                    slot -- int, socket
                    other_slot -- int, other socket
                """
                if other == None: return
                self.__connect(other, slot, other_slot)
                other.__connect(self, other_slot, slot)


            def __connect(self, other, slot, other_slot):
                """ Connects a socket of this beam to another socket in another beam.

                Args:
                    other -- Beam, other beam
                    slot -- int, socket
                    other_slot -- int, other socket
                """
                force = lambda : self.get_force(other, slot, other_slot)
                self.forces.append(force)


            def get_position(self, slot):
                """ Returns the position of the base acting as a socket between
                two beams.

                Args:
                    slot -- int, socket

                Returns:
                    Vector -- position
                """
                angle = self.get_angle(slot)
                rot = mathutils.Quaternion(self.helical_axisParam, angle)
                n = mathutils.Vector(self.east)
                n.rotate(rot)
                pos = self.positions[slot] + n
                return pos


            def get_force(self, other, slot, other_slot):
                """ Returns the spring force experienced by one socket due to the
                beam it is connected to.

                Args:
                    other -- Beam, other beam
                    slot -- int, socket
                    other_slot -- int, other socket

                Returns:
                    float -- a signed angular force
                """
                a_rel = other.get_position(other_slot) - self.positions[slot]
                n = (a_rel - a_rel.dot(self.helical_axisParam) * self.helical_axisParam).normalized()
                angle = n.angle(self.east)

                #print(other.get_position(other_slot), self.get_position(slot), "\n--")
                #print(n)

                if n.dot(self.north) > 0: # bottom
                    angle = 2 * math.pi - angle

                dif = angle - self.get_angle(slot)
                if dif > math.pi: dif = dif - 2 * math.pi
                elif dif < -math.pi: dif = 2 * math.pi + dif
                #print("Difference: ", dif)
                return dif


            def get_angle(self, slot):
                """ Returns the angle of a socket in regards to the east
                direction of this beam.

                Args:
                    slot -- int, socket

                Returns:
                    float -- angle; [0, 2pi]
                """
                #print(slot)
                angle = self.slots[slot] + self.rotation
                angle -= int(angle / 2 / math.pi) * 2 * math.pi
                if angle < 0: angle += 2 * math.pi
                #print("Angle: ", angle)
                return angle

            def relax(self, force_multiplier = 1 / 10, order = 1.0):
                """ Relaxes the orientation of this beam according to the forces
                acting upon it by all other connected beams.

                Args:
                    force_multiplier -- float, defines the turn rate
                    order -- float, the spring force order

                Returns:
                    float -- total torque of this beam
                """
                #print(self.forces)
                tot = 0
                forces = []
                for force in self.forces:
                    #print(f())
                    f = force()
                    forces.append(f)
                    tot += abs(f) ** order
                for f in forces:
                    self.rotation += (f * abs(f) ** (order - 1) / tot + (random.random() - 0.5) * 1.52)  * force_multiplier
                #print("rotation:" , self.rotation)
                self.rotation -= int(self.rotation / 2 / math.pi) * 2 * math.pi
                if self.rotation < 0: self.rotation += math.pi * 2
                #print(":: ", tot)
                return abs(tot)



    def translate_to_world(self, structure, world):
        """ Translates all bases from local coordinates to world coordinates

        Args:
            structure -- list<char, tuple<Vertex, Vertex>>
            world -- the matrix of transformation
        """
        #TODO
        return structure
        #print(structure)
        for v in structure:
            v.co = world * v.co


def get_structure(mesh, rst, eg_index):
    """ Returns a routing of a mesh based on a spanning tree.

    Args:
        mesh -- Blender mesh
        rst -- bool, use random spanning tree rather than edge groups
        eg_index -- int, active edge group id, i.e., spanning tree

    Returns:
        list<char, tuple<Vertex, Vertex>> -- array of edges and their types as either stem "e"dges or "p"seudoknots
    """

    #dict = mesh_to_dict(mesh)
    spanning_tree, pseudo_knots = get_trees(mesh, rst, eg_index)
    start = spanning_tree[0].verts[0]
    structure = traverse(start, spanning_tree, pseudo_knots)
    #print("----")
    #print(structure)
    return structure


def get_trees(mesh, rnd, eg_index):
    """ Returns a spanning tree of the mesh.

    Args:
        mesh -- Blender mesh
        rnd -- bool, use random spanning tree rather than edge groups
        eg_index -- int, active edge group id, i.e., spanning tree

    Returns:
        list<Edge> -- edges in spanning tree
        list<Edge> -- edges in pseudonots
    """
    spanning_tree = []
    pseudo_knots = []
    if rnd:
        rst = get_rst(0, mesh)
    else:
        egl = mesh.edges.layers.string.get("edge_groups")
    for e in mesh.edges:
        if not rnd and blenderEdgeGroups.is_in_edge_group(e[egl], eg_index):
            spanning_tree.append(e)
        elif rnd and e in rst:
            spanning_tree.append(e)
        else:
            pseudo_knots.append(e)
    if not is_mst(mesh, spanning_tree):
        raise Exception("Invalid spanning tree")
    return spanning_tree, pseudo_knots


def is_mst(mesh, spanning_tree):
    """ Returns whether the given edges form a spanning tree of the mesh.

    Returns:
        bool
    """
    return True
    span = set()
    all = set()
    for e in mesh.edges:
        all.add(e.verts[0])
        all.add(e.verts[1])
    for e in spanning_tree:
        span.add(e.verts[0])
        span.add(e.verts[1])
    expected = len(mesh.verts) - 1
    return len(spanning_tree) == expected and span == all

def sort_edges(vertex, previous):
    """ Sorts edges leading out from a vertex according to the zig-zag selection algorithm.
    If the vertex defines a convex corner, uses its normal vector to define the poles of the
    zig-zag algorithm. Otherwise uses the incoming edge to define the poles.

    Args:
        vertex -- bmesh.Vertex
        previous -- bmesh.Edge, previous edge visited

    Returns:
        list<Edge>
    """
    #print("----")
    #print(vertex)
    #print(previous)
    #print("----")

    dir = vertex.normal.normalized()
    dir = 2 * len(vertex.link_edges) * vertex.co
    for x in vertex.link_edges:
        dir -= x.verts[0].co + x.verts[1].co
    if dir.length == 0: dir = 2 * vertex.co - previous.verts[0].co - previous.verts[1].co
    dir = dir.normalized()
    ort = dir.orthogonal().normalized()
    nor = dir.cross(ort).normalized()
    #print(dir, ort, nor)
    basis = mathutils.Matrix((dir, ort, nor)).transposed().inverted()
    edges = vertex.link_edges

    angles = []
    for e in edges:
        dir_t = e.verts[0].co + e.verts[1].co - 2 * vertex.co
        dir2 = basis * dir_t
        #print(":::", dir2.x * dir + dir2.y * ort + dir2.z * nor)
        phi = math.atan2(dir2.y, dir2.z)
        if phi < 0: phi += 2 * math.pi
        #print(phi, dir_t)
        theta = math.atan2(dir2.x, dir2.z)
        angles.append((e, phi, theta))
    sorted_edges = [x[0] for x in sorted(angles, key=itemgetter(1))]
    #print(sorted_edges)
    for i, e in enumerate(sorted_edges):
        if e == previous:
            sorted_edges = sorted_edges[i+1:] + sorted_edges[:i]
    return sorted_edges


def traverse(root, spanning_tree, pseudo_knots, prev=None):
    """Traverses twice around a graph based on its spanning tree.

    Args:
        spanning_tree -- list<Edge>
        pseduo_knots -- list<Edge>
        prev -- bmesh.Edge, previously visited edge

    Returns:
        structure
    """
    #print("Root:" , root)
    structure = []
    neighbors = []
    if prev == None:
        prev = root.link_edges[0]
        neighbors.append(prev)
    neighbors.extend(sort_edges(root, prev))
    # traverse
    for edge in neighbors:
        if edge.verts[0] == root:
            l = (root, edge.verts[1])
        else:
            l = (root, edge.verts[0])
        if edge in pseudo_knots:
            e = ("p", l)
            structure.append(e)
        elif edge in spanning_tree:
            e1 = ("e", l)
            e2 = ("e", (l[1], l[0]))
            sub = traverse(l[1], spanning_tree, pseudo_knots, prev=edge)
            structure.append(e1)
            structure.extend(sub)
            structure.append(e2)
        else:
            raise Exception("Disconnected edge")
    return structure


def get_rst(start, mesh, seed = 1):
    """Returns a random minimum spanning tree of the given graph.

    Args:
        start -- int, index of the start edge.
        mesh -- Blender Mesh
        seed -- int, random seed

    Returns
        set<Edge> -- A set of edges
    """
    mesh.edges.ensure_lookup_table()
    mesh.verts.ensure_lookup_table()
    rng = random.Random(seed)
    closed = set()
    edges = set()
    stack = [mesh.edges[start]]
    while stack:
        e = rng.choice(stack)
        v1 = e.verts[0]
        v2 = e.verts[1]
        stack.remove(e)
        if(v1 in closed and v2 in closed):
            continue
        closed.add(v1)
        closed.add(v2)
        edges.add(e)
        for e in set(v2.link_edges).union(v1.link_edges):
            stack.append(e)
    assert(len(edges) == len(mesh.verts) - 1)
    return edges



def mesh_to_dict(mesh):
    """Returns a dictionary mapping vertices to an ordered list of neighbors.

    Args:
        bmesh -- Blender mesh
    Returns:
        dict<Vertex, list<Vertex>>
    """
    dict = {}
    for f in mesh.faces:
        for i in range(len(f.verts)):
            l = dict.setdefault(f.verts[i], [])
            v0 = f.verts[i - 1]
            vc = f.verts[i]
            v1 = f.verts[(i + 1) % len(f.verts)]
            e0 = mesh.edges.get((vc, v0))
            e1 = mesh.edges.get((vc, v1))
            l.append(((v0, e0), (v1, e1)))
    for v in dict:
        #print(v, ":", dict[v], "\n")
        order = []
        _max = 100 # in case infinite loop. Should never happen
        while len(order) < len(v.link_edges) and _max > 0:
            _max -= 1
            assert(_max > 0)
            for t in dict[v]:
                if len(order) < 1:
                    order.append(t[0])
                else:
                    try:
                        if t[0] not in order:
                            order.insert(order.index(t[1]), t[0])
                        if t[1] not in order:
                            order.insert(order.index(t[0]) + 1, t[1])
                    except ValueError as e:
                        continue
        dict[v] = order
    return dict
