import mbuild as mb
import numpy as np
import parmed as pmd
from foyer import Forcefield
import MDAnalysis as mda

def get_player_info():
    print("I am the main player")

class builder():
    """ builder tool to build graphene, cdc, and graphene + cdc.
        Made reference to Ray's Porebuilder package
    Parameters
    ----------
    pore_length : int, default=4
        dimensions of graphene sheet length in nm
    pore_depth : int, default=4
        dimensions of graphene sheet depth in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    pore_width: int, default=1
        width of slit pore in nm
    slit_pore_dim : int, default=1
        dimension slit pore, default is in the y-axis
    Attributes
    ----------
    see mbuild.Compound
    """
    # def __init__(self):
    #     # super(GraphenePore, self).__init__()
    #     # self.pore_length = pore_length
    #     # self.pore_depth = pore_depth
    #     # self.n_sheets = n_sheets
    #     print('working')

    def make_graphene(self, ff_file_path, pore_length, pore_depth, n_sheets):
        self.pore_length = pore_length
        self.pore_depth = pore_depth
        self.n_sheets = n_sheets
        factor = np.cos(np.pi/6)
        # Estimate the number of lattice repeat units
        replicate = [int(self.pore_length/0.2456), (self.pore_depth/0.2456)*(1/factor)]
        if all(x <= 0 for x in [self.pore_length, self.pore_depth]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        carbon = mb.Compound()
        carbon.name = 'C'
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {carbon.name: carbon_locations}
        lattice_spacing = [0.2456, 0.2456, 0.335]
        angles = [90.0, 90.0, 120.0]
        graphene_lattice = mb.Lattice(lattice_spacing=lattice_spacing,
                                      angles=angles, lattice_points=basis)

        graphene = graphene_lattice.populate(compound_dict={carbon.name: carbon},
                                             x=replicate[0], y=replicate[1],
                                             z=self.n_sheets)
        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.box.Lx
        boundary = [graphene.box.Lx*10, graphene.box.Ly*10 * factor, graphene.box.Lz*10, 90, 90, 90]
        graphene = graphene.to_parmed()
        graphene.box = boundary 
        gra_ff = Forcefield(forcefield_files=ff_file_path)
        graphene = gra_ff.apply(graphene, residues = 'gra')
        for res in graphene.residues:
            res.name = 'gra'
        # graphene.box = boundary
        self._graphene = graphene
        return self._graphene
        

    def cdc(self, box_lengths, ff_file_path, trj_file_path, n_atoms):
        """
        Parameters
        --------
        box_lengths: list, float
            box lengths in x, y, z direction in nm
        ff_file_path:
            path to force field file
        trj_file_path:
            path to trjectory file like xtc file
        n_atoms: int
            the number of carbon atoms 
        """
        carbon = mb.Compound()
        carbon.name = "C"
        box = mb.Box(lengths = box_lengths)
        box_system = mb.fill_box(carbon, n_atoms, box = box)
        cdc_ff = Forcefield(forcefield_files=ff_file_path)
        cdc = cdc_ff.apply(box_system, residues = 'cdc')
        for res in cdc.residues:
            res.name = 'cdc'
        universe = mda.Universe(trj_file_path)
        positions = universe.trajectory[0].positions
        i = 0 
        for atom in cdc.atoms:
            atom.xx = positions[i][0]
            atom.xy = positions[i][1]
            atom.xz = positions[i][2]
            i += 1
        self._cdc = cdc
        return self._cdc

    def carbon_slit_pore(self, ff_file_path, pore_size, pore_length, pore_depth, n_sheets):
        """
        Parameters
        --------
        ff_file_path: string
            path to force field
        pore_size: float, int, nm
            the size of pore between graphene block
        pore_length: float, nm
            the length of the pore
        pore_depth: float, nm
            the depth of the pore
        n_sheets: int
            the number of sheets in each grahene block
        """
        pore_length = pore_length/2
        self.make_graphene(ff_file_path, pore_length, pore_depth, n_sheets)
        graphene = self.rotate(axis = 2, target='graphene')
        graphene_copy = graphene.__copy__()
        graphene_vertial_negative = self.rotate(axis = 0 , target = 'graphene')
        
        ### set up line equation for plane
        z_values = list(set(graphene_copy.coordinates[:,2]))
        y_values = list(set(graphene_vertial_negative.coordinates[:,1]))
        coefficients = np.polyfit(y_values, z_values, 1)
        plane = [2, '>', 1, 1, -0.5]
        graphene_copy = self.remove_partial(graphene_copy, plane)
        plane = [2, '<', 1, 1, 0.0]
        graphene_vertial_negative = self.remove_partial(graphene_vertial_negative, plane)
        graphene_right = graphene_vertial_negative + graphene_copy

        ### make symmetric graphene_right to get left graphene
        distance = np.max(graphene_right.coordinates[:, 1])
        graphene_left = self.symmetric(graphene_right, 1, distance, edge = 0.3)
        graphene_slit = graphene_left + graphene_right
        factor = np.cos(np.pi/6) * 2.456
        # factor = 2.456
        number = 7
        graphene_slit = self.remove_part(graphene_slit, 2, ">", factor * number)
        graphene_slit_below = self.symmetric(graphene_slit, 2, 0, -1000)
        graphene_slit_below = self.translate(graphene_slit_below, 2, factor * number)
        pore_size = pore_size * 10 # angstrom
        graphene_slit = self.translate(graphene_slit, 2, factor * number + pore_size)

        graphene_slit_pore = graphene_slit + graphene_slit_below
        graphene_slit_pore.box[1] = np.max(graphene_slit_pore.coordinates[:,1]) * 2
        graphene_slit_pore.box[2] = 2 * factor * number + pore_size 
        max_limit = np.max(graphene_slit_pore.coordinates[:, 2])
        graphene_slit_pore = self.remove_part(graphene_slit_pore, 2, ">", max_limit-0.1)
        return graphene_slit_pore



    def rotate(self, axis, target):
        """  Rotate is for graphene for now
        Parameters
        --------
        structure: parmed.structure
        axis: int
            0: x axis
            1: y axis
            2: z axis
        """
        if target == 'graphene':
            if axis == 0:
                for atom in self._graphene.atoms:
                    xy = atom.xy
                    xz = atom.xz
                    atom.xy = xz
                    atom.xz = xy
                box = self._graphene.box
                self._graphene.box = [box[0], box[2], box[1], 90, 90, 90]
            if axis == 2:
                for atom in self._graphene.atoms:
                    xx = atom.xx
                    xy = atom.xy
                    atom.xx = xy
                    atom.xy = xx
                box = self._graphene.box
                self._graphene.box = [box[1], box[0], box[2], 90, 90, 90]
                # box = self._graphene.box
            return self._graphene


    def remove_partial_cdc(self, axis, limit):
        """ remove the part of cdc which is below limit in axis direction

        Parameters
        --------
        structure: parmed.structure
        axis: int
            0: x axis
            1: y axis
            2: z axis
        limit: float
            below this limit, atoms in this structure will be removed
        """
        count = np.where(self._cdc.coordinates[:, axis] < limit)
        n = len(count[0])
        if axis == 1:
            for i in range(n):
                for atom in self._cdc.atoms:
                    if atom.xy < limit:
                        self._cdc.atoms.remove(atom)
        if axis == 0:
            for i in range(n):
                for atom in self._cdc.atoms:
                    if atom.xx < limit:
                        self._cdc.atoms.remove(atom)
        if axis == 2:
            for i in range(n):
                for atom in self._cdc.atoms:
                    if atom.xz < limit:
                        self._cdc.atoms.remove(atom)
        return self._cdc

    def remove_part(self, structure, axis, dir, limit):
        """ remove the part of cdc which is below limit in axis direction

        Parameters
        --------
        structure: parmed.structure
        axis: int
            0: yz plane
            1: xz plane
            2: xy plane
        dir: "<" or ">"
            "<": remove below limit
            ">": remove above limit
        limit: float, angstrom
            below this limit, atoms in this structure will be removed
        """

        count = np.where(structure.coordinates[:, axis] < limit)
        n = len(count[0])
        if axis == 1 and dir == "<":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xy < limit:
                        structure.atoms.remove(atom)
        if axis == 1 and dir == ">":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xy > limit:
                        structure.atoms.remove(atom)
        if axis == 0 and dir == "<":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xx < limit:
                        structure.atoms.remove(atom)
        if axis == 0 and dir == ">":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xx > limit:
                        structure.atoms.remove(atom)
        if axis == 2 and dir == "<":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xz < limit:
                        structure.atoms.remove(atom)
        if axis == 2 and dir == ">":
            for i in range(n):
                for atom in structure.atoms:
                    if atom.xz > limit:
                        structure.atoms.remove(atom)
        return structure

    def remove_partial(self, structure, plane):
        """ remove the part of structure which is below or above the plane
            currently, only support z = f(y)

        Parameters
        --------
        structure: parmed.structure
        limit: float
            below this limit, atoms in this structure will be removed
        plane: array_like， For example, z = a*y + b, plane = [2, '>', a, 1, b]
            plane[0] : dependent variable 
                0 : x
                1 : y
                2 : z
            plane[1] : operator
            plane[2] : coefficient a
            plane[3] : independent variable
            plane[4] : coefficient b
        """
        if plane[1] == ">":
            count = np.where(structure.coordinates[:, plane[0]] > (plane[2]* structure.coordinates[:, plane[3]] + plane[4]))
            n = len(count[0])
            if plane[0] == 2 and plane[3] == 1:
                for i in range(n):
                    for atom in structure.atoms:
                        if atom.xz > (plane[2]*atom.xy + plane[4]):
                            # print(structure)
                            structure.atoms.remove(atom)
        if plane[1] == "<":
            count = np.where(structure.coordinates[:, plane[0]] < (plane[2]* structure.coordinates[:, plane[3]] + plane[4]))
            n = len(count[0])
            if plane[0] == 2 and plane[3] == 1:
                for i in range(n):
                    for atom in structure.atoms:
                        if atom.xz < (plane[2]*atom.xy + plane[4]):
                            structure.atoms.remove(atom)
        return structure

    def symmetric(self, structure, plane, distance, edge):
        """
        Parameters
        --------
        structure: parmed.structure
        plane: int
            0: yz plane
            1: xz plane
            2: xy plane
        distance: float
            sign like +,- means direction; the distance to plane
        edge: float
            within the edge, the atoms for sysmetry are not considered
        """
        structure_copy = structure.__copy__()
        if plane == 2:
            for atom in structure_copy.atoms:
                distance_to_middle_line = distance - atom.xz
                if distance_to_middle_line > edge:
                    atom.xz += 2 * distance_to_middle_line
        if plane == 1:
            for atom in structure_copy.atoms:
                distance_to_middle_line = distance - atom.xy
                if distance_to_middle_line > edge:
                    atom.xy += 2 * distance_to_middle_line
        return structure_copy

    def move_structure(self, axis, distance, target):
        """
        Parameters
        --------
        structure: parmed.structure
        axis: int
            0: x axis
            1: y axis
            2: z axis
        distance: float
            sign like +,- means direction
        target: string
            'graphene'
            'cdc'
            'all'
            the objective that you want to move
        """
        if target == 'graphene':
            structure = self._graphene
        if target == 'cdc':
            structure = self._cdc
        
        if axis == 0:
            for atom in structure.atoms:
                atom.xx += distance
        if axis == 1:
            for atom in structure.atoms:
                atom.xy += distance
        if axis == 2:
            for atom in structure.atoms:
                atom.xz += distance
        return structure

    def translate(self, structure_, axis, distance):
        """
        Parameters
        --------
        structure: parmed.structure
        axis: int
            0: x axis
            1: y axis
            2: z axis
        distance: float
            sign like +,- means direction
        target: string
            'graphene'
            'cdc'
            'all'
            the objective that you want to move
        """
        structure = structure_.__copy__()
        if axis == 0:
            for atom in structure.atoms:
                atom.xx += distance
        if axis == 1:
            for atom in structure.atoms:
                atom.xy += distance
        if axis == 2:
            for atom in structure.atoms:
                atom.xz += distance
        return structure

    def combine(self):
        """ combine cdc and graphene
        """
        return self._cdc + self._graphene

    def calculate_n(self, target, MW, density, volume, sol_ratio = 0):
        """calculate the estimated number of molecules for a certain volume

        Parammeters
        --------
        target : list of string
            Must be solvent first and then solutes
            Target molecules or combinations like acn, litfsi, emimtfsi  
        MW : dictionary
            molecular weight, g/mol; the dictionary should have the items in target
        density : dictionary
            the dictionary should have the items in target,  g/cm^(3)
        volume : float
            the volume for the solution, nm^3
        sol_ratio: int
            the ratio for solvent/ionic_liquid
        """
        import scipy.constants as constants
        import sympy as sym
        # How many nm^3 for 1 particle, and unit has been changed to nm^3/particle for V_per_molecule or unit_V
        unit_V = []
        for i in range(len(target)):
            V_per_molecule = MW[target[i]]/density[target[i]] * 10 ** (21) / (constants.Avogadro)
            unit_V.append(V_per_molecule)

        if sol_ratio == 0:
            x =  sym.symbols('x') # x is the number of molecules in a box
            solution = sym.solve([ 
                            sym.Eq(x * unit_V[0], volume)
                            ])
        else:
            x, y = sym.symbols('x y') # x is the number of solvent molecules in a box, y is the number of solutes in a box
            solution = sym.solve([ sym.Eq(sol_ratio, x/y),
                                sym.Eq(x * unit_V[0]+ y * unit_V[1], volume)
                                ])
        solution_list = list(solution.values())
        n_list = [int(num) for num in solution_list]
        return n_list