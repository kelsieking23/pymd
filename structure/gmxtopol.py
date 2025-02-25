import os
import numpy as np
from pymd.structure.structure_file import StructureFile

'''
We are going to need:
1) atom data
2) bond data 
   - what is bonded to what topol.top / [ bonds ]
   
'''
class Topology:

    def __init__(self, inp, crd=None) -> None:
        self.inp = inp
        if crd is not None:
            self.sf = StructureFile(crd)
        else:
            self.sf = None
        self.included = self.get_includes()
        self.molecules = self.get_molecules()
        # self.atoms = self.get_atoms()
        # self.bonds = self.get_bonds()

    def get_includes(self, inp=None, include_files=None):
        '''
        Get #included files in the topology. 
        '''
        if include_files is None:
            include_files = []
        if inp is None:
            inp = self.inp
        with open(inp, 'r') as f:
            for line in f:
                if line.startswith(';'):
                    continue
                if (line.startswith('#')) and ('include' in line):
                    filename = line.strip().split('include')[1].strip('"').strip("'").strip()[1:]
                    filepath = os.path.join(os.path.dirname(inp), filename)
                    if os.path.isfile(filepath):
                        include_files.append(filepath)
                        include_files = self.get_includes(filepath, include_files)
                    else:
                        raise ValueError(f'Error in topology {inp}: No such file or directory {filepath}')
        return include_files

    def get_moleculetype(self, inp):
        at_directive = False
        with open(inp, 'r') as f:
            for line in f:
                if line.startswith(';'):
                    continue
                elif line.startswith('[ '):
                    if line.startswith('[ moleculetype ]'):
                        at_directive = True
                    else:
                        at_directive = False
                else:
                    if at_directive:
                        if len(line.strip()) > 1:
                            return line.split()[0]
                    else:
                        continue
        return None

    def get_molecules(self):
        molecules = {}
        with open(self.inp, 'r') as f:
            at_directive = False
            for line in f:
                if line.startswith(';'):
                    continue
                elif line.startswith('[ '):
                    if line.startswith('[ molecules ]'):
                        at_directive = True
                    else:
                        at_directive = False
                else:
                    if at_directive:
                        if len(line.strip()) > 1:
                            name, n = tuple(line.strip().split())
                            molecule = Molecule(name, n, self)
                            molecules[name] = molecule
                    else:
                        continue
        return molecules

    def get_molecule(self, molecule_name):
        return self.molecules[molecule_name]
    
    @property
    def molecule_names(self):
        return [molecule.name for molecule in self.molecules.values()]

    @property
    def atom_types(self):
        return [atom.type for atom in self.atoms]

    def get_atom(self, idx):
        '''
        Get atom object from its 1-indexed position in the topology. 
        '''
        return self.atoms[idx-1]
    
    def get_atoms_by(self, strict=True, **kwargs) -> list:
        '''
        get atom(s) by an index, name, or other attribute(s)
        Arguments:
        * strict (bool, default=True): if True, return only atoms that satisfy all the listed criteria. 
                                       if false, return any atom that meets any of the listed criteria.
        * kwargs: any key-value pair, where the key is a valid attribute name of the Atom class. 
        '''
        atoms = []
        for k in kwargs.keys():
            if k not in self.atoms[0].__dict__.keys():
                raise ValueError(f'class Atom has no such attribute: {k}')
        for atom in self.atoms:
            meets_criteria = True
            for k, v in kwargs.items():
                if strict:
                    if atom.__dict__[k] != v:
                        meets_criteria = False
                else:
                    if atom.__dict__[k] == v:
                        meets_criteria = True
            if meets_criteria:
                atoms.append(atom)
        return atoms

    def set_atom(self, idx, new): #type: ignore
        self.atoms[idx-1] = new

    def get_atoms(self):
        '''
        Get entries in topology [ atoms ] directive. 
        '''
        atoms = self.parse_atoms(self.inp)
        for file in self.included:
            atoms = atoms + self.parse_atoms(file)
        return atoms

    def parse_atoms(self, inp, check_moleculetype=True):
        '''
        Gets atom entries in [ atoms ] directive of a given input (topol.top, *inp)
        '''
        atoms = []
        lines = self.parse_itp(inp, 'atoms', check_moleculetype=check_moleculetype)
        for line in lines:
            if len(line.strip().split()) == 0:
                continue
            if not line.strip().split()[0].isnumeric():
                continue
            else:
                args = line.strip().split()[:7]
                args.append(self)
                atom = Atom(*args)
                atoms.append(atom)
        return atoms

    def parse_itp(self, inp, directive, check_moleculetype=True):
        '''
        Get lines from .top/.itp associated with a specified directive
        '''
        directive_lines = []
        if check_moleculetype:
            moleculetype = self.get_moleculetype(inp)
            if moleculetype not in self.molecule_names:
                return []
            else:
                pass
        with open(inp, 'r') as f:
            at_directive = False
            for line in f:
                if line.startswith(';'):
                    continue
                elif line.startswith('[ '):
                    if line.startswith(f'[ {directive} ]'):
                        at_directive = True
                    else:
                        at_directive = False
                else:
                    if at_directive:
                        if len(line.strip().split()) == 0:
                            continue
                        else:
                            directive_lines.append(line)
        return directive_lines
    
    def get_bonds(self):
        '''
        Get entries in topology [ bonds ] directive. 
        '''
        bonds = self.parse_bonds(self.inp)
        for file in self.included:
            bonds = bonds + self.parse_bonds(file)
        return bonds 
    
    def parse_bonds(self, inp):
        bonds = []
        lines = self.parse_itp(inp, 'bonds')
        for line in lines:
            if len(line.strip().split()) == 0:
                continue
            if not line.strip().split()[0].isnumeric():
                continue
            else:
                ii, ij, funct = tuple(map(int, line.strip().split()[0:3]))
                ai = self.get_atom(ii)
                aj = self.get_atom(ij)
                if aj not in ai.bonded_to:
                    ai.bonded_to.append(aj)
                if ai not in aj.bonded_to:
                    aj.bonded_to.append(ai)
                self.set_atom(ii, ai)
                self.set_atom(ij, aj)
                bonds.append(Bond(ii, ij, ai, aj, funct))
        self.bonds = bonds
        self.get_bond_parameters()
        return self.bonds
    
    def get_bond_parameters(self):
        files = []
        for file in self.included:
            base, _ = os.path.splitext(os.path.basename(file))
            if base.endswith('ffbonded'):
                files.append(file)
        ff_bonded = self.get_ff_bondtypes(files)
        for idx, bond in enumerate(self.bonds):
            i = bond.ai.type
            j = bond.aj.type
            b0 = ff_bonded[i][j]['b0']
            kb = ff_bonded[i][j]['kb']
            bond.b0 = float(b0)
            bond.kb = float(kb)
            self.bonds[idx] = bond

    def get_ff_bondtypes(self, inps):
        bondtypes = {}
        '''
        structure of bondtypes should be:
        bondtypes = {
            ai:{
                aj: {
                    b0:float,
                    kb:float
                }
            }
        }
        '''
        for inp in inps:
            lines = self.parse_itp(inp, 'bondtypes', check_moleculetype=False)
            if lines == []:
                continue
            for line in lines:
                i, j, func, b0, kb = tuple(line.strip().split())
                if (i in self.atom_types) and (j in self.atom_types):
                    func = int(func)
                    b0 = float(b0)
                    kb = float(kb)
                    if i not in bondtypes.keys():
                        bondtypes[i] = {}
                    if j not in bondtypes.keys():
                        bondtypes[j] = {}
                    if j not in bondtypes[i].keys():
                        bondtypes[i][j] = {}
                    if i not in bondtypes[j].keys():
                        bondtypes[j][i] = {}
                    bondtypes[i][j]['b0'] = b0
                    bondtypes[j][i]['b0'] = b0
                    bondtypes[i][j]['kb'] = kb
                    bondtypes[j][i]['kb'] = kb
        return bondtypes
            
    def get_bond_stretch_e(self):
        return sum([bond.e for bond in self.bonds])



class Molecule:

    def __init__(self, name: str, n: int, parent: Topology):
        self.name = name
        self.n = int(n)
        self.parent = parent
        self.include = self._get_include()
        self.atoms = []
        if self.include is not None:
            self.atoms = self._get_atoms()

    def _get_include(self):
        for include in self.parent.included:
            if (self.parent.get_moleculetype(include) == self.name):
                return include
    
    def _get_atoms(self):
        return self.parent.parse_atoms(self.include, check_moleculetype=False)

    



class Atom:

    def __init__(self, *args) -> None:
        self.number = int(args[0])
        self.idx = int(args[0])
        self.type = str(args[1])
        self.residue_number = int(args[2])
        self.residue_name = str(args[3])
        self.name = str(args[4])
        self.cgnr = int(args[5])
        self.charge = float(args[6])
        # self.mass = float(args[7])
        self.parent = args[7]
        self.bonded_to = []
    
    def get_position(self):
        if self.parent.sf is not None:
            for line in self.parent.sf.read():
                if (line.atom_name == self.name):
                    return line.x/10, line.y/10, line.z/10
        return None
            

class Bond:

    def __init__(self, ii: int, ij: int, ai: Atom, aj: Atom, funct: int) -> None:
        self.ii = int(ii)
        self.ij = int(ij)
        self.ai = ai
        self.aj = aj
        self.funct = int(funct)
        self.b0 = 0.0
        self.kb = 0.0

    @property
    def e(self):
        return 0.5 * self.kb * (self.length - self.b0)**2

    @property
    def length(self):
        return self.distance(self.ai, self.aj)

    @staticmethod
    def distance(ai, aj):
        ix, iy, iz = ai.get_position()
        jx, jy, jz = aj.get_position()
        return np.sqrt((ix-jx)**2 + (iy - jy)**2 + (iz - jz)**2)


class Angle:

    def __init__(self) -> None:
        pass
        
class Dihedral:

    def __init__(self) -> None:
        pass

class BondedFF:

    def __init__(self):
        self.bondtypes = []
        self.constrainttypes = []
        self.angletypes = []
        self.dihedraltypes = []

if __name__ == '__main__':
    top = Topology('D:/Work/silcsbio.2023.1/examples/cgenff/cgenff_parameters_example/PBZ_gmx.top',
                'D:/Work/silcsbio.2023.1/examples/cgenff/cgenff_parameters_example/PBZ_gmx.pdb')
    # TODO: need angle bend energy