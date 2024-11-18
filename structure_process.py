from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from copy import deepcopy

class Structure_process:
    # Initialization
    def __init__(self, stru_file):
        self.stru = Structure.from_file(stru_file)
        self.element = self.stru.species
        self.atom_info = []

        # Structure Analysis
        sg_analyzer = SpacegroupAnalyzer(self.stru, symprec=0.1, angle_tolerance=5.0)
        syms = sg_analyzer.get_symmetrized_structure()
        self.eq_ids = syms.equivalent_indices

    def replace_atom(self, init_atom, new_atom):
        ids_atom = self.stru.indices_from_symbol(init_atom)
        eq_ids_atom = [[i for i in A if i in ids_atom] for A in self.eq_ids]
        eq_ids_atom = [sublist for sublist in eq_ids_atom if sublist]
        replaced_structure = []
        for atom_list in eq_ids_atom:
            atom_id = atom_list[0]
            atom_site = self.stru[atom_id]
            atom_coords = atom_site.frac_coords

            stru_copy = self.stru.copy()
            stru_copy.replace(atom_id, species = new_atom)
            replaced_structure.append((stru_copy, atom_id))

        # with open('defect_info.txt', 'w') as f:
        #     f.write(','.join(map(str, self.atom_info)))
        return replaced_structure
    
    def remove_atom(self, init_atom):
        ids_atom = self.stru.indices_from_symbol(init_atom)
        eq_ids_atom = [[i for i in A if i in ids_atom] for A in self.eq_ids]
        eq_ids_atom = [sublist for sublist in eq_ids_atom if sublist]
        removed_structure = []
        for atom_list in eq_ids_atom:
            atom_id = atom_list[0]
            atom_site = self.stru[atom_id]
            atom_coords = atom_site.frac_coords
            self.atom_info.append(atom_coords)

            stru_copy = self.stru.copy()
            stru_copy.remove_sites([atom_id])
            removed_structure.append((stru_copy, atom_id))

        # with open('defect_info.txt', 'w') as f:
        #     f.write(','.join(map(str, self.atom_info)))
        return removed_structure
