from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from copy import deepcopy
import numpy as np

def calculate_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

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
        for id, atom_list in enumerate(eq_ids_atom):
            atom_id = atom_list[0]
            atom_site = self.stru[atom_id]
            atom_coords = atom_site.frac_coords

            stru_copy = self.stru.copy()
            stru_copy.replace(atom_id, species = new_atom)
            defect_name = f'{new_atom}_{init_atom}{id}'
            replaced_structure.append((stru_copy, atom_id, defect_name))

        # with open('defect_info.txt', 'w') as f:
        #     f.write(','.join(map(str, self.atom_info)))
        return replaced_structure
    
    def remove_atom(self, init_atom):
        ids_atom = self.stru.indices_from_symbol(init_atom)
        eq_ids_atom = [[i for i in A if i in ids_atom] for A in self.eq_ids]
        eq_ids_atom = [sublist for sublist in eq_ids_atom if sublist]
        removed_structure = []
        for id, atom_list in enumerate(eq_ids_atom):
            atom_id = atom_list[0]
            atom_site = self.stru[atom_id]
            atom_coords = atom_site.frac_coords
            self.atom_info.append(atom_coords)

            stru_copy = self.stru.copy()
            stru_copy.remove_sites([atom_id])
            defect_name = f'V_{init_atom}{id}'
            removed_structure.append((stru_copy, atom_id, defect_name))

        # with open('defect_info.txt', 'w') as f:
        #     f.write(','.join(map(str, self.atom_info)))
        return removed_structure

    def multi_defects(self, defects_list: list):
        defects_structure = []

        defect_conf_0 = defects_list[0]

        init_atom = defect_conf_0[1]
        ids_atom = self.stru.indices_from_symbol(init_atom)
        eq_ids_atom = [[i for i in A if i in ids_atom] for A in self.eq_ids]
        eq_ids_atom = [sublist for sublist in eq_ids_atom if sublist]

        for id, atom_list in enumerate(eq_ids_atom):
            atom_id = atom_list[0]

            stru_copy = self.stru.copy()
            if defect_conf_0[0] == 'remove':
                stru_copy.remove_sites([atom_id])
                defect_name = f'V_{init_atom}{id}'
            elif defect_conf_0[0] == 'replace':
                new_atom = defect_conf_0[2]
                stru_copy.replace(atom_id, species = new_atom)
                defect_name = f'{new_atom}_{init_atom}{id}'
            defects_structure.append((stru_copy, atom_id, defect_name))

        for id, defect_conf in enumerate(defects_list):
            if id == 0:
                continue

            init_atom = defect_conf[1]
            ids_atom = self.stru.indices_from_symbol(init_atom)
            eq_ids_atom = [[i for i in A if i in ids_atom] for A in self.eq_ids]
            eq_ids_atom = [sublist for sublist in eq_ids_atom if sublist]

            tmp_structure = []

            for (stru, pre_atom, defect_name) in defects_structure:
                for id, atom_list in enumerate(eq_ids_atom):
                    sorted_atom_list = sorted(
                        atom_list,
                        key = lambda atom_id: calculate_distance(self.stru[atom_id].frac_coords, self.stru[pre_atom].frac_coords)
                    )
                    atom_id = sorted_atom_list[0]
                    
                    stru_copy = stru.copy()
                    if defect_conf[0] == 'remove':
                        stru_copy.remove_sites([atom_id])
                        defect_name_new = defect_name + f"V_{init_atom}{id}"
                    elif defect_conf[0] == 'replace':
                        new_atom = defect_conf[2]
                        stru_copy.replace(atom_id, species = new_atom)
                        defect_name_new = defect_name + f"{new_atom}_{init_atom}{id}"

                    tmp_structure.append((stru_copy, atom_id, defect_name_new))
            defects_structure = tmp_structure

        return defects_structure

if __name__ == "__main__":
    stru = Structure_process("POSCAR")
    defects_list = [["replace", "O", "N"], ["remove", "Ga"]]
    new_stru = stru.multi_defects(defects_list)
    
    test_stru = new_stru[2]
    test_stru[0].to("POSCAR_test", fmt = 'POSCAR')
    print(len(new_stru))