from pymatgen.io.vasp.outputs import Oszicar, Outcar, Locpot
from pymatgen.core import Structure
import os
from spinney.structures.defectivesystem import DefectiveSystem
from spinney.io.vasp import extract_potential_at_core_vasp
from spinney.tools.formulas import count_elements

import ase.io

# Post process for files
# Energy, Potential

class Post_process:
    # Initialization
    def __init__(self, dir):
        self.dir = dir
    
    def get_energy(self):
        self.outcar = Outcar(os.path.join(self.dir,"OUTCAR"))
        return self.outcar.final_energy
    
    def get_potential(self, axis, position):
        self.locpot = Locpot.from_file(os.path.join(self.dir,"LOCPOT"))
        potential_data = self.locpot.get_average_along_axis(axis)
        print(self.locpot.get_axis_grid(axis))
        if position == 0:
            return potential_data[0]
        elif position == 1:
            return potential_data[-1]
        
    def calculate_formation_energy_system(self, system_dir, correction_type, write_file, vbm, dielectric_tensor, chemical_potentials: dict):
        """
        The directory should be in the format:
        - data_defects
            - defect_name1
                - charge1
                - charge2
                ...
            - defect_name2
            ...
        - pristine
        """
        # Define directory
        defective_system = DefectiveSystem(system_dir, 'vasp')

        # Define data
        defective_system.vbm = vbm
        defective_system.dielectric_tensor = dielectric_tensor
        defective_system.chemical_potentials = chemical_potentials

        # Define correction
        defective_system.correction_scheme = correction_type

        # Do the calculation
        defective_system.calculate_energies(verbose = True)

        # Write the file
        if write_file == True:
            defective_system.write_formation_energies(f'{system_dir}/transition_levels.txt.txt')

        return defective_system.data
    
    def calculate_transition_level_system(self, system_dir, correction_type, write_file, vbm, dielectric_tensor, chemical_potentials: dict, gap_range):
        """
        The directory should be in the format:
        - data_defects
            - defect_name1
                - charge1
                - charge2
                ...
            - defect_name2
            ...
        - pristine
        """
        # Define directory
        defective_system = DefectiveSystem(system_dir, 'vasp')

        # Define data
        defective_system.vbm = vbm
        defective_system.dielectric_tensor = dielectric_tensor
        defective_system.chemical_potentials = chemical_potentials
        defective_system.gap_range = gap_range

        # Define correction
        defective_system.correction_scheme = correction_type

        # Do the calculation
        defective_system.calculate_energies(verbose = True)

        # Calculate transition levels
        transition_levels = defective_system.diagram.transition_levels

        # Write the file
        if write_file == True:
            with open(f"{system_dir}/transition_levels.txt", "w") as f:
                f.write(transition_levels.to_string())

        return transition_levels
