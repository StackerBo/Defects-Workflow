import os
import shutil
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPStaticSet

class Calculation_process:
    # Initilization
    def __init__(self, structure: Structure, dir, atom_info):
        self.stru = structure
        self.elements = structure.species
        self.dir = dir
        atom_site = self.stru[atom_info]
        self.atom_info = atom_site.frac_coords

    def scf_pbe(self, queue, charge):
        # create the dir
        cal_dir = os.path.join("calc",self.dir)
        if not os.path.exists(cal_dir):
            os.makedirs(cal_dir)
        vasp_input_set = MPStaticSet(self.stru)
        vasp_input_set.write_input(cal_dir)

        # Delete INCAR and KPOINTS
        os.remove(os.path.join(cal_dir,"INCAR"))
        os.remove(os.path.join(cal_dir,"KPOINTS"))

        shutil.copy("INCAR/SCF_PBE", os.path.join(cal_dir,"INCAR"))
        shutil.copy("KPOINTS/KPOINTS", cal_dir)

        # change the charge
        with open(os.path.join(cal_dir,"INCAR"), "r") as file:
            lines = file.readlines()

        with open(os.path.join(cal_dir,"INCAR"), "w") as file:
            for line in lines:
                if line.strip().startswith("NELECT"):
                    line = f"NELECT = {charge}\n"
                file.write(line)

        submit_file = os.path.join(cal_dir, "submit.sh")
        if "f96" in queue:
            shutil.copy("SUBMIT/f96.sh", submit_file)
        elif "std2" in queue:
            shutil.copy("SUBMIT/std2.sh", submit_file)
        elif "std40" in queue:
            shutil.copy("SUBMIT/std40.sh", submit_file)

        # Enter the path
        original_cwd = os.getcwd()

        os.chdir(cal_dir)
        with open('position.txt', 'w') as f:
            f.write(' '.join(map(str, self.atom_info)))
        # os.system("qsub submit.sh")
        # print("Successfully submit to the queue!")

        os.chdir(original_cwd)


    def scf_hse06(self, queue, charge):
        # create the dir
        cal_dir = os.path.join("calc",self.dir)
        if not os.path.exists(cal_dir):
            os.makedirs(cal_dir)
        vasp_input_set = MPStaticSet(self.stru)
        vasp_input_set.write_input(cal_dir)

        # Delete INCAR and KPOINTS
        os.remove(os.path.join(cal_dir,"INCAR"))
        os.remove(os.path.join(cal_dir,"KPOINTS"))

        shutil.copy("INCAR/SCF_HSE06", os.path.join(cal_dir,"INCAR"))
        shutil.copy("KPOINTS/KPOINTS", cal_dir)

        # change the charge
        with open(os.path.join(cal_dir,"INCAR"), "r") as file:
            lines = file.readlines()

        with open(os.path.join(cal_dir,"INCAR"), "w") as file:
            for line in lines:
                if line.strip().startswith("NELECT"):
                    line = f"NELECT = {charge}\n"
                file.write(line)

        submit_file = os.path.join(cal_dir, "submit.sh")
        if "f96" in queue:
            shutil.copy("SUBMIT/f96.sh", submit_file)
        elif "std2" in queue:
            shutil.copy("SUBMIT/std2.sh", submit_file)
        elif "std40" in queue:
            shutil.copy("SUBMIT/std40.sh", submit_file)

        # Enter the path
        original_cwd = os.getcwd()

        os.chdir(cal_dir)
        with open('position.txt', 'w') as f:
            f.write(' '.join(map(str, self.atom_info)))
        # os.system("qsub submit.sh")
        # print("Successfully submit to the queue!")

        os.chdir(original_cwd)

    def opt_pbe(self, queue):
        # create the dir
        cal_dir = os.path.join("calc",self.dir)
        if not os.path.exists(cal_dir):
            os.makedirs(cal_dir)
        vasp_input_set = MPStaticSet(self.stru)
        vasp_input_set.write_input(cal_dir)

        # Delete INCAR and KPOINTS
        os.remove(os.path.join(cal_dir,"INCAR"))
        os.remove(os.path.join(cal_dir,"KPOINTS"))

        shutil.copy("INCAR/OPT_PBE", os.path.join(cal_dir,"INCAR"))
        shutil.copy("KPOINTS/KPOINTS", cal_dir)

        submit_file = os.path.join(cal_dir, "submit.sh")
        if "f96" in queue:
            shutil.copy("SUBMIT/f96.sh", submit_file)
        elif "std2" in queue:
            shutil.copy("SUBMIT/std2.sh", submit_file)
        elif "std40" in queue:
            shutil.copy("SUBMIT/std40.sh", submit_file)

        # Enter the path
        original_cwd = os.getcwd()

        os.chdir(cal_dir)
        with open('position.txt', 'w') as f:
            f.write(' '.join(map(str, self.atom_info)))
        # os.system("qsub submit.sh")
        # print("Successfully submit to the queue!")

        os.chdir(original_cwd)