from structure_process import Structure_process
from calculation_process import Calculation_process
from post_process import Post_process
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Outcar
import json
import argparse
import os
import warnings
import time

warnings.filterwarnings("ignore")

def wait_for_complete():
    check = True
    while check:
        os.system("qstat -a > pbs_now.log")
        with open("pbs_now.log", "r") as file:
            lines = file.readlines()
        
        cnt = 0
        for i in range(5, len(lines)):
            if lines[i][99] != 'C':
                cnt += 1
        
        if cnt == 0:
            check = False
        else:
            time.sleep(30)

def submit_jobs(base_directory, waiting = True):
    cur_dir = os.getcwd()
    for root, dirs, files in os.walk(base_directory):
        if not files:
            continue
        # submitted
        if 'submitted' in files:
            continue

        check = True

        while check:
            os.system("qstat -a > pbs_now.log")
            with open("pbs_now.log", "r") as file:
                lines = file.readlines()
            
            cnt = 0
            for i in range(5, len(lines)):
                if lines[i][99] != 'C':
                    cnt += 1
            
            if cnt < 6:
                check = False
            else:
                time.sleep(5)

        os.chdir(root)
        os.system("qsub submit.sh")
        with open("submitted", "w") as file:
            pass
        print(f"Job {root} submitted!")
        os.chdir(cur_dir)

    if waiting:
        wait_for_complete()

def read_input(input_file):
    print("Reading input parameters...")
    with open(input_file, "r") as file:
        input_param = json.load(file)

    print("Successfully Read input parameters!")
    return input_param

def main(input_param):
    # Define data
    vbm = 3.404998
    e_r = [[10.2, 0, -0.1399],[0, 10.21, 0],[-0.1399, 0, 13.12]]
    chem_pot_O = -4.92997145
    chem_pot_N = -8.4823571
    chem_pot_Ga = -7.5189547
    chem_pot_P = -3

    # Read initial structure
    structure_file = os.path.join(input_param['structure']['path'], input_param['structure']['file'])
    stru = Structure_process(structure_file)
    print("Successfully Read initial structures!")

    # Read defects
    for defects_conf in input_param['defects']:

        # Get defects structures
        if defects_conf['type'] == 'replace':
            init_atom = defects_conf['init_atom']
            new_atom = defects_conf['new_atom']

            new_stru = stru.replace_atom(init_atom, new_atom)
            defect_name = f"{new_atom}_{init_atom}"
        elif defects_conf['type'] == 'remove':
            init_atom = defects_conf['init_atom']

            new_stru = stru.remove_atom(defects_conf['init_atom'])
            defect_name = f"V_{init_atom}"

        num_stru = len(new_stru)
        print(f"Created {num_stru} new defects structures!")

        # Structure optimization
        for id, (stru, atom_id) in enumerate(new_stru):
            # Structure opt (PBE) -> SCF (PBE) -> Formation energy and transition levels -> SCF (HSE06)
            opt_dir = os.path.join(f'{input_param['system_name']}', 'opt_defects', defect_name + str(id))
            # Structure optimization
            opt_calculator = Calculation_process(stru, opt_dir, atom_id)
            opt_calculator.opt_pbe(input_param['submit_queue'])
            print(f"No. {id} Structure optimization files created!")

        # Waiting for the calculation completed
        submit_jobs(os.path.join("calc", input_param['system_name'], 'opt_defects'), True)
        print(f"Structure optimization completed!")

        # SCF (PBE)
        if input_param['calculation']['scf'] == 'PBE':
            for id, (stru, atom_id) in enumerate(new_stru):
                opt_dir = os.path.join(f'{input_param['system_name']}', 'opt_defects', defect_name + str(id))
                # Extract the charge
                opt_outcar = Outcar(os.path.join("calc", opt_dir, "OUTCAR"))
                total_electrons = int(opt_outcar.nelect)        

                # SCF (PBE)
                scf_dir = os.path.join(f'{input_param['system_name']}', 'data_defects', defect_name + str(id))
                scf_stru = Structure.from_file(os.path.join("calc", opt_dir,'CONTCAR'))

                for charge in defects_conf['charge']:
                    charge_dir = os.path.join(scf_dir, str(charge))
                    scf_calculator = Calculation_process(scf_stru, charge_dir, atom_id)
                    scf_calculator.scf_pbe(input_param['submit_queue'], total_electrons - charge)
                print(f"No. {id} SCF files created!")
                submit_jobs(os.path.join("calc", input_param['system_name'], 'data_defects'), False)
        elif input_param['calculation']['scf'] == 'HSE06':
            for id, (stru, atom_id) in enumerate(new_stru):
                opt_dir = os.path.join(f'{input_param['system_name']}', 'opt_defects', defect_name + str(id))
                # Extract the charge
                opt_outcar = Outcar(os.path.join("calc", opt_dir, "OUTCAR"))
                total_electrons = int(opt_outcar.nelect)        

                # SCF (PBE)
                scf_dir = os.path.join(f'{input_param['system_name']}', 'data_defects', defect_name + str(id) + "_HSE06")
                scf_stru = Structure.from_file(os.path.join("calc", opt_dir,'CONTCAR'))

                for charge in defects_conf['charge']:
                    charge_dir = os.path.join(scf_dir, str(charge))
                    scf_calculator = Calculation_process(scf_stru, charge_dir, atom_id)
                    scf_calculator.scf_hse06(input_param['submit_queue'], total_electrons - charge) 
                print(f"No. {id} SCF files created!")
                submit_jobs(os.path.join("calc", input_param['system_name'], 'data_defects'), False)

        wait_for_complete()
        print("SCF completed!")
        # Calculate formation energy
        chemical_potentials = {"N": chem_pot_N, 'O': chem_pot_O, 'Ga': chem_pot_Ga, 'P': chem_pot_P}
        post_dir = os.path.join("calc", f'{input_param['system_name']}')
        postprocess = Post_process(post_dir)
        postprocess.calculate_formation_energy_system(post_dir, 'ko', True, vbm, e_r, chemical_potentials)

        # Calculate transition levels
        postprocess.calculate_transition_level_system(post_dir, 'ko', True, vbm, e_r, chemical_potentials, (input_param['post_process']['gap_range'][0], input_param['post_process']['gap_range'][1]))
        print("Post process completed!")



        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="inputs")

    # Define parameters
    parser.add_argument("-i", "--input", type = str, help = "directory of input parameters json file, required.", required=True)

    # Resolve
    args = parser.parse_args()

    input_param = read_input(args.input)

    main(input_param)

    print("Completed!")