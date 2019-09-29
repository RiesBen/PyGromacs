"""
FUNCTIONLIB:            generalFunctions_simulations
Description:
   This lib contains funtions for cphmd-simulations in a workflow.

Author: Benjamin Schroeder
"""
import glob
import os
import random

from generalutilities.function_libs.gromacs import top as top
from generalutilities.function_libs.utils import bash_commands as bash
from generalutilities.function_libs.gromacs import parser as parser
from generalutilities.function_libs.gromacs import gro as gro


#########################
# Important Lambda_files
# _____________________________

# lambda_groups.dat
def build_lambda_file(output_folder, gro_file, initial_lambdas= [], options=False, barrier=False, his_atoms=False, h3o_atoms=False, residues=False, his_3state=True ):
    """
    Args:
        output_folder: STRING path to folder
        gro_file: STRING path to gro file
        initial_lambdas: LIST(float) of initial lambda numbers
        options: DICT
        barrier: FLOAT setting the simulation energy barrier
        his_atoms: LIST(int)
        h3o_atoms: LIST(int)
        residues: LIST(int)
        his_3state:
    """
    print("Building Lambda file")
    if(options):
        his_atoms = options["his_atoms"]
        h3o_atoms = options["h3o_atoms"]
        residues = options["residues"]
        barrier = options["cphmd_barrier"]

    if(barrier and his_atoms and h3o_atoms and residues):
            raise Exception("There are not all necessary arguments set for build_lambda_file.")

    print("importatnt his: " + str(his_atoms))
    gro_lines = parser.read_gro(gro_file + ".gro")
    gro_coord = gro.clean_fields(gro_lines["coordinates"])
    lambda_file = open(output_folder + "/lambda_groups.dat", "w")

    print("start_writing lambda_groups.dat")
    print("got: "+str(initial_lambdas))
    if len(initial_lambdas) == 0:
        initial_lambdas = ["0" for x in range(2 * len(residues))]    # init with double protonated all residues
        initial_lambdas.append("1") #set H3O to opposite of hisitidines

    print("init_l_values:" + str(initial_lambdas))
    init_lambda_ndx = 0
    for residue in residues:
        atoms_list = gro.get_residue_atoms(residue, "HIS", gro_coord, his_atoms)
        atoms = ""
        for x in atoms_list:
            atoms += str(x) + " "
        print(residue + " found: " + atoms)
        print("index "+ str(init_lambda_ndx))
        print("hu"+ str(initial_lambdas[init_lambda_ndx]))
        HIA_line = ("name = HI2A\n"
                    "residue_number = " + str(residue) + "\n"
                    "initial_lambda = " + str(initial_lambdas[init_lambda_ndx]) + "\n"
                    "barrier = " + str(barrier) + "\n"
                    "points = 0\n"
                    "interval_adapt_barr = 0\n"
                    "number_of_atoms = 9\n"
                    "" + atoms + "\n\n")
        init_lambda_ndx += 1
        lambda_file.write(HIA_line)
        print("wrote_HIA")
        if(his_3state):
            HIB_line = ("name = HI2B\n"
                        "residue_number =  " + str(residue) + "\n"
                         "initial_lambda = " + str(initial_lambdas[init_lambda_ndx]) + "\n"
                         "barrier = " + str(barrier) + "\n"
                         "points = 0\n"
                         "interval_adapt_barr = 0\n"
                         "number_of_atoms = 9\n"
                         "" + atoms + "\n\n")
            lambda_file.write(HIB_line)
            print("wrote_HIB")
        init_lambda_ndx += 1

    print("wrote_lines")
    atoms_list = gro.get_residue_atoms(False, "H3O", gro_coord, h3o_atoms)
    atoms = ""
    for x in atoms_list:
        atoms += str(x) + " "
    Buffer_line = (  # add buffer
        "name = H3OB\n"
        "residue_number = " + str(int(len(atoms_list) / 6)) + "\n"
        "initial_lambda = " + str(initial_lambdas[init_lambda_ndx]) + "\n"
        "barrier = 0\n"
        "points = 0\n"
        "interval_adapt_barr = 0\n"
        "number_of_atoms = " + str(len(atoms_list)) + "\n"
        "" + atoms + "\n")
    lambda_file.write(Buffer_line)
    lambda_file.close()

def backup_lambda_file(simulation_folder):
    """
    Args:
        simulation_folder:
    """
    print("search for l_groups: ")
    l_dyn_files = glob.glob(simulation_folder + "/l_dynamics_group_*.dat")
    if (len(l_dyn_files) > 0):
        print("Found l_files!\n")
        initial_lambdas = []
        for lx in range(1, (2 * 4) + 2):
            lamb_group = simulation_folder + "/l_dynamics_group_" + str(lx) + ".dat"
            step_num = sum([1 for x in os.listdir(simulation_folder) if ("l_dynamics_group_" + str(lx) + "_step_" in x)])
            print(lamb_group)
            print("found steps: " + str(step_num))
            # SHOULD ONLY DO ONE IMPROVE HEREW!
            if (os.path.exists(lamb_group)):
                print("found " + lamb_group)
                # redef_lambda
                file = open(lamb_group, "r")
                tmp_line = ""
                lines = file.readlines()
                if (len(lines) != 0):
                    print("Found lines: " + str(len(lines)))
                    for line in lines:
                        tmp_line = line
                    print(tmp_line)
                    initial_lambdas.append(tmp_line.split(" ")[1])
                    file.close()
                    # backup
                    print("backup under step number: " + str(step_num))
                    bash.mv_file(lamb_group, simulation_folder + "/l_dynamics_group_" + str(lx) + "_step_" + str(step_num) + ".dat")

                else:
                    print("found-no data in l_dyn_group.")
            elif (os.path.exists(simulation_folder + "/l_dynamics_group_" + str(lx) + "_step_" + str(step_num - 1) + ".dat")):
                print("found " + "l_dynamics_group_" + str(lx) + "_step_" + str(step_num - 1))
                # redef_lambda
                file = open(simulation_folder + "/l_dynamics_group_" + str(lx) + "_step_" + str(step_num - 1) + ".dat", "r")
                tmp_line = ""
                lines = file.readlines()
                if (len(lines) != 0):
                    for line in lines:
                        tmp_line = line
                    print(tmp_line)
                    initial_lambdas.append(tmp_line.split(" ")[1])
                    file.close()
            else:
                print("somethings strange with lambda Files.")
                print(initial_lambdas)
                print("didn't find " + lamb_group)
                print("didn't find  " + "l_dynamics_group_" + str(lx) + "_step_" + str(step_num - 1))
                raise Exception("init lambda_problem")
    return initial_lambdas

# HI2B.dat
def build_HI2B(input_file, mapping_file, output_folder, residues):
    """
    Args:
        input_file:
        mapping_file:
        output_folder:
        residues:
    """
    protein_file = parser.read_top(input_file)
    mapping_file = parser.read_top(mapping_file)
    # Input Files:
    mapped_atoms = []
    mapped_atoms.append("nr  type  charge   mass  typeB    chargeB      massB   \n")
    if (not "atoms" in mapping_file):
        raise Exception("mapping file has no mapping section!")
    print(mapping_file)
    map_list = top.dictionize(mapping_file["atoms"])

    protein_atoms = protein_file["atoms"]
    for atom in protein_atoms:
        if ("HIS" in atom and not str(atom).startswith(";") and (" " + residues[0] in atom)):
            atom_fields = top.clean_fields([atom, ])[0]
            for map_atom in map_list:
                map = map_list[map_atom]
                #print(map)
                if (map["atom"] == atom_fields[4] and map["residue"] in atom_fields[3]):
                    new_atom = "{0:>6} {1:>2}     {2:>9}      {3:^5}     {4:>2}  {5:<11}      {6:^5}  ;\n".format(atom_fields[0], atom_fields[1], map["charge"], map["mass"], map["typeB"], map["chargeB"], map["massB"])
                    # "{0:^7}{1:^11}{2:^11}{3:^11}{4:^5}{5:^11}{6:^11}\n"
                    # "{0:^7}{1:^11}{2:^7}{3:^7}{4:^7}{5:^7}{6:^11}{7:^11}{8:^5}{9:^11}{10:^11} ;\n".format(atom_fields[0], atom_fields[1], atom_fields[2], atom_fields[3], atom_fields[4], atom_fields[5], map["charge"], map["mass"], map["typeB"], map["chargeB"], map["massB"])
                    # {0:^7}{1:^11}{2:^7}{3:^7}{4:^7}{5:^7}{6:^11}{7:^11}{8:^5}{9:^11}{10:^11
                    mapped_atoms.append(new_atom)
                    break

                    #  protein_file["atoms"] = mapped_atoms
    # Output
    output_file = open(output_folder + "/HI2B.top", "w")
    output_file.writelines(mapped_atoms)
    output_file.close()


#########################
# Updating according to cpHMD simulation - Decision
# _____________________________

def get_lambda_from_his_states(his_states, options):
    """translate his states to lambdas :param his_states:

    Args:
        his_states:
        options:
    """
    histidines = options["His_lambda_pairs"]
    lambda_values = {}
    buffer_part=1.0/len(histidines)
    buffer_num = 0
    for his in his_states:
        if his in histidines:
            if(his_states[his] == "HIP"):
                l_dims = histidines[his]
                lambda_values.update({int(l_dims[0][1])-1: 0.0})
                lambda_values.update({int(l_dims[1][1])-1: 0.0})
                buffer_num += buffer_part
            elif (his_states[his] == "HIE"):
                l_dims = histidines[his]
                lambda_values.update({int(l_dims[0][1])-1: 1.0})
                lambda_values.update({int(l_dims[1][1])-1: 1.0})
            elif (his_states[his] == "HID"):
                l_dims = histidines[his]
                lambda_values.update({int(l_dims[0][1])-1: 1.0})
                lambda_values.update({int(l_dims[1][1])-1: 0.0})
            else:
                raise Exception("Unknown state here! in his_states->lambda")

    lambda_values.update({(2*len(histidines)): buffer_num}) #buffersite!
    print("generated lambdas: "+str(lambda_values))
    return lambda_values

def translate_decision(decision_dict, options):
    """translate lambdas to his states :param decision_dict: :param options:
    :return:

    Args:
        decision_dict:
        options:
    """
    histidines = options["His_lambda_pairs"]
    residues = options["residues"]
    states = {}
    print("start translating: " + str(histidines))
    for i, x in enumerate(histidines):
        if (decision_dict[histidines[x][0]] == 0):  # HIP state
            states.update({residues[i]: "HIP"})
        else:  # not HIP
            if (decision_dict[histidines[x][1]] == 0):  # HID state
                states.update({residues[i]: "HID"})
            elif (decision_dict[histidines[x][1]] == 1):  # HIE state
                states.update({residues[i]: "HIE"})
    print("found states:" + str(states))
    return states


def top_set_states(top_file, out_top_file, his_states, HIP_HID_top, HIP_HIE_top):
    """
    Args:
        top_file:
        out_top_file:
        his_states: dict with histidine res number to -> state
        HIP_HID_top:
        HIP_HIE_top:
    """
    print(his_states)
    HIP = []
    HID = []
    HIE = []

    print("sorting states")
    for x in his_states:
        if (his_states[x] == "HIP"):
            HIP.append(x)
        elif (his_states[x] == "HID"):
            HID.append(x)
        elif (his_states[x] == "HIE"):
            HIE.append(x)
        else:
            raise Exception("Unknown state!")

    print("map states")
    if (len(HIP) > 0):
        top.map_residue(top_file, HIP_HID_top, out_top_file, HIP, states="A")
    if (len(HID) > 0):
        top.map_residue(out_top_file, HIP_HID_top, out_top_file, HID, states="B")
    if (len(HIE) > 0):
        top.map_residue(out_top_file, HIP_HIE_top, out_top_file, HIE, states="B")


#########################
# gro file modification
# _____________________________
def get_h3o_atom_shift(atom_O_line, h3O_O_line):
    """
    Args:
        atom_O_line:
        h3O_O_line:
    """
    if (len(atom_O_line) == 5 or len(atom_O_line) == 8):
        x_shift = float(atom_O_line[2]) - float(h3O_O_line[3])
        y_shift = float(atom_O_line[3]) - float(h3O_O_line[4])
        z_shift = float(atom_O_line[4]) - float(h3O_O_line[5])
    elif (len(atom_O_line) == 6 or len(atom_O_line) == 9):
        print("SIZE is 9" + str(atom_O_line))
        x_shift = float(atom_O_line[3]) - float(h3O_O_line[3])
        y_shift = float(atom_O_line[4]) - float(h3O_O_line[4])
        z_shift = float(atom_O_line[5]) - float(h3O_O_line[5])
    else:
        raise Exception("Could not calculate atom shift!")
    print(str(x_shift) + "  " + str(y_shift) + "   " + str(z_shift))
    return x_shift, y_shift, z_shift


def move_h3o_on_water(water_line, insert_coordinates, molecules):
    """
    Args:
        water_line:
        insert_coordinates:
        molecules:
    """
    coordinates_H3O = []
    x_shift, y_shift, z_shift = get_h3o_atom_shift(water_line, insert_coordinates[0])
    for atom in insert_coordinates:
        atom[0] = str(molecules) + "H30"
        if (atom[1] == "OW"):  # map O onto water O
            if (len(water_line) <= 5):  # check format of input lines (molten atomtype to atomnum?)
                atom[3] = water_line[2]
                atom[4] = water_line[3]
                atom[5] = water_line[4]
            elif (len(water_line) == 6):
                atom[3] = water_line[3]
                atom[4] = water_line[4]
                atom[5] = water_line[5]
            elif (len(water_line) == 8):  # you have to add velocieties
                print("velocity line =8")
                atom[3] = water_line[2]
                atom[4] = water_line[3]
                atom[5] = water_line[4]
                atom += water_line[5:]
            elif (len(water_line) >= 9):
                print("velocity line =9")
                atom[3] = water_line[3]
                atom[4] = water_line[4]
                atom[5] = water_line[5]
                atom += water_line[6:]

        else:  # adjust Hs to O atom
            if (len(water_line) < 8):  # Do I need to add velocities?
                print("velocity line <8 ")
                print(atom)
                atom[3] = str(round(float(atom[3]) + x_shift, 3))
                atom[4] = str(round(float(atom[4]) + y_shift, 3))
                atom[5] = str(round(float(atom[5]) + z_shift, 3))

            elif (len(water_line) == 8):
                print("velocity line =8 ")
                print(atom)
                atom[3] = str(round(float(atom[3]) + x_shift, 3))
                atom[4] = str(round(float(atom[4]) + y_shift, 3))
                atom[5] = str(round(float(atom[5]) + z_shift, 3))
                atom += water_line[5:]
                print(str(atom) + "\n")
            elif (len(water_line) >= 9):
                print("velocity line >9 ")
                print(atom)
                atom[3] = str(round(float(atom[3]) + x_shift, 3))
                atom[4] = str(round(float(atom[4]) + y_shift, 3))
                atom[5] = str(round(float(atom[5]) + z_shift, 3))
                atom += water_line[6:]
                print(str(atom) + "\n")

        coordinates_H3O.append(atom)
    return coordinates_H3O


def insert_h3O(selection_file, protein_file, insert_mol, output, selection_name="waterSel", insert_times=1):
    # input
    """
    Args:
        selection_file:
        protein_file:
        insert_mol:
        output:
        selection_name:
        insert_times:
    """
    selection_dict = parser.read_ndx_file(selection_file)
    protein_gro = parser.read_gro(protein_file)  #
    protein_coordinates = gro.clean_fields(protein_gro["coordinates"])
    insert_gro = parser.read_gro(insert_mol)

    # get random SOL residues from selection
    # getting precise selection keay
    selection_key = ""
    for x in selection_dict.keys():
        if (selection_name in x):
            selection_key = x
            break

    positions = len(selection_dict[selection_key])
    replace_residues = []
    replace_positions = []
    random.seed()

    for i in range(insert_times):
        x = random.randint(0, positions - 1)
        print("for pos " + str(x) + "\t\t" + str(selection_dict[selection_key][x]))
        replace_positions.append(x)
        replace_residues.append(selection_dict[selection_key][x])
    print("Replace position: \t" + str(replace_positions) + " from total pos: " + str(positions))
    print("\nReplace Residue: \t" + str(replace_residues))

    # get residues, to be replaced
    replace_molecule = []
    for i, line_prot in enumerate(protein_coordinates):
        for x in replace_residues:
            if (i == int(x)):
                print("found: " + str(x) + "\t" + str(line_prot))
                replace_residues.remove(x)
                replace_molecule.append(line_prot[0])

    print("\nreplace mol: " + str(replace_molecule))
    # insert molecules - split coordinates in parts
    print("Map Coordinates on H2Os")
    new_coordinates_H3O = []
    new_coordinates_NA = []
    new_coordinates_CL = []
    new_coordinates_SOL = []
    new_coordinates_rest = []
    molecule = 1
    for line_prot in protein_coordinates:
        if (line_prot[0] in replace_molecule):
            if ("OW" in line_prot[1]):
                print("Found: " + str(line_prot))
                insert_coordinates = gro.clean_fields(insert_gro["coordinates"])
                new_coordinates_H3O += move_h3o_on_water(line_prot, insert_coordinates, molecule)
                molecule += 1
            else:
                print("ignore: " + str(line_prot))
                continue
        elif ("SOL" in line_prot[0]):
            new_coordinates_SOL.append(line_prot)
        elif ("NA" in line_prot[0]):
            new_coordinates_NA.append(line_prot)
        elif ("CL" in line_prot[0]):
            new_coordinates_CL.append(line_prot)
        else:
            new_coordinates_rest.append(line_prot)

    new_coordinates = new_coordinates_rest + new_coordinates_H3O + new_coordinates_NA + new_coordinates_CL + new_coordinates_SOL
    print(" Length of new coordinates: " + str(len(new_coordinates)))
    protein_gro["atoms"] = str(len(new_coordinates))
    protein_gro["coordinates"] = gro.make_clean_gro_atom_lines(new_coordinates)
    gro.write_gro_file(protein_gro, output)

def get_lambda_from_decision(folder):
    """
    Args:
        folder:
    """
    decision_dict = parser.read_decision_file(folder+"/decision.csv")
    lambdas = [decision_dict[x] for x in decision_dict]
    return lambdas

if __name__ == "__main__":
    options = {"residues": ["899", "909", "924", "930"], "cphmd_barrier": "3.0", "his_atoms": ["CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2"], "h3o_atoms": ["OW", "HW1", "HW2", "H31", "H32", "H33"]}
    gro_file = "/home/benjamin/Downloads/1VPR_Sys/output/step_0/1cphmd/0build/build.gro"
    output_folder = "/home/benjamin/Downloads/1VPR_Sys/output/step_0/1cphmd/3cphmd/"
    build_lambda_file(output_folder, options, gro_file)

#______________________________________________________________
#Titration Analysis - own Files pH-Files
#______________________________________________________________
def read_l_dyn_group(path, l_dim="l"):
    """
    Args:
        path:
        l_dim:
    """
    input_file = open(path, "r")
    result_dict = {"t": [], str(l_dim): [], "dVd"+str(l_dim): []}
    for x in input_file.readlines():
        if(len(x.split(" ")) == 3):
            columns = x.split(" ")
            result_dict["t"].append(float(columns[0]))
            result_dict[str(l_dim)].append(float(columns[1]))
            result_dict["dVd"+str(l_dim)].append(float(columns[2].strip()))
        else:
            print("unrecognized line: \n"+x)
    return result_dict


def read_phsim_file(path):
    """reads in the Files from the pHdependent simulations analysis :param
    data_folder: :return:

    Args:
        path:
    """
    ph_file = {}
    # readin CSV
    tmp = open(path, "r")
    lines = tmp.readlines()
    lines = [x.strip("\n").strip().split("\t") for x in lines]
    key_dict = "general"
    ph_file.update({key_dict: {}})
    for x in lines:
        if(len(x) == 1 and not (x[0]== "")):
            key_dict = x[0]

        if len(x) >= 2:
            if(not key_dict in ph_file):
                ph_file.update({key_dict: {}})
            ph_file[key_dict].update({x[0]: x[1]})
    return ph_file

def read_phsim_new_files(data_folder):
    """reads in the Files from the pHdependent simulations analysis :param
    data_folder: :return:

    Args:
        data_folder:
    """
    ph_files = {}

    files = [x for x in os.listdir(data_folder) if x.startswith("ph_") and x.endswith(".csv")]
    # readin CSV
    for ph_file in files:
        print(ph_file)
        tmp = open(data_folder + "/" + ph_file, "r")
        lines = tmp.readlines()
        lines = [x.strip("\n").split("\t") for x in lines]
        file_dict = {}
        tmp_dict = {}
        tmp_name = "general"
        for x in lines:
            if len(x) == 1 and not ("" in x or " " in x):
                if(len(tmp_dict.keys()) > 0):
                    file_dict.update({tmp_name: tmp_dict})
                tmp_dict = {}
                tmp_name = x[0].strip()
            elif len(x) >= 2:
                tmp_dict.update({x[0]: x[1]})
        file_dict.update({tmp_name: tmp_dict})
        ph_files.update({"pH" + file_dict["general"]["pH"]: file_dict})
    return ph_files



def read_sim_out(file_path):
    """
    Args:
        file_path:
    """
    l1_file = open(file_path, "r")
    t_l_dvdl = {"t": [], "l": [], "dVdl": []}

    # read in Files and translate to dictionaries
    for x in l1_file.readlines():
        line = x.strip("\n").split(" ")
        line = [float(x) for x in line]
        t_l_dvdl["t"] += [line[0]]
        t_l_dvdl["l"] += [line[1]]
        t_l_dvdl["dVdl"] += [line[2]]
    l1_file.close()

    return t_l_dvdl

def read_titration_files(data_folder):
    """reads in the analysed Files from a whole simulation . :param data_folder:
    :return:

    Args:
        data_folder:
    """
    csv_files = []
    titration_folders = [x for x in os.listdir(data_folder) if x.startswith("Results_") and os.path.isdir(data_folder + "/" + x)]
    for x in titration_folders:
        csv_files += [data_folder + "/" + x + "/" + y for y in os.listdir((data_folder + "/" + x)) if
                      y.endswith("analysis.csv")]
    print(csv_files)
    # readinCSVs
    titrations_dict = {}
    for csv in csv_files:
        tmp = open(csv, "r")
        lines = tmp.readlines()
        lines = [x.strip("\n").split("\t") for x in lines]
        file_dict = {}
        found_table = False
        table_name = ""
        table_header = True
        for x in lines:
            if len(x) >= 2 and not found_table:
                file_dict.update({x[0]: x[1]})
            elif found_table:
                if x[0] == "</Table:" + table_name + ">":
                    found_table = False
                    continue
                elif table_header:
                    table_header = False
                    file_dict[table_name].update({"header": x})
                    continue
                else:
                    file_dict[table_name].update({x[0]: x[1:]})
            elif len(x) >= 1 and (x[0].startswith("<Table:") and x[0].endswith(">")):
                found_table = True
                table_name = str(x[0]).replace("<Table:", "").replace(">", "")
                file_dict.update({table_name: {}})

        titrations_dict.update({"Simulation_" + file_dict["simulation"]: file_dict})
    return titrations_dict, csv_files, titration_folders




def fetch_values_from_ph_simdict(ph_files, residue="res2"):
    """fetch Values from the lines. used in analyse2DpHtitration &
    analyse2DpHtitrations :param ph_files: :return:

    Args:
        ph_files:
        residue:
    """
    # dataPoints -  visualise:
    data = {"x_ph": []}
    data.update({"y_state_l1": [], "sem_l1": []})
    data.update({"y_state_l2": [], "sem_l2": []})
    data.update({"y_state_l1_with_l2_1": [], "sem_state_l1_with_l2_1": []})
    data.update({"y_state_l1_with_l2_0": [], "sem_state_l1_with_l2_0": []})
    data.update({"y_l1_transitions": [], "y_l2_transitions": []})
    data.update({"HIE_ratios": [], "HID_ratios": [], "HIP_ratios": []})
    data.update({"unprot_states": [], "prot_states": []})
    data.update({"l2_0_states": [], "l2_1_states": []})
    data.update({"total_time": []})
    data.update({"ba_err_est_l1": [], "ba_err_est_l2":[]})

    # fetch Values
    for x in sorted(ph_files):
        data["x_ph"].append(float(ph_files[x]["general"]["pH"]))
        data["total_time"].append(float(ph_files[x]["general"]["total_time"])) #"total_time"])

        data["y_state_l1"].append(float(ph_files[x]["l1"]["state_mean"]))
        data["sem_l1"].append(float(ph_files[x]["l1"]["l_state_sem"]))
        data["ba_err_est_l1"].append(float(ph_files[x]["l1"]["ba_err_est"]))
        data["y_state_l2"].append(float(ph_files[x]["l2"]["state_mean"]))
        data["sem_l2"].append(float(ph_files[x]["l2"]["l_state_sem"]))
        data["ba_err_est_l2"].append(float(ph_files[x]["l2"]["ba_err_est"]))
        data["y_l1_transitions"].append(float(ph_files[x]["l1"]["transitions_total"]))
        data["y_l2_transitions"].append(float(ph_files[x]["l2"]["transitions_total"]))

        data["unprot_states"].append(float(ph_files[x]["l1"]["time_in_state_1"]))
        data["prot_states"].append(float(ph_files[x]["l1"]["time_in_state_0"]))
        data["l2_1_states"].append(float(ph_files[x]["l2"]["time_in_state_1"]))
        data["l2_0_states"].append(float(ph_files[x]["l2"]["time_in_state_0"]))

        data["HIP_ratios"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["HIP_ratio"]))
        data["HID_ratios"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["HID_ratio"]))
        data["HIE_ratios"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["HIE_ratio"]))
        data["y_state_l1_with_l2_1"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["l1_if_l2_1_mean"]))
        data["sem_state_l1_with_l2_1"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["l1_if_l2_1_sem"]))
        data["y_state_l1_with_l2_0"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["l1_if_l2_0_mean"]))
        data["sem_state_l1_with_l2_0"].append(float(ph_files[x]["His_"+residue+"_l1_l2"]["l1_if_l2_0_sem"]))

    return data

def read_decision_file(decision_path):
    """
    Args:
        decision_path:
    """
    file = open(decision_path, "r")
    states = {}
    for line in file.readlines():
        if(line.startswith(";")):
            continue
        else:
            columns = line.split("\t")
            states.update({columns[0]: int(columns[1])})
    return states

#______________________________________________________________
#Titration Analysis - Volume_commandline
#______________________________________________________________
def parse_input_determine_volume(arguments, depend, parser):
    """
    Args:
        arguments:
        depend:
        parser:
    """
    help_text = "python determineVolume.py <file.pdb/gro/xtc> -<target> <task_num> -o <prefix>\n" \
                "e.g.: python determineVolume.py 1VPR.gro -1VPR 12 -o outputdir"
    #######################################################################################################################
    # vars
    pdb_structure = "undefined"
    output_folder = "undefined"
    target = "undefined"
    reference_structure = "undefined"
    protocols = []
    #######################################################################################################################
    # Check Input:
    try:
        # check if all depenencies are available
        bash.check_dependencies(depend)

        # check commandline!
        # check amount of args
        if len(arguments) > 1 and arguments[1] in ["-h", "--h", "-help", "--help"]:
            print(help_text)
            for x in dict(parser).keys():
                if "-" + str(x) in arguments:
                    path = depend["protocol_location"] + "/" + parser[x]["protocol_folder"]
                    print("Available tasks: \n")
                    i = 0
                    for y in os.listdir(path):
                        if x == "ligand" and str(y).startswith("empty"):
                            continue
                        print(str(i) + " " + y)
                        i += 1
            exit()
        elif len(arguments) >= 6:
            # check input-pdb file.
            if os.path.isfile(arguments[1]):
                if str(arguments[1]).endswith(".pdb") or str(arguments[1]).endswith(".gro") or str(arguments[1]).endswith(".xtc"):
                    pdb_structure = arguments[1]
                    print("Found .pdb structure: " + pdb_structure)
                else:
                    raise Exception("The file is not a .pdb file by file ending.")
            else:
                raise Exception("The .pdb-File could not be found.")

            # check for flags:
            flags = ["o"] + [key for key in parser.keys()]
            for i, x in enumerate(arguments[1:]):
                flag = str(x).replace("-", "")
                if flag in flags:
                    if flag == "o":
                        output_folder = arguments[i + 2]
                        if not os.path.isdir(os.path.dirname(output_folder)):
                            raise Exception("Could not find a path to Output Folder!\n")
                    else:
                        target = str(flag)
                        print("Found target: " + flag)
                        j = i + 2
                        tasks = str(arguments[j]).replace(" ", "")
                        print("Tasks: " + str(tasks))
                        path = str(depend["protocol_location"]) + "/" + str(parser[flag]["protocol_folder"])
                        translate_vol_prot(tasks, path)

                        if target == "1VPR" and os.path.isfile(parser[target]["Reference_struct"]):
                            reference_structure = parser[target]["Reference_struct"]
                            print("found Reference Structure: \n" + reference_structure)

            if target == "undefined":
                raise Exception("Could not find target! maybe uknonw?\n")
        else:
            raise Exception("No arguments or not the right ammount were provided. please look at the help text.\n")

    except Exception as err:
        print("Input Problem: \n" + str(err.args))
        print(help_text)
        exit()

    return pdb_structure, target, protocols, output_folder, reference_structure

def translate_vol_prot(tasks, path):
    """
    Args:
        tasks:
        path:
    """
    protocols = []
    for task in tasks:
        print(task)
        for y in os.listdir(path):
            if str(y).startswith("empty"):
                continue
            elif str(y).startswith(str(task)):
                tasks = tasks.replace(str(task), "")
                protocols.append(path + "/" + y)
                break

    if tasks is not "":
        raise Exception("Could not find all tasks!\n")

    return protocols
