import generalutilities.function_libs.gromacs.parser as parser
import generalutilities.function_libs.utils.bash_commands as gf
import random
import os
import re


def write_gro_file(gro_dict, path):
    """
    Args:
        gro_dict:
        path:
    """
    output_file = open(path, "w")
    print("Writing: " + gro_dict["header"] + "contains atoms:"+gro_dict["atoms"])
    output_file.write(gro_dict["header"])
    output_file.write(gro_dict["atoms"].strip("\n") + "\n")
    output_file.writelines(gro_dict["coordinates"])
    output_file.write(gro_dict["term"])
    output_file.close()
    return path

def clean_up_gro_file(gro_file):
    """
    Args:
        gro_file:
    """
    if (not "header" in gro_file):
        gro_file.update({"header": "A formidable header!"})
    if ("atoms" in gro_file):
        gro_file["atoms"] = "{: >5}\n".format(str(gro_file["atoms"].strip().strip("\n")))
    if ("coordinates" in gro_file):
        gro_fields = clean_fields(gro_file["coordinates"])
        gro_file["coordinates"] = make_gro_atom_lines(gro_fields)
    if ("term" in gro_file):
        if (len(gro_file["term"].split("   ")) == 4):
            y = []
            for x in gro_file["term"].split("  "):
                x = x.strip().strip("\n")
                y.append(x)
            gro_file["term"] = "{: >10.5}{: >10.5}{: >10.5}\n".format(y[1], y[2], y[3])
        elif (len(gro_file["term"].split("   ")) == 10):
            y = []
            for x in gro_file["term"].split("  "):
                x = x.strip().strip("\n")
                y.append(x)
            gro_file["term"] = "{: >10.5}{: >10.5}{: >10.5}{: >10.5}{: >10.5}{: >10.5}{: >10.5}{: >10.5}{: >10.5}\n".format(y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9])
        else:
            gro_file["term"] = "{: >10}{: >10}{: >10}\n".format("0.00000", "0.00000", "0.00000")
    return gro_file

def clean_fields(list):
    """
    Args:
        list:
    """
    columns=[5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8]
    result = []
    for x in list:
        row = []
        field_str = ""
        column_ind = 0
        column=columns[column_ind]
        for i, y in enumerate(x):
            field_str += y
            if i == (column -1):
                column_ind += 1
                column += columns[column_ind]
                row.append(field_str.replace(" ", "").strip("\n"))
                field_str=""
        result.append(row)#
    return result

def make_gro_atom_lines(gro_ar, renumber=False):
    """
    Args:
        gro_ar:
        renumber:
    """
    gro_lines = []

    res_num = gro_ar[0][0]
    new_res_num = int(gro_ar[0][0])
    new_atom_num = 0
    for x in gro_ar:
        if(renumber):
            if(int(res_num) != int(x[0])):
                res_num = x[0]
                new_res_num += 1
            x[0] = str(new_res_num)

            if (int(new_atom_num) == 10000-1):
                new_atom_num = 0
            new_atom_num +=1
            x[3] = str(new_atom_num)

        if (len(x) == 7):
            gro_lines.append("{: >5}{: <5}{: >5}{: >5}{: >8.3f}{: >8.3f}{: >8.3f}\n".format(str(x[0]), x[1], x[2], str(x[3]), float(x[4]), float(x[5]), float(x[6])))
        elif (len(x) == 10):
            gro_lines.append("{: >5}{: <5}{: >5}{: >5}{: >8.3f}{: >8.3f}{: >8.3f}{: >8.4f}{: >8.4f}{: >8.4f}\n".format(str(x[0]), x[1], x[2], str(x[3]), float(x[4]), float(x[5]), float(x[6]), float(x[7]), float(x[8]), float(x[9])))
        else:
            raise Exception("No known line length! \n" + str(x))
    return gro_lines

def map_ligand_atom_names_to_new_ligand_pos(ligands_gro_path, ligand_template_gro_path, res_name="LIG"):
    """maps all atom names of a template which has the same atom order to all
    atoms of another molecule :param ligands_gro_path: :param
    ligand_template_gro_path: :param res_name: :return:

    Args:
        ligands_gro_path:
        ligand_template_gro_path:
        res_name:
    """

    #INPUT:
    print("start Mapping Process:\n Cleaning old_pos")
    ligand_positions = parser.read_gro(ligands_gro_path)
    clean_new_pos = clean_fields(ligand_positions["coordinates"])  # new complexfile
    # print(clean_ligand)
    print("\n Cleaning new_pos")
    ligand_template = parser.read_gro(ligand_template_gro_path)
    clean_ligand = clean_fields(ligand_template["coordinates"])

    # print(clean_new_pos)

    #Initial Analyse
    #count_residuest to map.
    residue_count, resn_numbers, atoms_per_res = count_residue(clean_new_pos, res_name)
    print("found new residue_positions: "+str(residue_count))
    print("each residue and its atom count: \n"+str(atoms_per_res))
    print("start Mapping")

    #filter lines for residue and seperate them:
    lines_to_map = {}
    for residue in resn_numbers:
        lines_to_map.update({residue:[line for line in clean_new_pos if (line[0] == residue)]})

    print("total num of residues: "+str(len(resn_numbers)))
    for residue in resn_numbers:
        print("res: "+str(residue)+"\t"+str(len(lines_to_map[residue]))+"\n")
    #DO
    #vars:
    mapped_residues_lines = []

    # map lines for each residue seperate:
    for residue in resn_numbers:
        tmp_ligand = clean_ligand
        map_list = lines_to_map[residue]
        atom_count = 0
        mapped_sub = []
        print("mapping residue: "+str(residue))
        #Mapping process
        for y in tmp_ligand:
            for x in map_list:
                if(re.search(x[2], y[2])):
                    mapped_line = [x[0]] + y[1:4] + x[4:]
                    mapped_sub.append(mapped_line)
                    map_list.remove(x)
                    atom_count += 1
                    break
        if(len(map_list)==0):
            mapped_residues_lines += mapped_sub
        else:
            print("skipping this molecule! it has more atoms!\n"+str(map_list))

        print("len of mapped_list: "+str(len(mapped_residues_lines)) +" added lines: "+str(atom_count))

    #Result
    gro_lines = make_gro_atom_lines(mapped_residues_lines, renumber=True)
    ligand_positions["atoms"] = str(len(gro_lines))
    ligand_positions["coordinates"] = gro_lines
    mapped_path = os.path.splitext(ligands_gro_path)[0]+"_mapped_lines.gro"
    write_gro_file(ligand_positions, mapped_path)
    return mapped_path

def append_structures(protein_gro, ligand_gro):
    """
    Args:
        protein_gro:
        ligand_gro:
    """
    ligand_gro = clean_up_gro_file(ligand_gro)
    protein_gro = clean_up_gro_file(protein_gro)
    lines = [protein_gro["header"]] + [str(int(protein_gro["coordinates"]) + int(ligand_gro["coordinates"])) + "\n"] + protein_gro["coordinates"] + ligand_gro["coordinates"] + [protein_gro["term"]]
    return lines

def append_structure_files(protein_gro_path, ligand_gro_path, output_path=False):
    """
    Args:
        protein_gro_path:
        ligand_gro_path:
        output_path:
    """
    gro_file_prot = parser.read_gro(protein_gro_path)
    gro_file_lig = parser.read_gro(ligand_gro_path)

    ligand_gro = clean_up_gro_file(gro_file_lig)
    protein_gro = clean_up_gro_file(gro_file_prot)
    combined_coordinates = (protein_gro["coordinates"] + ligand_gro["coordinates"])
    print(combined_coordinates)
    compbined_gro = {"header": protein_gro["header"], "atoms": str(int(len(protein_gro["coordinates"])) + int(len(ligand_gro["coordinates"]))) + "\n", "coordinates": combined_coordinates , "term": protein_gro["term"]}

    if(not output_path):
        output_path= os.path.splitext(protein_gro_path)[0]+"_concat.gro"

    write_gro_file(compbined_gro,output_path)
    return output_path

def get_atom_shift(atom_O_line, h3O_O_line):
    """
    Args:
        atom_O_line:
        h3O_O_line:
    """
    if(len(atom_O_line) == 7 or len(atom_O_line) == 10):
        x_shift = float(atom_O_line[4]) - float(h3O_O_line[4])
        y_shift = float(atom_O_line[5]) - float(h3O_O_line[5])
        z_shift = float(atom_O_line[6]) - float(h3O_O_line[6])
    else:
        raise Exception("Could not calculate atom shift!")
    #   print(str(x_shift)+"  "+str(y_shift)+"   "+str(z_shift))
    return x_shift, y_shift, z_shift

def move_h3o_on_water(water_line, insert_coordinates, molecules):
    """
    Args:
        water_line:
        insert_coordinates:
        molecules:
    """
    coordinates_H3O = []
    x_shift, y_shift, z_shift = get_atom_shift(water_line, insert_coordinates[0])
    for atom in insert_coordinates:
        atom[0] = str(molecules)
        if (atom[1] == "OW"):  # map O onto water O
            if(len(water_line) == 7 or len(water_line) == 10):
                atom[4] = water_line[4]
                atom[5] = water_line[5]
                atom[6] = water_line[6]
                if(len(water_line) == 10):
                    atom = atom[:7] + water_line[7:]
            else:
                raise Exception("line unrecogniseable: \n"+str(water_line))
        else:  # adjust Hs to O atom
            if(len(water_line) == 7 or len(water_line) == 10):
                atom[4] = str(round(float(atom[4]) + x_shift, 3))
                atom[5] = str(round(float(atom[5]) + y_shift, 3))
                atom[6] = str(round(float(atom[6]) + z_shift, 3))
                if(len(water_line) == 10):
                    atom = atom[:7] + water_line[7:]
            else:
                raise Exception("line unrecogniseable: \n"+str(water_line))
        coordinates_H3O.append(atom)
    return coordinates_H3O

def insert_h3O(selection_file, protein_file, insert_mol, output, selection_name ="waterSel", insert_times=1):
    #input
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
    protein_coordinates = clean_fields(protein_gro["coordinates"])
    insert_gro = parser.read_gro(insert_mol)
    insert_coordinates = clean_fields(insert_gro["coordinates"])

    #get random SOL residues from selection
    #getting precise selection keay
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
        x = random.randint(0, positions-1)
        print("for pos "+str(x)+"\t\t"+str(selection_dict[selection_key][x]))
        replace_positions.append(x)
        replace_residues.append(selection_dict[selection_key][x])
    print("Replace position: \t" + str(replace_positions) + " from total pos: " + str(positions))
    print("\nReplace Residue: \t" + str(replace_residues))

    #get residues, to be replaced
    replace_molecule = []
    for i,line_prot in enumerate(protein_coordinates):
        for x in replace_residues:
            if (i == int(x)):
                print("found: "+str(x)+"\t"+str(line_prot))
                replace_residues.remove(x)
                replace_molecule.append(line_prot[0])

    print("\nreplace mol: " + str(replace_molecule))

    #insert molecules - split coordinates in parts
    print("Map Coordinates on H2Os")
    new_coordinates_H3O = []
    new_coordinates_NA = []
    new_coordinates_CL = []
    new_coordinates_SOL = []
    new_coordinates_rest = []
    molecule = 1
    for line_prot in protein_coordinates:
        if (str(line_prot[0]) in replace_molecule):
            if ("OW" in line_prot[2]):
                print("H3O replace: " + str(molecule))
                print("Found: "+str(line_prot))
                new_coordinates_H3O += move_h3o_on_water(line_prot, insert_coordinates, molecule)
                molecule += 1
                continue
            else:
                print("ignore: "+str(line_prot))
                continue
        elif ("SOL" in line_prot[1]):
            new_coordinates_SOL.append(line_prot)
        elif ("NA" in line_prot[1]):
            new_coordinates_NA.append(line_prot)
        elif ("CL" in line_prot[1]):
            new_coordinates_CL.append(line_prot)
        else:
            new_coordinates_rest.append(line_prot)

    new_coordinates = new_coordinates_rest + new_coordinates_H3O + new_coordinates_SOL + new_coordinates_NA + new_coordinates_CL
    print(" Length of new coordinates: "+str(len(new_coordinates)))
    protein_gro["atoms"] = str(len(new_coordinates))
    protein_gro["coordinates"] = make_gro_atom_lines(new_coordinates, renumber=True)
    write_gro_file(protein_gro, output)

def get_residue_atoms(residue, residue_name, gro_coord, atomtype):    #mnot cool
    """
    Args:
        residue:
        residue_name:
        gro_coord:
        atomtype:
    """
    atoms = []
    if(not residue):
        for x in gro_coord:
            if (residue_name == x[1] and x[2] in atomtype):
                print(x)
                atoms.append(x[3])

    else:
        for x in gro_coord:
            if (residue == x[0] and residue_name == x[1] and x[2] in atomtype):
                print(x)
                atoms.append(x[3])
    return atoms

def update_velocities(gro, vel_gro):
    """
    Args:
        gro:
        vel_gro:
    """
    gro_file = parser.read_gro(gro)
    vel_file = parser.read_gro(vel_gro)

    gro_coord = clean_fields(gro_file["coordinates"])
    vel_coord = clean_fields(vel_file["coordinates"])

    write = True
    new_coord =[]
    for i,x in enumerate(gro_coord):
        if (len(x) >= 8):
            print("Found velocities in gro - file! \n"+str(x)+"\n Abort transfer")
            write = False
            break
        else:
            if(len(vel_coord[i]) == 9 ):
                new_coord.append(x+vel_coord[i][6:])
            elif(len(vel_coord[i]) == 8 ):
                new_coord.append(x+vel_coord[i][5:])

    if(write):
        gro_file["coordinates"] = make_gro_atom_lines(new_coord)
        write_gro_file(gro_file, gro)

def H3O_to_H2O(protein_file, output, replace_times=1):

    #input
    """
    Args:
        protein_file:
        output:
        replace_times:
    """
    protein_gro = parser.read_gro(protein_file)  #
    protein_coordinates = clean_fields(protein_gro["coordinates"])

    #insert molecules - split coordinates in parts
    new_coordinates_H3O = []
    new_coordinates_NA = []
    new_coordinates_CL = []
    new_coordinates_SOL = []
    new_coordinates_rest = []
    molecule = 1
    if(replace_times == "all"):
        replace_times_atoms = len(protein_coordinates)
    else:
        replace_times_atoms = 6*replace_times

    print("going to transform H3O -> H2O molecules: \tmols: "+str(replace_times)+"\t atoms: "+str(replace_times_atoms))
    for line_prot in protein_coordinates:
        if ("H3O" in line_prot[1] and replace_times_atoms > 0):
            if (not "H31" in line_prot[2] and not "H32" in line_prot[2] and not "H33" in line_prot[2]):
                print("H3O replace: " + str(molecule))
                print("Found: "+str(line_prot))
                line_prot = [line_prot[0]]+["SOL"]+line_prot[2:]
                new_coordinates_SOL.append(line_prot)

                molecule += 1
                replace_times_atoms -= 1
                print("remove  atom "+str(replace_times_atoms)+" of "+str(replace_times))
                continue
            else:
                print("ignore: "+str(line_prot))
                replace_times_atoms -= 1
                continue
        elif ( "H3O" in line_prot[1]):
            new_coordinates_H3O.append(line_prot)
        elif ("SOL" in line_prot[1]):
            new_coordinates_SOL.append(line_prot)
        elif ("NA" in line_prot[1]):
            new_coordinates_NA.append(line_prot)
        elif ("CL" in line_prot[1]):
            new_coordinates_CL.append(line_prot)
        else:
            new_coordinates_rest.append(line_prot)

    new_coordinates = new_coordinates_rest + new_coordinates_H3O + new_coordinates_SOL + new_coordinates_NA + new_coordinates_CL
    print(" Length of new coordinates: "+str(len(new_coordinates)))
    protein_gro["atoms"] = str(len(new_coordinates))
    protein_gro["coordinates"] = make_gro_atom_lines(new_coordinates)
    write_gro_file(protein_gro, output)

def count_residue(coordinates, residue_name="LIG"):
    """
    Args:
        coordinates:
        residue_name:
    """
    print("Counting residue "+str(residue_name)+" in coordinates")
    residue_nums = []
    atom_nums = {}
    atom_num = -1
    residue_num = -1
    for line in coordinates:
        if(line[1] == residue_name and not line[0] in residue_nums):
            residue_num=line[0]
            residue_nums.append(residue_num)
            atom_nums.update({residue_num : 1})
        elif(line[1] == residue_name and line[0] in residue_nums):
            atom_nums[residue_num] += 1

    residue_count =len(residue_nums)
    if(residue_num == -1 ):
        print("couldn't find any residue: \t"+residue_name+"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n STILL CONTINUING ;)")
    #else:
        #print("Residues: \t"+str(residue_count))
        #print("Residue_numbers: \t"+str(residue_nums))
        #print(atom_nums)
    return residue_count, residue_nums, atom_nums

def filter_file(input_gro_file, output_prefix=False, residue_name=False, atom_name=False, residue_number=False, atom_number=False):
    """
    Args:
        input_gro_file:
        output_prefix:
        residue_name:
        atom_name:
        residue_number:
        atom_number:
    """
    output_name = output_prefix if (output_prefix) else os.path.splitext(input_gro_file)[0]
    gro_file = parser.read_gro(input_gro_file)
    filtered_lines = []
    removed_lines = []

    #get criteria for filtering:
    filter_dict = {}
    if(residue_number):
        filter_dict.update({0: residue_number})
    if(residue_name):
        filter_dict.update({1: residue_name})
    if(atom_name):
        filter_dict.update({2: atom_name})
    if (atom_number):
        filter_dict.update({3: atom_number})

    #filter:
    lines = clean_fields(gro_file["coordinates"])
    for line in lines:
        if(all([True if (line[key]==filter_dict[key]) else False for key in filter_dict])):
            removed_lines.append(line)
        else:
            filtered_lines.append(line)

    filterd_name=output_name+"_filtered.gro"
    gro_file["atoms"] = str(len(filtered_lines))
    gro_file["coordinates"] = make_gro_atom_lines(filtered_lines)
    write_gro_file(gro_file, filterd_name)

    removed_name=output_name+"_removed.gro"
    gro_file["atoms"] = str(len(removed_lines))
    gro_file["coordinates"] = make_gro_atom_lines(removed_lines)
    write_gro_file(gro_file, removed_name)

    return  filterd_name, removed_name

def protonate_iterativley(gro_file_path, res_name="LIG", protonation_ph=7.4):
    """
    Args:
        gro_file_path:
        res_name:
        protonation_ph:
    """
    gro_file = parser.read_gro(gro_file_path)
    lines = clean_fields(gro_file["coordinates"])
    resn_num, residue_nums, atom_nums = count_residue(lines)
    print("protonation found residues : "+str(resn_num))
    tmp_res=-1
    tmp_res_lines = []
    protonated_list = []
    for x in lines:
        if(res_name == x[1]):
            if(tmp_res!=int(x[0])):
                if(len(tmp_res_lines)):
                    protonated_list += protonate_lines(tmp_res_lines, gro_file_path, protonation_ph)
                resn_num -= 1
                tmp_res = int(x[0])
                tmp_res_lines = [x]
            else:
                tmp_res_lines.append(x)
    if(len(tmp_res_lines)):
        protonated_list += protonate_lines(tmp_res_lines,gro_file_path, protonation_ph)
    gro_file["atoms"]= str(len(protonated_list))
    gro_file ["coordinates"] = make_gro_atom_lines(protonated_list)
    gro_file_protonated_path = os.path.splitext(gro_file_path)[0]+"_protonated.gro"
    write_gro_file(gro_file, gro_file_protonated_path)
    return gro_file_protonated_path#

def protonate_lines(coordinate_lines, output_path, protonation_ph=7.4):
    """
    Args:
        coordinate_lines:
        output_path:
        protonation_ph:
    """
    lines = make_gro_atom_lines(coordinate_lines)
    tmp_gro = {"header":"tmp_file\n", "atoms": str(len(coordinate_lines)), "coordinates": lines, "term": "\n"}
    tmp_path =  os.path.dirname(output_path)+"/tmp.gro"
    tmp_path = write_gro_file(tmp_gro, tmp_path)
    #THINK ABOUT HOW TO CHECK OPENBABEL!
    protonation_command = "module load openbabel/2.4.1\n module load gcc/5.3.0\n obabel -igro "+tmp_path+" -ogro -O "+tmp_path+" -p "+str(protonation_ph)+"\n"
    gf.execute_bash(protonation_command)

    tmp_gro_file = parser.read_gro(tmp_path)
    protonated_lines = clean_fields(tmp_gro_file["coordinates"])
    return protonated_lines
