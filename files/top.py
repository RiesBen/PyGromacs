import os
from PyGromacs.files import parser
from PyGromacs.files import gro

def clean_fields(list):
    """
    Args:
        list:
    """
    result = []
    for x in list:
        row = []
        field_str = ""
        for i, y in enumerate(x):
            if (not y == " "):
                field_str += y
                if (i == len(x) - 1):
                    row.append(field_str.strip("\n"))
                    field_str = ""
            else:
                if (len(field_str) > 0):
                    row.append(field_str.strip("\n"))
                    field_str = ""
        result.append(row)
    return result

def dictionize(list):
    """
    Args:
        list:
    """
    categories = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass', 'typeB', 'chargeB', 'massB']
    fields = clean_fields(list)
    atoms_dict = {}
    for atom_line in fields:
        atom_dict = {}
        for i, category in enumerate(categories):
            if (i < len(atom_line)):
                atom_dict.update({category: atom_line[i]})
            else:
                atom_dict.update({category: "undefined"})
        atoms_dict.update({str(atom_dict["nr"]): atom_dict})
    return atoms_dict

def include_line(name, path):
    """
    Args:
        name:
        path:
    """
    return ("; Include " + name + "\n"
            "#include \"" + path + ".itp\"\n")

def posres_line(name, condition):
    """
    Args:
        name:
        condition:
    """
    return ("; Include posres " + str(os.path.basename(name).split(".")[0]) + "\n"
            "#ifdef " + condition + "\n"
            "#include \"" + name + ".itp\"\n"
            "#endif \n\n")

def extract_gafft_types(segments):
    """
    Args:
        segments:
    """
    lines = []
    lines += "; These atom types were extracted form " + segments["name"] + "\n;in file: " + segments["path"] + "\n;with script complex_topol.py\n\n"
    lines += "[ atomtypes ]\n"
    lines += segments["atomtypes"]
    return lines

# just get all molecul props
def reduce_top_to_molecule_properties(segments):
    """
    Args:
        segments:
    """
    lines = []
    lines += "; These properties were extracted form " + segments["name"] + "\n;in file: " + segments["path"] + "\n;with script complex_topol.py\n\n"
    lines += ["[ moleculetype ]\n"] + segments["moleculetype"] + ["\n"]
    lines += ["[ atoms ]\n"] + segments["atoms"] + ["\n"]
    lines += ["[ bonds ]\n"] + segments["bonds"] + ["\n"]
    lines += ["[ pairs ]\n"] + segments["pairs"] + ["\n"]
    lines += ["[ angles ]\n"] + segments["angles"] + ["\n"]
    if ("dihedrals" in segments):
        lines += ["[ dihedrals ]\n"] + segments["dihedrals"] + ["\n"]
    elif ("dihedrals-propers" in segments):
        lines += ["[ dihedrals ] ;propers\n"] + segments["dihedrals-propers"] + ["\n"]
        lines += ["[ dihedrals ] ;impropers\n"] + segments["dihedrals-impropers"] + ["\n"]
    else:
        raise Exception("Couldnt find dihedrals!")
    return lines

def build_system_topology(molecules, force_field, water_model, output_file, atom_types_path=""):
    """
    Args:
        molecules:
        force_field:
        water_model:
        output_file:
        atom_types_path:
    """
    super_output_file_lines = ["this is a script generated top file ;)\n"]
    system_title = ""
    include_forcefield = "; include force_field\n#include \"" + force_field + ".ff/forcefield.itp\"\n\n"

    if (atom_types_path):
        include_gafft_file = "; include gafft Atom types\n#include \"" + atom_types_path + "\"\n\n"

    include_molecule = ""
    for molecule in molecules:
        include_molecule += "; include " + molecule + " \n#include \"topologies/" + os.path.basename(molecules[molecule]["top"]) + "\"\n\n"
        for x in molecules[molecule]["pos_res"]:
            include_molecule += posres_line("topologies/" + molecule + "_" + x[0], x[1])

    include_water_model = "; include water\n#include \"" + force_field + ".ff/" + water_model + ".itp\"\n\n"
    include_ion_model = "; include ions\n#include \"" + force_field + ".ff/ions.itp\"\n\n"

    # Part: [ system ]
    protein_top = parser.read_top(molecules["protein"]["top"])
    print(protein_top.keys())
    print(molecules["protein"]["top"])
    if (not "system" in protein_top):
        protein_top["system"] = "Protein in a box\n"
    system_lines = "[ system ] \n" + str(protein_top["system"]).strip("[").strip("]").strip("\n") + "\n"

    # Part: [ molecules ]
    molecule_lines = "[ molecules ] \n"
    add_later = ""
    if ("molecules" in protein_top):
        for x in protein_top["molecules"]:
            if (not str(x).startswith(";")):
                if (str(x).startswith("SOL") or str(x).startswith("CL") or str(x).startswith("NA")):
                    add_later += x
                else:
                    molecule_lines += x
    else:
        molecule_lines += "protein          1\n"
    if (len(molecules) > 1):
        for top in molecules[1:]:
            molecule_top = parser.read_top(molecules[top]["top"])
            for x in molecule_top["molecules"]:
                if (not str(x).startswith(";")):
                    molecule_lines += x
    molecule_lines += add_later

    super_output_file_lines.append(include_forcefield)
    if (atom_types_path):
        super_output_file_lines.append(include_gafft_file)
    super_output_file_lines.append(include_molecule)
    super_output_file_lines.append(include_water_model)
    super_output_file_lines.append(include_ion_model)
    super_output_file_lines.append(system_lines)
    super_output_file_lines.append(molecule_lines)

    file = open(output_file, "w")
    file.writelines(super_output_file_lines)
    file.close()

    return super_output_file_lines

def map_residue(input_file, mapping_file, output_file, residues, states="both"):
    """
    Args:
        input_file:
        mapping_file:
        output_file:
        residues:
        states:
    """
    protein_file = parser.read_top(input_file)
    mapping_file = parser.read_top(mapping_file)
    # Input Files:
    mapped_atoms = []
    if (not "atoms" in mapping_file):
        raise Exception("mapping file has no mapping section!")
    map_list = dictionize(mapping_file["atoms"])

    protein_atoms = protein_file["atoms"]
    for atom in protein_atoms:
        if ("HIS" in atom and not str(atom).startswith(";") and any([True if (" "+x in atom) else False for x in residues])):
            atom_fields = clean_fields([atom, ])[0]
            not_found = True
            for map_atom in map_list:
                not_found = True
                map = map_list[map_atom]
                if (map["atom"] == atom_fields[4] and map["residue"] in atom_fields[3] and map["typeB"] != ";"):
                    if(states == "both"):#atom num
                        print("Go to both states for: "+str(residues)+" ")
                        new_atom = "{0:>6}{1:>11}{2:>7}{3:>7}{4:>7}{5:>7}{6:>11}{7:>11}{8:>5}{9:>11}{10:>11} ;mapped\n".format(atom_fields[0], atom_fields[1], atom_fields[2], atom_fields[3], atom_fields[4], atom_fields[5], map["charge"], map["mass"], map["typeB"].replace(" ", ""), map["chargeB"], map["massB"])
                    elif (states == "B"):
                        print("Go to state B for: "+str(residues)+" ")
                        new_atom = "{0:>6}{1:>11}{2:>7}{3:>7}{4:>7}{5:>7}{6:>11}{7:>11}   ;mapped\n".format(atom_fields[0], atom_fields[1], atom_fields[2], atom_fields[3],  atom_fields[4], atom_fields[5], map["chargeB"], map["massB"])
                    elif (states == "A"):
                        print("Go to state A for: "+str(residues)+" ")
                        new_atom = "{0:>6}{1:>11}{2:>7}{3:>7}{4:>7}{5:>7}{6:>11}{7:>11}   ;mapped\n".format(atom_fields[0], atom_fields[1], atom_fields[2], atom_fields[3],  atom_fields[4], atom_fields[5], map["charge"], map["mass"])
                    else:
                        raise Exception("don't know this state! in topology_general_functions")
                    not_found = False
                    print(new_atom)
                    # ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass', 'typeB', 'chargeB', 'massB'
                    mapped_atoms.append(new_atom)
                    break

            if(not_found):
                print("no_map")
                new_atom = "{0:>6}{1:>11}{2:>7}{3:>7}{4:>7}{5:>7}{6:>11}{7:>11}   ;\n".format(atom_fields[0], atom_fields[1], atom_fields[2], atom_fields[3], atom_fields[4], atom_fields[5], atom_fields[6], atom_fields[7])
                print(new_atom)
                mapped_atoms.append(new_atom)
        else:
            mapped_atoms.append(atom)
    protein_file["atoms"] = mapped_atoms

    # Output
    output_file = open(output_file, "w")
    output_file.writelines(reduce_top_to_molecule_properties(protein_file))
    output_file.close()

def write_super_top_file(top_segments, output_top, coarse=False):
    """
    Args:
        top_segments:
        output_top:
        coarse:
    """
    print("Write lines to super_top")
    print(top_segments.keys())

    if(coarse):
        top_lines = "".join(top_segments["include"])
        top_lines += "[ system ] \n"+"".join(top_segments["system"])
        top_lines += "[ molecules ] \n"+"".join(top_segments["molecules"])
    else:
        top_lines = "".join(top_segments["header"])
        top_lines += "".join(top_segments["include_force_field"])
        top_lines += "".join(top_segments["include_protein"])
        #top_lines += "; include protein posres\n"+"".join(top_segments["include_protein_posres"])
        top_lines += "; offtopic \n"+"".join(top_segments["offtopic"])
        top_lines += "".join(top_segments["include_water"])
        top_lines += "".join(top_segments["include_ions"])
        top_lines += "[ system ] \n"+"".join(top_segments["system"])
        top_lines += "[ molecules ] \n"+"".join(top_segments["molecules"])
    #print(top_lines)
    new_top_file = open(output_top, "w")
    new_top_file.write(top_lines)
    new_top_file.close()

def insert_h3O(protein_top, insert_top, insert_res_top, output_path, insert_times=1):
    """
    Args:
        protein_top:
        insert_top:
        insert_res_top:
        output_path:
        insert_times:
    """
    print("parse: ")
    protein_dict = parser.read_top(protein_top)
    insert_dict = parser.read_top(insert_top)
    prefix = "topologies/H3O"
    h3o_out = output_path+"/"+prefix

    posres_H3O = posres_line(insert_res_top, "H3ORES") #generate a top include line
    molecules_seg = combine_molecules_segment(protein_dict["molecules"], insert_dict["molecules"], ["H3O", insert_times])   #count and sort molecules part

    print("build_top_file lines")
    print_top_content(protein_dict)
    sys_lines = []
    sys_lines+=["; :[ "]+protein_dict["header"]
    sys_lines+= protein_dict["include_force_field"]
    sys_lines+= protein_dict["include_protein"]
    sys_lines+= protein_dict["offtopic"]
    sys_lines+=["\n"+include_line(path=prefix, name="H3O")+"\n"+posres_H3O]
    sys_lines+= protein_dict["include_water"]
    sys_lines+= protein_dict["include_ions"]
    sys_lines+= ["\n[ system ]\n"]+[protein_dict["system"][0].strip() +" and "+insert_dict["system"][0]]
    sys_lines+=[molecules_seg]

    categories_h3o = ["moleculetype", "atoms", "constraints", "virtual_sites3", "exclusions", "settles"]
    h3O_lines=attribute_segments_to_string(categories_h3o, insert_dict)

    print("write output")
    output_sys = open(protein_top, "w")
    output_sys.writelines(sys_lines)
    output_sys.close()

    output_h3o = open(h3o_out+".itp", "w")
    output_h3o.writelines(h3O_lines)
    output_h3o.close()

def combine_molecules_segment(molecules1, molecules2, adjust_mol):
    """
    Args:
        molecules1:
        molecules2:
        adjust_mol:
    """
    output_str = "\n[ molecules ]\n"
    mol_list=[]
    for x in molecules1:
        element=[]
        for y in x.split(" "):
            if(y.strip() != ""):
                element.append(y.strip())
        mol_list.append(element)

    for x in molecules2:
        element=[]
        for y in x.split(" "):
            if(y.strip() != ""):
                element.append(y.strip())
        mol_list.append(element)

    #remove duplicate:
    mol_clean = []
    for pos,mol in enumerate(mol_list):
        doubles = [i for i, x in enumerate(mol_list) if mol[0] == x[0]]
        if( sum(doubles) > pos):
            mol_name = mol[0]
            mol_count = 0
            for x in doubles:
                mol_count += int(mol_list[x][1])
            offset = 0
            for x in doubles[1:]:
                mol_list.remove(mol_list[x-offset])
                offset += 1
            mol_clean.append([mol_name, str(mol_count)])
        elif(mol[0] in adjust_mol):
            mol[1] = adjust_mol[1]
            mol_clean.append(mol)
        else:
            mol_clean.append(mol)

    #sort mols
    mols=len(mol_clean)
    mol_pos = []
    counter = 3
    for x in mol_clean:
        if(";" in x):
             mol_pos.append([1, x])
        elif("NA" in x):
             mol_pos.append([mols-2, x])
        elif("CL" in x):
            mol_pos.append([mols-1, x])
        elif("SOL" in x):
            x[1] = str(int(x[1])-adjust_mol[1])
            mol_pos.append([mols, x])
        elif("protein" in x):
            mol_pos.append([2, x])
        else:
            mol_pos.append([counter, x])
            counter += 1
    print(mol_pos)
    for count in range(mols+1):
        for pos in mol_pos:
            if(pos[0] == count):
                if(len(pos[1])==3):
                    output_str += "{} {:<11}{:>10}\n".format(pos[1][0], pos[1][1], pos[1][2])
                else:
                    output_str += "{:<13}{:>10}\n".format(pos[1][0], pos[1][1])
    return output_str

def attribute_segments_to_string(categories, dictionary):
    """
    Args:
        categories:
        dictionary:
    """
    lines=[]
    for x in categories:
        lines += ["\n[ " + x + " ]\n"] + dictionary[x] + ["\n"]
    return lines

def print_top_content(top_segments):
    """
    Args:
        top_segments:
    """
    for x in top_segments:
        if(top_segments[x] is (str or int or float)):
            print(x + "\n" +str(top_segments[x]) + "\n")
        elif(top_segments[x] is list or dict):
            print(x+"\n"+"".join(top_segments[x])+"\n")
        else:
            print("skipped: "+str(x)+"\n"+str(top_segments[x])+" \n"+str(type(top_segments[x])))

def update_system_mol_nums(gro_file, top_file):
    """
    Args:
        gro_file:
        top_file:
    """
    print("read in gro:")
    protein = parser.read_gro(gro_file)
    print("read in top:")
    top_segments = parser.read_top(top_file, coarse=True)

    print("original molecules" + str(top_segments["molecules"]))
    new_sys = []
    for mols in top_segments["molecules"]:
        if (";" in mols):
            continue
        if ("protein" in mols or "Protein" in mols):
            new_sys.append(mols)
            continue
        else:
            name = mols.split("  ")[0].strip()
            residue_count, residue_nums, atom_nums = gro.count_residue(gro.clean_fields(protein["coordinates"]), residue_name=name)
            if(residue_count == 0):
                continue
            else:
                new_sys.append("{:<13}{:>10}\n".format(name, str(residue_count)))

    print("counted: \n" + str(new_sys))
    top_segments["molecules"] = new_sys
    print("write updated top file")
    print_top_content(top_segments)
    write_super_top_file(top_segments, top_file, coarse=True)
