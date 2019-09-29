"""
WORKFLOWMODULE: workflow_parser - module
Description:
    this parser, checks the input to a workflow.
Input:
Output:
Dependencies outside of the code:
Author: Benjamin Schroeder
"""
import os
import generalutilities.function_libs.utils.bash_commands as bash

#______________________________________________________________
# Workflow_parsing_commandline
#______________________________________________________________
def do_workflow_parse(arguments, help_text, md_time="", cphmd_barrier=""):
    """
    Args:
        arguments:
        help_text:
        md_time:
        cphmd_barrier:
    """
    output_folder = "undefined"
    comand_line_options = {}
    try:
        try:
            for i,x in  enumerate(arguments):
                if(x == "-o" and len(arguments) > i+1 and not str(arguments[i+1]).startswith("-")):
                    output_folder = arguments[i+1]
                elif(str(x).startswith("-") and len(arguments) > i + 1 and not str(arguments[i + 1]).startswith("-")):
                    if(x.replace("-", "") in comand_line_options):
                        raise Exception("command_line option is given twice")
                    else:
                        comand_line_options.update({x.replace("-", ""): arguments[i + 1]})
                elif ("-f" in arguments):
                    comand_line_options.update({"Ff": ""})  # flages
        except Exception as err:
            print(bash.increment_error_level("IO-Problem", err))
    except Exception as err:
        print(bash.increment_error_level("IO-Problem", err))
        print("\n\n\n"+help_text)
        exit()
    return output_folder, comand_line_options

#______________________________________________________________
# GMX FIles
#______________________________________________________________
def read_top(path, coarse=False):
    """
    Args:
        path:
        coarse:
    """
    print("open path: "+str(path))
    topol_file = open(path, "r")
    topol_section_lines = {}

    topol_section_lines.update({"path": path})
    topol_section_lines.update({"name": os.path.basename(path)})
    topol_section_lines.update({"offtopic": []})

    print("start reading file:")
    if(coarse):
        segment = "include"
        topol_section_lines = {segment: []}
        for x in topol_file.readlines():
            if ("[" and "]" in x):
                segment = x.replace(";", "-").replace("[", "").replace("]", "").strip().replace(" ", "")
                topol_section_lines.update({segment: []})
            elif (("include" in x or "Include" in x) and not "#" in x and not "posres" in x):
                topol_section_lines["include"].append(x)
            else:
                topol_section_lines[segment].append(x)
    else:
        segment = "header"
        topol_section_lines.update({segment: []})
        known_include= ["force_field", "protein", "water", "ions"]
        for x in topol_file.readlines():
            if("[" and "]" in x):
                segment = x.replace(";", "-").replace("[","").replace("]","").strip().replace(" ","")
                topol_section_lines.update({segment:[]})
            elif(("include" in x or "Include" in x) and not "#" in x and not "posres" in x):
                print(x)
                if(any([True for y in known_include if (y in x)])):
                    if(len(x.split(" ")) <= 3):
                        segment = (x.split(" ")[1]+"_"+x.split(" ")[2]).strip().replace(";", "")
                    elif(len(x.split(" ")) > 3):
                        segment = (x.split(" ")[1]+"_"+"".join(x.split(" ")[2:])).strip().replace(";", "")
                    else:
                        segment = (x.split(" ")[1]).strip().replace(";", "")
                else:
                    segment="offtopic"
                topol_section_lines.update({segment:[x]})
            else:
                topol_section_lines[segment].append(x)
    return topol_section_lines

def check_gro_term(line):
    """
    Args:
        line:
    """
    for x in line.split(" "):
        if(len(x) > 0 and any(y.isalpha() for y in x )):
            return False
    return True

def read_gro(path):
    """
    Args:
        path:
    """
    print("Parsing_file: "+path)
    try:
        gro_file = open(path, "r")
    except Exception as err:
        raise Exception("could not find File: "+path+"\n "+str(err.args))

    gro_lines = gro_file.readlines()

    gro_segments = {}
    gro_segments.update({"path": path})
    gro_segments.update({"name": os.path.basename(path)})
    estimate_atoms = False
    pointer= 0
    if(len(gro_lines[pointer].split("   "))<3):
        gro_segments.update({"header": gro_lines[pointer]})
        print("FoundHeader: "+gro_lines[pointer])
        pointer += 1
    else:
        gro_segments.update({"header": "None"})
    if( gro_lines[pointer].strip().isalnum() ):
        if(int(gro_lines[pointer].strip())==len(gro_lines)-3):
            gro_segments.update({"atoms": gro_lines[pointer]})
            print("Found atom num: "+gro_lines[pointer])
            pointer += 1
        else:
            raise Exception("Not the same ammount of atoms!")
    else:
        estimate_atoms = True
    term_pointer = len(gro_lines)-1

    if(check_gro_term(gro_lines[term_pointer])):
        gro_segments.update({"coordinates": gro_lines[pointer:term_pointer]})
        gro_segments.update({"term": gro_lines[term_pointer]})
    else:
        gro_segments.update({"coordinates": gro_lines[pointer:]})
        gro_segments.update({"term": "\n"})

    if(estimate_atoms):
        print("no AtomNumber in  line "+str(pointer)+" - Estimate:"+str(len(gro_segments["coordinates"])))
        gro_segments.update({"atoms": str(len(gro_segments["coordinates"]))})
    else:
        if(int(gro_segments["atoms"])!=len(gro_segments["coordinates"])):
            print("atoms_segment: "+gro_segments["atoms"]+" found_atoms: "+str(len(gro_segments["coordinates"])))
            raise Exception("Wrong amount of coordinate-lines and total atom number")
    print("Done Parsing!")
    return gro_segments

def read_ndx_file(path):
    """
    Args:
        path:
    """
    file = open(path, "r")
    file_lines = file.readlines()
    file_dict = {}
    key = "undefined"
    #write dict
    for x in file_lines:
        if(x.startswith("[") and x.strip().endswith("]")):
            key = x.replace("[", "").replace("]", "").strip()
            file_dict.update({key: []})
        else:
                file_dict[key] += x.strip("\n").split(" ")
    #clean
    for key in file_dict:
        for x in file_dict[key]:
            if(x == "" or x ==" " or x == "\n"):
                file_dict[key].remove(x)
    return file_dict
