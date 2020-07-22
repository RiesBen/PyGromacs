"""
FUNCTIONLIB:            generalFunctions_analysis
Description:
   This lib contains funtions for analysisin a gromacs simulation.

Author: Benjamin Schroeder
"""

#imports
import os
from PyGromacs.utils import bash

#analysis:
def rmsd_fit(trajectory, reference_structure, fit_to, index_file=False, gmx="gmx",  output=False):
    """
    Args:
        trajectory:
        reference_structure:
        fit_to:
        index_file:
        gmx:
        output:
    """
    if (not output):
        output = os.path.dirname(trajectory) + "/selection_" + os.path.basename(trajectory)

    if(not reference_structure):
        reference_structure=trajectory

    if (index_file):
        command = gmx + " select -f " + trajectory + " -s " + reference_structure + " -n " + index_file + " -on " + output + " -select "+fit_to+" \n"
    else:
        command = "cd "+os.path.dirname(output)+"\n "+ gmx + " select -f " + trajectory + " -s " + reference_structure + " -select "+fit_to+" -on " + output + " \n"

    bash.execute(command)
    return output

def rmsd_residual(trajectory, reference_structure=False, gmx="gmx",  output=False, group_selection="1"):
    """
    Args:
        trajectory:
        reference_structure:
        gmx:
        output:
        group_selection:
    """

    if (not output):
        output = os.path.dirname(trajectory) + "/rmsd_residual_" + os.path.basename(trajectory)

    if (reference_structure):
        reference_structure = trajectory

    command = "cd "+os.path.dirname(output)+"\necho \""+group_selection+"\" | "+ gmx + " rmsf -f " + trajectory + " -s " + reference_structure + " -od " + output + " -res\n"
    bash.execute(command)
    return output

def rmsf_residual(trajectory, reference_structure=False, gmx="gmx",  output=False, group_selection="1"):
    """
    Args:
        trajectory:
        reference_structure:
        gmx:
        output:
        group_selection:
    """

    if (not output):
        output = os.path.dirname(trajectory + "/rmsf_residual_" + os.path.basename(trajectory))

    if (reference_structure):
        reference_structure = trajectory

    command = "echo \""+group_selection+"\" | "+ gmx + " rmsf -f " + trajectory + " -s " + reference_structure + " -o " + output + " -res\n"
    bash.execute(command)
    return output

def hbond(structure_file, tpr_file, group_1, group_2, output, index_file=False, gmx=" gmx ", num=False, dist=False, life=False, hbn=False, hbm=False, ac=False, hx=False):
    """
    Args:
        structure_file:
        tpr_file:
        group_1:
        group_2:
        output:
        index_file:
        gmx:
        num:
        dist:
        life:
        hbn:
        hbm:
        ac:
        hx:
    """

    selection = "echo \""+group_1+" "+group_2+"\" "
    executes = []
    ndx=""
    if(num):
        executes.append("num")
    if(dist):
        executes.append("dist")
    if(life):
        executes.append("life")
    if(hbn):
        executes.append("hbn")
    if(hbm):
        executes.append("hbm")
    if(ac):
        executes.append("ac")
    if(hx):
        executes.append("hx")
    if(index_file):
        ndx=" -n "+str(index_file)+" "

    if(os.path.isfile(output)):
        output_dir = os.path.dirname(output)
    else:
        output_dir = output

    for ex in executes:
        if(ex=="hbm"):
            hbond_command = "cd  " + selection + " \n" + output_dir + " | " + gmx + " hbond -f " + structure_file + " -s " + tpr_file + " " + ndx + " -" + ex + " " + output + "_" + ex + "\n"  # shows numbers of H-bonds over time
            hbond_command += "cd  "+output_dir+"\n gmx xpm2ps -f "+output+"_"+ex+".xpm -o "+output+"_"+"ex \n"
        else:
            hbond_command = "cd  " + output_dir + " \n" + selection + " | " + gmx + " hbond -f " + structure_file + " -s " + tpr_file + " " + ndx + " -" + ex + " " + output + "_" + ex + "\n"  # shows numbers of H-bonds over time

        bash.execute(hbond_command)

def dist(structure_file, index_file, selection, output, oav=True, oallstat=True, oh=True, gmx="gmx"):
    """
    Args:
        structure_file:
        index_file:
        selection:
        output:
        oav:
        oallstat:
        oh:
        gmx:
    """

    if (not os.path.isdir(output)):
        output_file = output
        output_dir = os.path.dirname(output)
    else:
        output_file = output + "/dist_" + selection
        output_dir = output

    executes = []
    if (oav):
        executes.append("oav ")
    if (oallstat):
        executes.append("oallstat ")
    if (oh):
        executes.append("oh ")

    for ex in executes:
        dist_command = "cd  " + output_dir + " \n" + gmx + " distance -f " + structure_file + " -n "+index_file+" -select " + selection + " -" + ex + " " + output_file + "_" + ex + "\n"  # shows numbers of H-bonds over time
        bash.execute(dist_command)

    return output

def rdf(structure_file, index_file, select_1, select_2, output, gmx="gmx"):
    """
    Args:
        structure_file:
        index_file:
        select_1:
        select_2:
        output:
        gmx:
    """

    if (not os.path.isdir(output)):
        output_file =output
        output_dir = os.path.dirname(output)
    else:
        output_file = output+"/rdf_"+select_1+"_"+select_2
        output_dir = output

    selection = "echo \""+select_1+"\" "
    rdf_command = "cd  " + output_dir + " \n" + selection + " | " + gmx + " rdf -f " + structure_file + " -n " + index_file + " -sel "+select_2+" -o "+output_file+" \n"  # shows numbers of H-bonds over time

    bash.execute(rdf_command)
    return output

def sasa(structure_file, ref_file, index_file, select_1, output, gmx="gmx", residual=False):
    """
    Args:
        structure_file:
        ref_file:
        index_file:
        select_1:
        output:
        gmx:
        residual:
    """

    if(residual):
        sasa_res = "-or "+str(residual)+"\n"
    else:
        sasa_res =""

    if (not os.path.isdir(output)):
        output_file =output
        output_dir = os.path.dirname(output)
    else:
        output_file = output+"/rdf_"+select_1
        output_dir = output

    selection = "echo \""+select_1+"\" "
    sasa_command = "cd  " + output_dir + " \n" + selection + " | " +gmx+" sasa -f "+structure_file +" -s "+ref_file+" -n "+index_file+" -o "+ output_file +" "+sasa_res+"\n"

    bash.execute(sasa_command)
    return output_file
