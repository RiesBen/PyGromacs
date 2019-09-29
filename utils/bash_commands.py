"""
.. module:            generalFunctions
Description:
   This lib contains bash-wrappers for general use in a workflow, with nice error messages and other features.
   Additionally there are some convenience functions for error managment and dependencies check

Author: Benjamin Schroeder
"""

import io
import os
import sys
import glob
import typing as t
import subprocess as sub

#################################
#   General functions:
## Error Managment
def increment_error_level(err_prefix, old_err):
    """
    Args:
        err_prefix: add new information in which section something crashed
        old_err: pass on old err, message.
    """
    old_messages = ""
    for x in old_err.args:
        old_messages += "\t" + str(x)
    new_ex = Exception(str(err_prefix) + "\t See line: " + str(sys.exc_info()[-1].tb_lineno) + " \n" + old_messages)
    return new_ex

def print_error(err: Exception):
    """
    Args:
        err (Exception):
    """
    text= "ERROR: \t ".join(err.args)
    print(text,  file=sys.stderr)

## dependency checks.
def check_dependencies(depend:t.Union[t.Dict[any, str], t.List[str]], verbose:bool=True)->int:
    """
    .. function: check_deppendencies

    :
        checks a list of dependencies if each path is present or not. throws an
        IOError, if an Path is not existing.

    :rtype:int

    Args:
        depend (t.Union[t.Dict[any, str]):
        verbose (bool):

    Returns:
        returns int
    """
    found_error = False
    missing = []
    if(verbose and type(depend) is list):
        print("\n\n==================\n\tCHECK dependencies\n")
        print("\n".join(list(map(lambda s: "Check "+str(s), depend))))
    elif(verbose and type(depend) is dict):
        print("\nCHECK dependencies")
        print("\n".join(list(map(lambda s: "Check " + str(s), [depend[x] for x in depend]))))

    if (type(depend) is dict):
        for x in depend:
            if verbose: print(x)
            if ( depend[x] is str and not os.path.exists(depend[x])):
                found_error = True
                missing.append(x)

    elif (type(depend) is list):
        for x in depend:
            if verbose: print(x)
            if (not os.path.exists(x)):
                found_error = True
                missing.append(x)

    if found_error:
        print("\n==================\nAUTSCH\n==================\n")
        missing_str="\n\t".join(missing)
        raise IOError("COULD NOT FIND all DEPENDENCY!\n\t Could not find path to: \n\t" + str(missing_str),"\n\n")
    elif verbose:
        print("All dependencies are correct!", "\n\n")
    return 0

#################################
#   bash wrapper:
def extract_tar(in_path:str, out_path:str)->(str, io.FileIO):
    """
    Args:
        in_path (str):
        out_path (str):
    """
    option="-"
    option+="xf"
    command = "cd " + os.path.dirname(out_path) + " && tar "+option+" " + out_path + " "+str(in_path)+" && rm -r "+in_path
    ret = execute(command)


def  compress_tar(in_path:str, out_path:str=None, gnuzip:bool=False)->(str, io.FileIO):

    """
    Args:
        in_path (str):
        out_path (str):
        gnuzip (bool):
    """
    option= ""
    if (out_path==None):
        out_path=in_path
    if(not out_path.endswith(".tar.gz") and gnuzip):
        out_path+=".tar.gz"
    elif(not out_path.endswith(".tar")):
        out_path+=".tar"

    if(gnuzip):
        option+="z"

    print(in_path)
    print(out_path)
    option+="cf"
    command = "cd " + os.path.dirname(out_path) + " && tar "+option+" " + out_path + " "+str(in_path)+" && rm -r "+in_path
    ret = execute(command)
    return out_path, ret

def concat(outFile:str, infiles:list, verbose:bool=False):
    """
    Args:
        outFile (str):
        infiles (list):
    """
    command = "cat "+ " ".join(infiles)+ " > "+outFile+" \n"
    if verbose: print("CONCAT: "+command)
    if(os.path.exists(os.path.dirname(outFile))):
        if(os.system(command)):
            raise Exception("could not concate Files:\n " + str(" ".join(infiles))+"\n to\n"+outFile)
    else:
        raise IOError("could not find folder for:\n "+outFile)
    return outFile

def make_folder(path, option="", verbose=False):
    """
    Args:
        path:
        option:
        verbose:
    """
    mk_folder = "mkdir " + option + " " + str(path) + "\n"
    if (not os.path.isdir(path)):
        if (os.system(mk_folder)):
            raise Exception("could not make folder:\n " + str(path))

    elif(verbose):
        print("Warning! Did not build already existing folder: "+path)
    return path

def remove_folder(path:str, option:str="", verbose:bool=False):
    """
    Args:
        path:
        option:
        verbose:
    """
    mk_folder = "rmdir " + option + " " + str(path) + "\n"
    if (os.path.isdir(path)):
        if (os.system(mk_folder)):
            raise Exception("could not remove folder:\n " + str(path))

    elif(verbose):
        print("Warning! Did not remove non existing folder: "+path)
    return path


def save_make_folder(path, option=""):
    """
    Args:
        path:
        option:
    """

    offset = 1
    dir_versions = list(filter(lambda x: not ".tar" in x or not ".gz" in x,  sorted(glob.glob(path+"*"))))
    print("dirversions:", dir_versions)
    if(len(dir_versions)>0):
        last_dir = dir_versions[len(dir_versions)-1]
        print("last:", last_dir)

        suffix = str(last_dir.replace(path+"_", ""))
        print(suffix)
        if(suffix != "" and suffix.isalnum()):
            offset = int(suffix)+1
        path += "_" + str(offset)

    mk_folder = "mkdir " + option + " " + str(path) + "\n"
    if (os.system(mk_folder)):
        raise Exception("could not make folder:\n " + str(path))

    return path

def copy_file(origin, target, option=""):
    """
    Args:
        origin:
        target:
        option:
    """
    copy_files = "cp " + str(option) + " " + str(origin) + " " + str(target) + "\n"
    if (os.system(copy_files)):
        raise Exception("could not copy:\n " + str(origin) + "\n \t to \n" + str(target) + "\n \t options: " + str(option))
    return target

def link_folder(origin, target, option=""):
    """
    Args:
        origin:
        target:
        option:
    """
    link_folders = "ln -s " + option + " " + origin + " " + target + "\n"
    if (not os.path.exists(target)):
        if (os.system(link_folders)):
            raise Exception("could not link:\n " + str(origin) + "\n \t to \n" + str(target) + "\n \t options: " + str(option))


def mv_file(origin, target, option="")->str:
    """
    Args:
        origin:
        target:
        option:
    """
    copy_files = "mv " + option + " " + origin + " " + target + "\n"
    if (os.system(copy_files)):
        raise Exception("could not copy:\n " + str(origin) + "\n \t to \n" + str(target) + "\n \t options: " + str(option))
    return target


def execute(command: (list or str))-> io.FileIO:
    """
    Args:
        command:
    """
    if(type(command) == list):
        command = " ".join(command)
    ret = -1
    try:
        ret = os.popen(command)
        os.wait()
    except:
       raise Exception("could not execute:\n " + str(command)+ "\n\tCommand returned: \t"+str(ret.read()))

    return ret


def execute_sub(command:(str or t.List[str]), verbose:bool=False)-> io.FileIO:
    """
        not very robust, crashes with gromos!

    Args:
        command:
        verbose (bool):
    """
    if(type(command)==str):
        command = command.split()
        pass

    if(verbose): print("EXECUTING: \t "+"".join(str(command)))
    try:
        proc = sub.Popen(command, stdout=sub.PIPE, stderr=sub.PIPE)
    except OSError as err:
        raise Exception("\n   bash.execute:Process opening Error!"+str(err.args)+"\n Output: ")
    proc.wait()

    if(proc.returncode != 0):
        raise Exception("\n   bash.execute:Could not execute(returnCode:"+str(proc.returncode)+"):\n \t" + str(" ".join(command)) + "\n\n")

    if(verbose): print("\nSTDOUT\n=======",proc.stdout.read().decode("utf-8"), "\nSTDERR\n=======","\n",proc.stderr.read().decode("utf-8") )
    return proc.returncode


def rm_file(target, options=""):
    """
    Args:
        target:
        options:
    """
    rm_command = "rm " + str(target) + " " + str(options)
    if (os.path.exists(target)):
        if (os.system(rm_command)):
            raise Exception("could not delete file/folder: " + str(target) + "\n options: " + str(options))

def sed(in_file, find_pattern, replace_pattern, out_file=False):
    """
    Args:
        in_file:
        find_pattern:
        replace_pattern:
        out_file:
    """
    command =""
    if out_file:
       command = "sed s/"+find_pattern+"/"+replace_pattern+"/g "+in_file+" > "+out_file + " \n"
    else:
        command = "sed s/" + find_pattern + "/"+replace_pattern + "/g " + in_file + " \n"

    os.system(command)

