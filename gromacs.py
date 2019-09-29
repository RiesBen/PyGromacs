"""
FUNCTIONLIB:            gromacs Object
Description:
   This lib contains gromacs objects for diff. gromacs versions!

Author: Benjamin Schroeder
"""

#imports
import os
import generalutilities.function_libs.utils.bash_commands as bash


def get_gmx(version):
    """
    Args:
        version:
    """
    if type(version) == str:
        if (version.startswith("4")):
            return _Gromacs_4(version)
        if (version.startswith("5")):
            return _Gromacs_5(version)
        if (version.startswith("2018")):
            return _Gromacs_2018(version)
        else:
            raise ValueError("Unknown Gromacs version! ")
    else:
        raise IOError("Please give verstion as a string! (e.g. \"4.6.7\")")


class _Gromacs:
    def tell_version(self):
        print("My Version is: "+self.version)

    #################################
    #   gromacs functions:
    # ---------------General
    def pdb2gmx(self, pdb_file, force_field, watermodel, output_path=False, gromacs="gmx", additional_options=" "):
        """
        Args:
            pdb_file:
            force_field:
            watermodel:
            output_path:
            gromacs: path to gromacs binary (default: gmx)
            additional_options:
        """
        print("converting pdb to gmx")
        output_path = output_path if (output_path) else os.path.splitext(pdb_file)[0]

        if(str(force_field).isdigit()):
            pdb2gmx_command = "cd " + os.path.dirname(output_path) + "\n  echo \"" + str(force_field) + "\" | " + gromacs + " pdb2gmx -f " + pdb_file + " -water " + watermodel + " -o " + output_path + ".gro -p " + output_path + " -n " + output_path + " " + additional_options
        else:        
            pdb2gmx_command = "cd " + os.path.dirname(output_path) + "\n " + gromacs + " pdb2gmx -f " + pdb_file + " -ff " + force_field + " -water " + watermodel + " -o " + output_path + ".gro -p " + output_path + " -n " + output_path + " " + additional_options
        bash.execute(pdb2gmx_command)
        return output_path

    def editconf(self, input_file, output_name=False, gromacs="gmx", additional_options=" "):
        """
        Args:
            input_file:
            output_name:
            gromacs: path to gromacs binary (default: gmx)
            additional_options:
        """

        output_name = output_name if (output_name) else os.path.splitext(input_file)[0]
        editconf_command = "cd " + os.path.dirname(
            output_name) + "\n " + gromacs + " editconf -f " + input_file + " -o " + output_name + " " + additional_options
        bash.execute(editconf_command)
        return output_name

    def set_box(self, gro_file, border_distance, box_form="dodecahedron", output_name=False, gromacs="gmx", additional_options=""):
        """
        Args:
            gro_file:
            border_distance:
            box_form:
            output_name:
            gromacs: path to gromacs binary (default: gmx)
            additional_options:
        """
        output_name = str(output_name) if (output_name) else gro_file
        # box command for gmx editconf
        editconf_command = gromacs + " editconf -f " + gro_file + " -o " + output_name + ".gro -bt " + box_form + " -d " + str(
            border_distance) + " " + additional_options + " &> " + output_name + ".log"
        bash.execute(editconf_command)
        return output_name

    def solvate_box(self, gro_file, top_file, output_name=False, watersNum="none", gromacs="gmx"):
        """
        Args:
            gro_file:
            top_file:
            output_name:
            watersNum:
            gromacs: path to gromacs binary (default: gmx)
        """
        output_name = str(output_name) if (output_name) else os.path.splitext(gro_file)[0]
        if (watersNum != "none"):
            solvate_command = gromacs + " solvate -cp " + gro_file + " -cs spc216.gro -o " + output_name + ".gro -maxsol " + str(
                watersNum) + " -p " + top_file + " &> " + output_name + ".log \n"
        else:
            solvate_command = gromacs + " solvate -cp " + gro_file + " -cs spc216.gro -o " + output_name + ".gro -p " + top_file + ".top &> " + output_name + ".log \n"

        bash.execute(solvate_command)
        return output_name

    def add_ions(self, gro_file, top_file, group_sel, output_name, conc="0.15", np=False, nn=False, additional_options="",
                 gromacs="gmx"):
        """
        Args:
            gro_file:
            top_file:
            group_sel:
            output_name:
            conc:
            np:
            nn:
            additional_options:
            gromacs: path to gromacs binary (default: gmx)
        """
        folder = os.path.dirname(output_name)
        ion_mdp_path = folder + "/ion"

        # generating an ion_protocol
        print("writing dummy protocol to: " + ion_mdp_path)
        ion_mdp = open(ion_mdp_path + ".mdp", "w")
        ion_mdp.write(";this is a generic em-protocol to execute grompp\n"
                      "integrator      = steep\n"
                      "emtol           = 250.0\n"
                      "nsteps          = 50000\n"
                      "nstenergy       = 1\n"
                      "energygrps      = System\n"
                      "nstlist         = 1\n"
                      "ns_type         = grid\n"
                      "coulombtype     = PME\n"
                      "rlist           = 1.0\n"
                      "rcoulomb        = 1.0\n"
                      "rvdw            = 1.0\n"
                      "constraints     = none\n"
                      "pbc             = xyz")
        ion_mdp.close()

        print("grompp")
        # salting command for gmx genion
        grompp = ("cd " + folder + "\n" +
                  gromacs + " grompp -f " + ion_mdp_path + " -c " + gro_file + " -p " + top_file + " -o " + output_name + ".tpr 1> " + output_name + "_grompp.log 2> " + output_name + "_grompp.err\n")
        bash.execute(grompp)

        print("salting")
        if (np and nn):
            add_salt = ("cd " + folder + "\necho -e " + str(
                group_sel) + " | " + gromacs + " genion -s " + output_name + ".tpr -o " + output_name + ".gro -nn " + str(
                nn) + " -np " + str(
                np) + " -p " + output_name + ".top " + additional_options + " 1> " + output_name + "_mdrun.log 2> " + output_name + "_mdrun.err\n")
        else:
            add_salt = ("cd " + folder + "\necho -e " + str(
                group_sel) + " | " + gromacs + " genion -s " + output_name + ".tpr -o " + output_name + ".gro -neutral -conc " + str(
                conc) + " -p " + output_name + ".top " + additional_options + " 1> " + output_name + ".log 2> " + output_name + "_mdrun.err\n")
        bash.execute(add_salt)
        return output_name
    def grompp(self, workdir, gro_file, top_file, ndx_file, in_protocol, gromacs="gmx", out_prefix=False):
        """This starts grompp to build a tpr file. :param workdir: absolute path
        to folder with step prefix as last part ("absolutePath/prefix") :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder) :return: stepcounter (updates the gro_file for next
        step)

        Args:
            workdir:
            gro_file:
            top_file:
            ndx_file:
            in_protocol:
            gromacs:
            out_prefix:
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(in_protocol))[0]  # prefix of the step
        output = workdir + "/" + step_name

        # Formulate commands
        grompp_command = gromacs + " grompp -f " + in_protocol + " -c " + gro_file + " -p " + top_file + " -n " + ndx_file + " -o " + output + ".tpr 1> " + output + "_grompp.log 2>" + output + "_grompp.err\n"
        try:
            # execute commands:
            #print("\tCopy_protocol: " + step_name)
            #bash.copy_file(in_protocol, output + ".mdp")

            print("\tgrompp command: " + step_name)
            bash.execute("cd " + workdir + "\n" + grompp_command)

        except Exception as err:
            raise bash.increment_error_level("grompp failed!\n", err)

        return output

    def mdrun(self, workdir, tpr_file, out_prefix=False, gromacs="gmx ", additional_options:str=" "):
        """This starts an initial run for a simulation "step". very crude and
        easy ;) :param workdir: absolute path to folder with step prefix as last
        part ("absolutePath/prefix") :param gromacs_run: contains the string,
        executing the gromacs command: "gmx mdrun" "gmx_mpi mdrun" ... :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder) :return: stepcounter (updates the gro_file for next
        step)

        Args:
            workdir:
            tpr_file:
            out_prefix:
            gromacs:
            additional_options (str):
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(tpr_file))[0]  # prefix of the step
        output = workdir + "/" + step_name

        # Formulate commands
        mdrun_command = gromacs + " mdrun -s " + tpr_file + " -deffnm " + output + " " +additional_options+" 1> " + output + "_mdrun.log 2>" + output + "_mdrun.err"
        try:
            # execute commands:
            print("\tCopy_protocol: " + step_name)
            if(not os.path.exists(output+".tpr")):
                bash.copy_file(tpr_file, output + ".tpr")

            print("\trun_sim: " + step_name)
            bash.execute("cd " + workdir + " && " + mdrun_command)
        except Exception as err:
            raise bash.increment_error_level("mdrun failed!\n", err)
        return output

    def md(self, workdir, gro_file, top_file, ndx_file, in_protocol, gromacs_run, gromacs_grompp, out_prefix=False):
        """
            IMPROVE by plugging two subcommands together: grompp and mdrun

        This starts an initial run for a simulation "step". very crude and
        easy ;) :param workdir: absolute path to folder with step prefix as last
        part ("absolutePath/prefix") :param gromacs_run: contains the string,
        executing the gromacs command: "gmx mdrun" "gmx_mpi mdrun" ... :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder) :return: stepcounter (updates the gro_file for next
        step)

        Args:
            workdir:
            gro_file:
            top_file:
            ndx_file:
            in_protocol:
            gromacs_run:
            gromacs_grompp:
            out_prefix:
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(in_protocol))[0]  # prefix of the step
        output = workdir + "/" + step_name

        # Formulate commands
        grompp_command = gromacs_grompp + " -f " + output + ".mdp -c " + gro_file + " -p " + top_file + " -n " + ndx_file + " -o " + output + ".tpr 1> " + output + "_grompp.log 2>" + output + "_grompp.err\n"
        mdrun_command = gromacs_run + " -s " + output + ".tpr -deffnm " + output + " 1> " + output + "_mdrun.log 2>" + output + "_mdrun.err\n"
        try:
            # execute commands:
            print("\tCopy_protocol: " + step_name)
            bash.copy_file(in_protocol, output + ".mdp")

            print("\tgrompp command: " + step_name)
            bash.execute("cd " + workdir + "\n" + grompp_command)

            print("\trun_sim: " + step_name)
            bash.execute("cd " + workdir + "\n" + mdrun_command)
        except Exception as err:
            raise bash.increment_error_level("simple sim failed!\n", err)
        return output

    def make_ndx(self, gro_file, gromacs="gmx", index_file=False, output_name=False, selection=False):
        """
        Args:
            gro_file:
            gromacs:
            index_file:
            output_name:
            selection:
        """
        if (not output_name):
            output_name = gro_file

        if (selection):
            selection_string = "\"" + str(selection) + " \n q\""
            print("Got selection " + selection_string)
        else:
            selection_string = "\"q\""
        if (index_file):
            gro_file += " -n " + index_file + " "
        print("making index file:\t" + output_name)
        make_ndx_command = gromacs + " make_ndx -f " + gro_file + " -o " + output_name + " <<< $(echo -e \""+selection_string+"\")"
        bash.execute(make_ndx_command)

        return output_name

    def make_posres(self, gro_file, output_name, selection="1", fc="1000", gromacs="gmx"):
        """
        Args:
            gro_file:
            output_name:
            selection:
            fc:
            gromacs:
        """

        print("making posrestriction file:\t" + output_name)
        make_posres_command = "echo \"" + selection + " 0\" | " + gromacs + " genrestr -f " + gro_file + " -o " + output_name + ".itp -fc " + fc + " \n"
        bash.execute(make_posres_command)
        return output_name

    def trajectory_convert(self, struct_file, s_file, gromacs="gmx", output_name=False, format=False, index_file=False,
                           frames_file=False, begin_frame=False, end_frame=False, time_unit=False, timestep=False,
                           split=False, pbc=False, unit_cell_rep=False, center=False, dt=False, fit=False,
                           additional_options=False, group_sel="0"):
        """
        Args:
            struct_file:
            s_file:
            gromacs: path to gromacs binary (default: gmx)
            output_name:
            format:
            index_file:
            frames_file:
            begin_frame:
            end_frame:
            time_unit:
            timestep:
            split:
            pbc:
            unit_cell_rep:
            center:
            dt:
            fit:
            additional_options:
            group_sel:
        """

        if (not output_name):
            output_name = struct_file
        if (format):
            output_name = os.path.splitext(output_name)[0] + "." + str(format)
        if (not group_sel):
            group_sel = 0

        options = " "
        options += " -n " + str(index_file) + " " if (index_file) else ""
        options += " -fr " + str(frames_file) + " " if (frames_file) else ""
        options += " -b " + str(begin_frame) + " " if (begin_frame) else ""
        options += " -e " + str(end_frame) + " " if (end_frame) else ""
        options += " -tu " + str(time_unit) + " " if (time_unit) else ""
        options += " -timestep " + str(timestep) + " " if (timestep) else ""
        options += " -dt " + str(dt) + " " if (dt) else ""
        options += " -split " + str(split) + " " if (split) else ""
        options += " -pbc " + str(pbc) + " " if (pbc) else ""
        options += " -ur " + str(unit_cell_rep) + " " if (unit_cell_rep) else ""
        options += " -center " if (center) else ""
        options += " -fit " + str(fit) + " " if (fit) else ""
        options += additional_options if (additional_options) else ""

        trjconv_command = "echo -e \"" + group_sel + "\" | " + gromacs + " trjconv -f " + struct_file + " -s " + s_file + " -o " + output_name + " " + options + "\n"
        bash.execute(trjconv_command)
        return output_name

class _Gromacs_4(_Gromacs):
    def __init__(self, version):
        """
        Args:
            version:
        """
        self.version = version
        print("Careful with gmx4! It has several independent gmx bins and is not ordered with the gmx prefix!")

    def md(self, workdir, gro_file, top_file, ndx_file, in_protocol, gromacs_run="mdrun_mpi", gromacs_grompp="grompp_mpi", out_prefix=False):
        """This starts an initial run for a simulation "step". very crude and
        easy ;) :param workdir: absolute path to folder with step prefix as last
        part ("absolutePath/prefix") :param gromacs_run: contains the string,
        executing the gromacs command: "gmx mdrun" "gmx_mpi mdrun" ... :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder)

        :return output_path

        Args:
            workdir:
            gro_file:
            top_file:
            ndx_file:
            in_protocol:
            gromacs_run:
            gromacs_grompp:
            out_prefix:
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(in_protocol))[0]  # prefix of the step
        output = workdir + "/" + step_name

        # Formulate commands
        grompp_command = gromacs_grompp + " -f " + output + ".mdp -c " + gro_file + " -p " + top_file + " -n " + ndx_file + " -o " + output + ".tpr 1> " + output + "_grompp.log 2>" + output + "_grompp.err\n"
        mdrun_command = gromacs_run + " -s " + output + ".tpr -o "+out_prefix+".trr -c "+out_prefix+".gro -e "+out_prefix+".edr 1> " + output + "_mdrun.log 2>" + output + "_mdrun.err\n"
        try:
            # execute commands:
            print("\tCopy_protocol: " + step_name)
            bash.copy_file(in_protocol, output + ".mdp")

            print("\tgrompp command: " + step_name)
            bash.execute("cd " + workdir + "\n" + grompp_command)

            print("\trun_sim: " + step_name)
            bash.execute("cd " + workdir + "\n" + mdrun_command)
            return output
        except Exception as err:
            raise bash.increment_error_level("simple sim failed!\n", err)

    def grompp(self, workdir, gro_file, top_file, ndx_file, in_protocol, gromacs_grompp="grompp_mpi", out_prefix=False):
        """This starts grompp to build a tpr file. :param workdir: absolute path
        to folder with step prefix as last part ("absolutePath/prefix") :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder) :return: stepcounter (updates the gro_file for next
        step)

        Args:
            workdir:
            gro_file:
            top_file:
            ndx_file:
            in_protocol:
            gromacs_grompp:
            out_prefix:
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(in_protocol))[0]  # prefix of the step
        output = workdir + "/" + step_name

        # Formulate commands
        grompp_command = gromacs_grompp + " -f " + output + ".mdp -c " + gro_file + " -p " + top_file + " -n " + ndx_file + " -o " + output + ".tpr 1> " + output + "_grompp.log 2>" + output + "_grompp.err\n"
        try:
            # execute commands:
            print("\tCopy_protocol: " + step_name)
            bash.copy_file(in_protocol, output + ".mdp")

            print("\tgrompp command: " + step_name)
            bash.execute("cd " + workdir + "\n" + grompp_command)

        except Exception as err:
            raise bash.increment_error_level("grompp failed!\n", err)

        return output

    def make_ndx(self, gro_file, gromacs="make_ndx_mpi", index_file=False, output_name=False, selection=False):
        """
        Args:
            gro_file:
            gromacs:
            index_file:
            output_name:
            selection:
        """

        if (not output_name):
            output_name = gro_file

        if (selection):
            selection_string = "\"" + str(selection) + " \n q\""
            print("Got selection " + selection_string)
        else:
            selection_string = "\"q\""
        if (index_file):
            gro_file += " -n " + index_file + " "

        print("making index file:\t" + output_name)
        make_ndx_command = gromacs + " -f " + gro_file + " -o " + output_name + " <<< $(echo -e \""+selection_string+"\")"
        bash.execute(make_ndx_command)

        return output_name

class _Gromacs_5(_Gromacs):
    def __init__(self, version):
        """
        Args:
            version:
        """
        self.version = version

class _Gromacs_2018(_Gromacs):
    def __init__(self, version):
        """
        Args:
            version:
        """
        self.version = version

    def md(self, workdir,  gro_file, top_file, ndx_file, in_protocol, out_prefix=False, restriction_gro=False, gromacs_run="gmx mdrun", gromacs_grompp="gmx grompp"):
        """This starts an initial run for a simulation "step". very crude and
        easy ;) :param workdir: absolute path to folder with step prefix as last
        part ("absolutePath/prefix") :param gromacs_run: contains the string,
        executing the gromacs command: "gmx mdrun" "gmx_mpi mdrun" ... :param
        gro_file:absolute path to .gro file input (build input or previous
        simulation step) :param top_file: absolute path to .top file input
        (build input folder) :param ndx_file: absolute path to .ndx file input
        (build input folder) :return: stepcounter (updates the gro_file for next
        step)

        Args:
            workdir:
            gro_file:
            top_file:
            ndx_file:
            in_protocol:
            out_prefix:
            restriction_gro:
            gromacs_run:
            gromacs_grompp:
        """
        if (out_prefix):
            step_name = out_prefix
        else:
            step_name = os.path.splitext(os.path.basename(in_protocol))[0]  # prefix of the step
        if(not restriction_gro):
            restriction_gro=gro_file
        output = workdir + "/" + step_name

        # Formulate commands
        grompp_command = gromacs_grompp + " -f " + output + ".mdp -c " + gro_file + " -p " + top_file + " -n " + ndx_file + " -r "+restriction_gro+" -o " + output + ".tpr 1> " + output + "_grompp.log 2>" + output + "_grompp.err\n"
        mdrun_command = gromacs_run + " -s " + output + ".tpr -deffnm " + output + " 1> " + output + "_mdrun.log 2>" + output + "_mdrun.err\n"
        try:
            # execute commands:
            print("\tCopy_protocol: " + step_name)
            bash.copy_file(in_protocol, output + ".mdp")

            print("\tgrompp command: " + step_name)
            bash.execute("cd " + workdir + "\n" + grompp_command)

            print("\trun_sim: " + step_name)
            bash.execute("cd " + workdir + "\n" + mdrun_command)
        except Exception as err:
            raise bash.increment_error_level("simple sim failed!\n", err)
        return output





