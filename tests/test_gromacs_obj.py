from PyGromacs import gromacs as gmx

try:
    gromacs = gmx.get_gmx(4)
except Exception as err:
    if(type(err)== IOError):
        print("input str check!")
    else:
        print("got wrong input!")

try:
    gromacs = gmx.get_gmx("4")
    print(gromacs.version)
except Exception as err:
    print("something should not have went wrong!"+str(err.args))
