import os
import importlib
import unittest


if __name__ == "__main__":
    root_dir = os.path.dirname(__file__)
    test_files = []
    print(root_dir)
    #gather all test_files
    for dir in os.walk(root_dir):
        test_files.extend([dir[0]+"/"+path for path in dir[2] if( path.startswith("test") and path.endswith(".py"))])

    print(test_files)

    #do tests:
    suite = unittest.TestSuite()
    modules = []
    for test in test_files:
        module_name = "generalutilities."+test.replace(os.path.dirname(root_dir), "").replace("/", ".").replace(".py", "")
        modules.append(module_name)

    print(modules)

    for test in modules:
        try:
            # If the module defines a suite() function, call it to get the suite.
            mod = __import__(test, globals(), locals(), ['suite'])
            suitefn = getattr(mod, 'suite')
            suite.addTest(suitefn())
        except (ImportError, AttributeError):
            # else, just load all the test cases from the module.
            suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test))

    try:
        unittest.TextTestRunner().run(suite)
        exit(0)
    except:
        exit(1)
