import subprocess
import os
import sys
import pkg_resources


def Check_Environments():
    print("Checking Python environments...")
    print("Python path: {}".format(sys.executable))
    print("Python version: {}".format(sys.version))
    print("Python venv: {}".format(sys.prefix))
    print("Working Dictionary: {}".format(os.path.dirname(os.path.abspath(__file__))))
    # check R environment
    try:
        print("Checking R environments...")
        os.system("R --version")
    except OSError:
        print("R ecosystem remains unsolved.")
        sys.exit(-1)


def Check_Requirements(required_packages):
    print("Checking basical requirements...")
    required = required_packages
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed

    if missing:
        print("Missing requirements:{}".format(missing))
        python = sys.executable
        print("Installing missing requirements...")
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        print("Missing requirements installed.")
        return 0
    else:
        print("All requirements satisfied.")
        return 0
