

from pharein.uniform_model import UniformModel
from pharein.diagnostics import FluidDiagnostics, ElectromagDiagnostics
from pharein.simulation import Simulation

from pharein.globals import objects

import os



def prepare_job():
    """
    prepare a simulation for a run:
        - creates the simulation directory [simulation.path] if it does not exist yet
        - writes an INI file in the directory [simulation.path]
        - for each diagnostic registered in the simulation:
            - creates the diagnostic directory 'diag['path']' under [simulation.path]/
              if it does not exist yet
    """

    if "Simulation" not in objects:
        raise RuntimeError("No Simulation created")


    sim = objects["Simulation"]

    print("preparing job...")
    if not os.path.exists(sim.path):
        print("mkdir "+sim.path)
        os.makedirs(sim.path)
    print("writing ini file in " + sim.path + os.path.sep)
    sim.write_ini_file()

    for diag in sim.diagnostics:
        path = diag.path
        full_path = os.path.join(sim.path, path)
        if not os.path.exists(full_path):
            print("mkdir " + full_path)
            os.makedirs(full_path)

