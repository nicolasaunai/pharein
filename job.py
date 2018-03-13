#!/usr/bin/env python

import sys
import miniphare.inputs.uniform_model as uniform_model
import miniphare.inputs.simulation as simulation
import miniphare.inputs.diagnostics as diagnostics
import miniphare.inputs.phare as phare

def main():
    simu = simulation.Simulation(time_step_nbr=1000, boundary_types="periodic", cells=80, dl=0.1, final_time=1.,
                                 path='simu1')

    fd1 = diagnostics.FluidDiagnostics(name="FluidDiagnostics1", diag_type="rho_s",
                                       write_every=10, compute_every=5, start_iteration=0,
                                       last_iteration=990, species_name="proton1")

    fd2 = diagnostics.FluidDiagnostics(name="FluidDiagnostics2",
                                       diag_type="flux_s", write_every=10, compute_every=5,
                                       start_iteration=0, last_iteration=990, species_name="proton1")

    ed1 = diagnostics.ElectromagDiagnostics(name="ElectromagDiagnostics1", diag_type="E",
                                            write_every=10, compute_every=5, start_teration=0, last_iteration=990)

    ed2 = diagnostics.ElectromagDiagnostics(name="ElectromagDiagnostics2", diag_type="B",
                                            write_every=10, compute_every=5, start_teration=0, last_iteration=990)

    simu.add_diagnostics(fd1)
    simu.add_diagnostics(fd2)
    simu.add_diagnostics(ed1)
    simu.add_diagnostics(ed2)

    model = uniform_model.UniformModel()
    model.add_fields()
    model.add_species("proton1")
    model.add_species("proton2", density=2.)

    simu.set_model(model)

    simu.refine([0, 1], (0.4, 0.6), (0, 3))

    simulation.prepare_job(simu)

    #phare.MiniPHARE().run(simu)


if __name__ == "__main__":

    main()

