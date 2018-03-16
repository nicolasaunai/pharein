#!/usr/bin/env python

import pharein as ph


ph.Simulation(
    time_step_nbr=1000,
    boundary_types="periodic",
    cells=80,
    dl=0.1,
    final_time=1.,
    path='test5'
    #, refinement = {"level":[0,1], "extent_ratio":[0.4, 0.6], "refinement_iterations":[0, 3]}
)


ph.UniformModel(proton1={}, proton2={"density":2, "vbulk":(1.,0.,0.)})


ph.FluidDiagnostics(
    name="FluidDiagnostics1",
    diag_type="rho_s",
    write_every=10,
    compute_every=5,
    start_iteration=0,
    last_iteration=990,
    species_name="proton1"
)


ph.FluidDiagnostics(
    name="FluidDiagnostics2",
    diag_type="flux_s",
    write_every=10,
    compute_every=5,
    start_iteration=0,
    last_iteration=990,
    species_name="proton1"
)


ph.ElectromagDiagnostics(
    name="ElectromagDiagnostics1",
    diag_type="E",
    write_every=10,
    compute_every=5,
    start_teration=0,
    last_iteration=990
)


ph.ElectromagDiagnostics(
    name="ElectromagDiagnostics2",
    diag_type="B",
    write_every=10,
    compute_every=5,
    start_teration=0,
    last_iteration=990
)


ph.ParticleDiagnostics(
        name = "ParticleDiagnostics1",
        compute_every=10,
        write_every=10,
        start_iteration=0,
        last_iteration=90,
        diag_type="space_box",
        extent=(2., 4.),
        species_name="proton1"
        )


ph.prepare_job()

#phare.MiniPHARE().run(simu)


