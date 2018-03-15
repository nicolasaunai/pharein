

import numpy as np
import configparser
import os

import pharein.phare_utilities as phare_utilities


def compute_dimension(cells):
    return len(cells)


# ------------------------------------------------------------------------------


def check_domain(**kwargs):
    """
    check parameters relative to the domain and raise exceptions if they are not valid
    return 'dl' and 'cells'
    """
    dom_and_dl = ('domain_size' in kwargs and 'dl' in kwargs and 'cells' not in kwargs)
    cells_and_dl = ('cells' in kwargs and 'dl' in kwargs and 'domain_size' not in kwargs)
    dom_and_cells = ('cells' in kwargs and 'domain_size' in kwargs and 'dl' not in kwargs)

    if not (dom_and_dl or dom_and_cells or cells_and_dl):
        raise ValueError("Error: Specify either 'domain_size' and 'dl' or 'cells' and 'dl' or 'domain_size' and 'dl'")

    if dom_and_dl:
        try:
            phare_utilities.check_iterables(kwargs['domain_size'], kwargs['dl'])
            phare_utilities.check_equal_size(kwargs['domain_size'], kwargs['dl'])
        except ValueError:
            raise ValueError("Error: 'domain_size' and 'dl' must have the same dimension")

        domain_size = phare_utilities.listify(kwargs['domain_size'])
        dl = phare_utilities.listify(kwargs['dl'])

        cells = [int(dom_size / float(mesh_size)) for dom_size, mesh_size in zip(domain_size, dl)]

        for s, d in zip(domain_size, dl):
            if s < d:
                raise ValueError("Error: mesh size should be smaller than domain_size")

        dl = [dom_size / n for dom_size, n in zip(domain_size, cells)]

    elif cells_and_dl:
        try:
            phare_utilities.check_iterables(kwargs['cells'], kwargs['dl'])
            phare_utilities.check_equal_size(kwargs['cells'], kwargs['dl'])
        except ValueError:
            raise ValueError("Error: 'cells' and 'dl' must have the same dimension")

        dl = phare_utilities.listify(kwargs['dl'])
        cells = phare_utilities.listify(kwargs['cells'])

    elif dom_and_cells:
        try:
            phare_utilities.check_iterables(kwargs['cells'], kwargs['domain_size'])
            phare_utilities.check_equal_size(kwargs['cells'], kwargs['domain_size'])
        except ValueError:
            raise ValueError("Error: 'cells' and 'domain_size' must have the same dimension")

        cells = phare_utilities.listify(kwargs['cells'])
        domain_size = phare_utilities.listify(kwargs['domain_size'])
        dl = [dom_size / float(n) for dom_size, n in zip(domain_size, cells)]

    for n, d in zip(cells, dl):
        if n != 0 and n < 10:
            raise ValueError("Error: number of cells in non-invariant directions should be >= 10")
        if n == 0 and d != 0:
            raise ValueError("Error: dl should be 0 for invariant directions")
        if n < 0:
            raise ValueError("Error: number of cells must be >= 0")
        if d < 0:
            raise ValueError("Error: mesh size must be >= 0")

    return dl, cells


# ------------------------------------------------------------------------------


def check_time(**kwargs):
    """
    check parameters relative to simulation duration and time resolution and raise exception if invalids
     return time_step_nbr and time_step
    """
    final_and_dt = ('final_time' in kwargs and 'time_step' in kwargs and 'time_step_nbr' not in kwargs)
    nsteps_and_dt = ('time_step_nbr' in kwargs and 'time_step' in kwargs and 'final_time' not in kwargs)
    final_and_nsteps = ('final_time' in kwargs and 'time_step_nbr' in kwargs and 'time_step' not in kwargs)

    if final_and_dt:
        time_step_nbr = int(kwargs['final_time'] / kwargs['time_step'])
        time_step = kwargs['final_time'] / time_step_nbr

    elif final_and_nsteps:
        time_step = kwargs['final_time'] / kwargs['time_step_nbr']
        time_step_nbr = kwargs['time_step_nbr']

    elif nsteps_and_dt:
        time_step = kwargs['time_step']
        time_step_nbr = kwargs['time_step_nbr']

    else:
        raise ValueError("Error: Specify either 'final_time' and 'time_step' or 'time_step_nbr' and 'time_step'" + \
                         " or 'final_time' and 'time_step_nbr'")

    return time_step_nbr, time_step


# ------------------------------------------------------------------------------


def check_interp_order(**kwargs):
    interp_order = kwargs.get('interp_order', 1)

    if interp_order not in [1, 2, 3, 4]:
        raise ValueError("Error: invalid interpolation order. Should be in [1,2,3,4]")

    return interp_order

# ------------------------------------------------------------------------------


def check_pusher(**kwargs):
    pusher = kwargs.get('particle_pusher', 'modified_boris')
    if pusher not in ['modified_boris']:
        raise ValueError('Error: invalid pusher ({})'.format(pusher))
    return pusher


# ------------------------------------------------------------------------------


def check_layout(**kwargs):
    layout = kwargs.get('layout', 'yee')
    if layout not in ('yee'):
        raise ValueError('Error: invalid layout ({})'.format(layout))
    return layout


# ------------------------------------------------------------------------------


def check_path(**kwargs):
    path = kwargs.get('path', '.' + os.path.sep)
    return path

# ------------------------------------------------------------------------------


def check_boundaries(dims, **kwargs):
    valid_boundary_types = ("periodic",)
    boundary_types = kwargs.get('boundary_types', ['periodic'] * dims)
    phare_utilities.check_iterables(boundary_types)

    if phare_utilities.none_iterable(boundary_types):
        bc_length = 1
        if boundary_types not in valid_boundary_types:
            raise ValueError("Error: '{}' is not a valid boundary type".format(boundary_types))
        else:
            boundary_types = phare_utilities.listify(boundary_types)
    else:
        bc_length = len(boundary_types)
        for bc in boundary_types:
            if bc not in valid_boundary_types:
                raise ValueError("Error: '{}' is not a valid boundary type".format(bc))

    if bc_length != dims:
        raise ValueError("Error- boundary_types should have length {} and is of length {}".format(dims, bc_length))

    return boundary_types


# ------------------------------------------------------------------------------


def check_origin(dims, **kwargs):
    origin = kwargs.get("origin", [0] * dims)
    return origin

# ------------------------------------------------------------------------------


def add_dict_to_config(the_dict,conf, section_key="section"):
    """
    adds the dictionnary 'the_dict' to a config 'conf' under the section 'section_key'
    if the section does not exist, it is created
    """
    sections = conf.sections()
    this_section = the_dict[section_key]
    if this_section not in sections:
        conf.add_section(this_section)
    for key in the_dict.keys():
        if key != section_key:
            conf.set(the_dict[section_key], key, str(the_dict[key]))

# ------------------------------------------------------------------------------


def checker(func):
    def wrapper(simulation_object, **kwargs):
        accepted_keywords = ['domain_size', 'cells', 'dl', 'particle_pusher', 'final_time',
                             'time_step', 'time_step_nbr', 'layout', 'interp_order', 'origin',
                             'boundary_types', 'refined_particle_nbr', 'path',
                             'diag_export_format']

        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise ValueError("Error: invalid arguments - " + " ".join(wrong_kwds))
        try:
            dl, cells = check_domain(**kwargs)
            time_step_nbr, time_step = check_time(**kwargs)
            interp_order = check_interp_order(**kwargs)
            pusher = check_pusher(**kwargs)
            layout = check_layout(**kwargs)
            path = check_path(**kwargs)
            dims = compute_dimension(cells)
            boundary_types = check_boundaries(dims, **kwargs)
            origin = check_origin(dims, **kwargs)

            refined_particle_nbr = kwargs.get('refined_particle_nbr', 2)  # TODO change that default value
            diag_export_format = kwargs.get('diag_export_format', 'ascii') #TODO add checker with valid formats

            return func(simulation_object, cells=cells, dl=dl, interp_order=interp_order,
                        time_step=time_step, time_step_nbr=time_step_nbr,
                        particle_pusher=pusher, layout=layout, origin=origin,
                        boundary_types=boundary_types, path=path, refined_particle_nbr=refined_particle_nbr,
                        diag_export_format=diag_export_format)

        except ValueError as msg:
            print(msg)

    return wrapper


# ------------------------------------------------------------------------------


class Simulation(object):
    """represents a simulation input"""

    @checker
    def __init__(self, **kwargs):
        """
        1D run example: Simulation(time_step_nbr = 100, boundary_types="periodic", cells=80)
        2D run example: Simulation(time_step_nbr = 100, boundary_types=("periodic","periodic"), cells=(80,20))
        3D run example: Simulation(time_step_nbr = 100, boundary_types=("periodic","periodic","periodic"), cells=(80,20,40))

        optional parameters:
        -------------------

        dl                   : grid spacing dx, (dx,dy) or (dx,dy,dz) in 1D, 2D or 3D
        domain_size          : size of the physical domain Lx, (Lx,Ly), (Lx,Ly,Lz) in 1D, 2D or 3D
        cells                : number of cells nx or (nx, ny) or (nx, ny, nz) in 1, 2 and 3D.
        final_time           : final simulation time. Must be set if 'time_step' is not
        time_step            : simulation time step. Must be specified if 'final_time' is not
        interp_order         : order of the particle/mesh interpolation. Either 1, 2, 3 or 4 (default=1)
        layout               : layout of the physical quantities on the mesh (default = "yee")
        origin               : origin of the physical domain, (default (0,0,0) in 3D)
        refined_particle_nbr : number of refined particles for particle splitting ( TODO default hard-coded to 2)
        particle_pusher      : algo to push particles (default = "modifiedBoris")
        path                 : path for outputs (default : './')
        boundary_types       : type of boundary conditions (default is "periodic" for each direction)
        diag_export_format   : format of the output diagnostics (default= "ascii")

        """

        self.cells = kwargs['cells']
        self.dims = compute_dimension(self.cells)
        self.origin = kwargs['origin']
        self.path = kwargs['path']
        self.layout = kwargs['layout']
        self.particle_pusher = kwargs['particle_pusher']
        self.dl = kwargs['dl']
        self.time_step = kwargs['time_step']
        self.time_step_nbr = kwargs['time_step_nbr']
        self.interp_order = kwargs['interp_order']
        self.boundary_types = kwargs['boundary_types']
        self.refined_particle_nbr = kwargs['refined_particle_nbr']
        self.diag_export_format = kwargs['diag_export_format']

        self.levels_to_refine = []
        self.extent_ratio = []
        self.refinement_iterations = []

        self.diagnostics = []
        self.model = None

    def final_time(self):
        return self.time_step * self.time_step_nbr

    def simulation_domain(self):
        return [dl * n + self.origin[0] for dl, n in zip(self.dl, self.cells)]

    def within_simulation_duration(self, time_period):
        return time_period[0] >= 0 and time_period[1] < self.time_step_nbr

    def within_simulation_extent(self, extent):
        domain = self.simulation_domain()
        if len(extent) == 2:
            # 1D case
            return extent[0] >= domain[0] and extent[1] <= domain[1]
        else:
            raise NotImplementedError("Error: 2D and 3D not implemented yet")

    def refine(self, levels, extent_ratio, refinement_iterations):
        if len(levels) != len(refinement_iterations):
            raise ValueError("Error: levels and refinement_iterations must have the same length")
        if extent_ratio[0] >= extent_ratio[1]:
            raise ValueError("Error: extent_ratio[0] must be < extent_ratio[1]")
        if extent_ratio[0] < 0:
            raise ValueError("Error: extent_ratio[0] must be > 0")
        if extent_ratio[1] > 1:
            raise ValueError("Error: extent_ratio[1] must be < 1")

        self.levels_to_refine = levels
        self.extent_ratio = extent_ratio
        self.refinement_iterations = refinement_iterations


# ------------------------------------------------------------------------------

    def add_diagnostics(self, diag):
        if diag.name in self.diagnostics:
            raise ValueError("Error: diagnostics {} already registered".format(diag.name))

        # check if the duration of the diagnostics is valid, given the duration of the simulation
        if not self.within_simulation_duration(diag.iteration_interval()):
            raise RuntimeError("Error: invalid diagnostics duration")

        # check whether the spatial extent of the diagnostics is valid, given the domain size
        if diag.extent() is not None and not self.within_simulation_extent(diag.extent()):
            raise RuntimeError("Error: invalid diagnostics spatial extent")

        self.diagnostics.append(diag)


# ------------------------------------------------------------------------------

    def write_ini_file(self, filename="phare.ini"):

        config = configparser.ConfigParser()
        config.add_section('Simulation')

        dim_names = ('x', 'y', 'z')
        for d in np.arange(self.dims):
            config.set('Simulation', 'nbr_cells_'+dim_names[d], str(self.cells[d]))
            config.set('Simulation', 'd'+dim_names[d],str(self.dl[d]))
            config.set('Simulation', 'origin_'+dim_names[d], str(self.origin[d]))
            config.set('Simulation', 'boundary_condition_'+dim_names[d], self.boundary_types[d])

        params = {"section": "Simulation",
                  "layout": str(self.layout),
                  "interp_order": str(self.interp_order),
                  "time_step": str(self.time_step),
                  "particle_pusher": self.particle_pusher,
                  "refined_particle_nbr": self.refined_particle_nbr,
                  "time_step_nbr": str(self.time_step_nbr),
                  "diag_export_format": self.diag_export_format}

        add_dict_to_config(params, config)

        if len(self.extent_ratio) != 0:
            amr = {"section": "amr",
                   "min_ratio": self.extent_ratio[0],
                   "ratio_ratio": self.extent_ratio[1],
                   "refine_at_iteration": ", ".join([str(i) for i in self.refinement_iterations]),
                   "levels_to_refine": ", ".join([str(l) for l in self.levels_to_refine]),
                   "patch_to_refine": ", ".join([str(p) for p in [0]*len(self.levels_to_refine)])}

            add_dict_to_config(amr, config)

        for diag in self.diagnostics:
            add_dict_to_config(diag.to_dict(), config, section_key="name")

        add_dict_to_config(self.model.to_dict(), config, section_key="model")

        with open(os.path.join(self.path, filename), "w") as config_file:
            config.write(config_file)


# ------------------------------------------------------------------------------

    def set_model(self, model):
        self.model = model

# ------------------------------------------------------------------------------


def prepare_job(sim):
    """
    prepare a simulation for a run:
        - creates the simulation directory [simulation.path] if it does not exist yet
        - writes an INI file in the directory [simulation.path]
        - for each diagnostic registered in the simulation:
            - creates the diagnostic directory 'diag['path']' under [simulation.path]/
              if it does not exist yet
    """

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

