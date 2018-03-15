
from pharein.globals import objects


class UniformModel(object):

    def __init__(self, b=(1., 0., 0.), e=(0., 0., 0.)):
        self.model = {"model": "model", "model_name": "uniform"}
        if "Simulation" not in objects:
            raise RuntimeError("A simulation must be declared before a model")

        if "Model" in objects:
            raise RuntimeError("A model is already created")

        else:
            objects["Model"] = self
            simulation = objects["Simulation"]
            simulation.set_model(self)

        if len(b) != 3 or not isinstance(b, tuple) or not isinstance(b, list):
            ValueError("invalid B")
        if len(e) != 3 or not isinstance(e, tuple) or not isinstance(e, list):
            ValueError("invalid E")

        self.model.update({"bx": b[0],
                           "by": b[1],
                           "bz": b[2],
                           "ex": e[0],
                           "ey": e[1],
                           "ez": e[2]})



# ------------------------------------------------------------------------------

    def nbr_species(self):
        """
        returns the number of species currently registered in the model
        """
        keys = self.model.keys()
        nbr = 0
        for k in keys:
            if k.startswith('species'):
                nbr += 1
        return nbr

#------------------------------------------------------------------------------

    def add_species(self, name,
                    charge=1,
                    mass=1,
                    nbr_part_per_cell=100,
                    density=1.,
                    vbulk=(0., 0., 0.),
                    beta=1.0,
                    anisotropy=1.):
        """
        add a species to the current model

        add_species(name,charge=1, mass=1, nbrPartCell=100, density=1, vbulk=(0,0,0), beta=1, anisotropy=1)

        Parameters:
        -----------
        name        : name of the species, str

        Optional Parameters:
        -------------------
        charge      : charge of the species particles, float (default = 1.)
        nbrPartCell : number of particles per cell, int (default = 100)
        density     : particle density, float (default = 1.)
        vbulk       : bulk velocity, tuple of size 3  (default = (0,0,0))
        beta        : beta of the species, float (default = 1)
        anisotropy  : Pperp/Ppara of the species, float (default = 1)
        """

        idx = str(self.nbr_species())

        new_species = {"speciesName"+idx: name,
                       "charge"+idx: charge,
                       "mass"+idx: mass,
                       "density"+idx: density,
                       "vx"+idx: vbulk[0],
                       "vy" + idx: vbulk[1],
                       "vz" + idx: vbulk[2],
                       "beta"+idx: beta,
                       "anisotropy"+idx: anisotropy,
                       "nbrParticlesPerCell"+idx: nbr_part_per_cell}

        keys = self.model.keys()
        if "speciesName"+idx in keys:
            raise ValueError("species already registered")

        self.model.update(new_species)
#------------------------------------------------------------------------------

    def to_dict(self):
        self.model['nbr_ion_populations'] = self.nbr_species()
        return self.model
