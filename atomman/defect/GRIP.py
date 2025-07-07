from typing import Optional
import secrets

import numpy as np

from .. import System
from . import GrainBoundary
from . import interstitial_site_finder

from yabadaba.record import Record

class GRIP(Record):
    """
    Class for managing input settings and building grain boundary
    configurations according to the GRIP algorithm.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'grip'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'grip'

    ####################### Define Values and attributes #######################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """
        # Shift input settings
        self._add_value('float', 'shift_delta', defaultvalue=0.05,
                        modelpath='input-parameter.shift.delta',
                        description=' '.join([
                            'The spacing to use between possible shift samples.',
                            'Default value is 0.05. Note shifts should range',
                            'from 0 to 1.']))

        # Temperature input settings
        self._add_value('float', 'temperature_min', defaultvalue=0,
                        modelpath='input-parameter.temperature.min',
                        description=' '.join([
                            'The minimum temperature to include in the sampling.']))
        self._add_value('float', 'temperature_max', defaultvalue=2000,
                        modelpath='input-parameter.temperature.max',
                        description=' '.join([
                            'The maximum temperature to include in the sampling.']))
        self._add_value('float', 'temperature_delta', defaultvalue=100,
                        modelpath='input-parameter.temperature.delta',
                        description=' '.join([
                            'The spacing to use between possible temperature',
                            'samples.  Default value is 100.']))
        self._add_value('str', 'temperature_sample_style', defaultvalue='uniform',
                        modelpath='input-parameter.temperature.sample_style',
                        allowedvalues=['uniform'],
                        description=' '.join([
                            'The sampling style to use for the temperature.',
                            '"uniform" will uniformly sample from the range provided',
                            'and is the only currently supported style.']))
        
        # Sizemults input settings
        self._add_value('int', 'sizemults_max1', defaultvalue=3,
                        modelpath='input-parameter.sizemults.max1',
                        description=' '.join([
                            'The maximum size multiplier to use for the first in-plane',
                            'box vector. The selected sample will be between 1 and this',
                            'value inclusively. Default value is 3.']))
        self._add_value('int', 'sizemults_max2', defaultvalue=3,
                        modelpath='input-parameter.sizemults.max2',
                        description=' '.join([
                            'The maximum size multiplier to use for the second in-plane',
                            'box vector. The selected sample will be between 1 and this',
                            'value inclusively. Default value is 3.']))
        self._add_value('str', 'sizemults_sample_style',
                        modelpath='input-parameter.sizemults.sample_style',
                        defaultvalue='exponential',
                        allowedvalues=['uniform', 'exponential'],
                        description=' '.join([
                            'Indicates the sample style to use for the sizemults.',
                            '"exponential" (default) will place higher sample weights',
                            'on smaller sizemults. "uniform" will equally sample from',
                            'the sizemults range.']))
        
        # Runsteps input settings
        self._add_value('float', 'runsteps_chance', defaultvalue=0.95,
                        modelpath='input-parameter.runsteps.chance',
                        description=' '.join([
                            'A probability rate between 0 and 1 for if the MD relaxation',
                            'will be performed: a value of 0 means MD is never done and',
                            'a value of 1 means MD is always done.  Default value is',
                            '0.95 (5% chance of no MD).']))
        self._add_value('int', 'runsteps_min', defaultvalue=5000,
                        modelpath='input-parameter.runsteps.min',
                        description=' '.join([
                            'The minimum runsteps value to include in the sampling.',
                            'Default value is 5000.']))
        self._add_value('int', 'runsteps_max', defaultvalue=500000,
                        modelpath='input-parameter.runsteps.max',
                        description=' '.join([
                            'The maximum runsteps value to include in the sampling.',
                            'Default value is 500000.']))
        self._add_value('int', 'runsteps_delta', defaultvalue=1000,
                        modelpath='input-parameter.runsteps.delta',
                        description=' '.join([
                            'The spacing to use between possible runsteps samples.',
                            'Default value is 1000.']))
        self._add_value('str', 'runsteps_sample_style',
                        modelpath='input-parameter.runsteps.sample_style',
                        defaultvalue='exponential',
                        allowedvalues=['uniform', 'exponential'],
                        description=' '.join([
                            'The sampling style to use for the runsteps.',
                            '"uniform" will uniformly sample from the range provided.',
                            '"exponential" (default) will weight smaller runsteps as',
                            'being more likely.']))
        
        # Grain boundary builder settings
        self._add_value('int', 'maxmult', defaultvalue=10,
                        modelpath='input-parameter.builder.maxmult',
                        description=' '.join([
                            'The maximum size multiplier to use in searching for',
                            'correspondence between the in-plane vectors of both',
                            'grains.  Both grains will be searched up to this max,'
                            "so only one grain's multiplier is guaranteed to be",
                            'maxmult or less.  Default value is 10.']))
        self._add_value('float', 'minwidth', defaultvalue=30,
                        modelpath='input-parameter.builder.minwidth',
                        description=' '.join([
                            'The minimum width for both grains perpendicular to',
                            'the grain boundary.  Default value is 30.']))
        
        # Atom deletion settings
        self._add_value('float', 'delete_min', defaultvalue=0.0,
                        modelpath='input-parameter.delete.min',
                        description=' '.join([
                            'Minimum fraction of atoms to delete.  Default value is 0.0.']))
        self._add_value('float', 'delete_max', defaultvalue=1.0,
                        modelpath='input-parameter.delete.max',
                        description=' '.join([
                            'Maximum fraction of atoms to delete.  Default value is 1.0.']))

        # Atom perturbation settings
        self._add_value('float', 'perturb_width', defaultvalue=10,
                        modelpath='input-parameter.perturb.width',
                        description=' '.join([
                            'The width from the grain boundary in each grain within which',
                            'atoms will be perturbed.  Default value is 10.']))
        self._add_value('float', 'perturb_max1', defaultvalue=0.3,
                        modelpath='input-parameter.perturb.max1',
                        description=' '.join([
                            'The maximum atomic perturbations that will be applied to',
                            'atoms in grain 1.  Default value is 0.3.']))
        self._add_value('float', 'perturb_max2', defaultvalue=0.0,
                        modelpath='input-parameter.perturb.max2',
                        description=' '.join([
                            'The maximum atomic perturbations that will be applied to',
                            'atoms in grain 2.  Default value is 0.0.']))

        # Interstitial site settings
        self._add_value('int', 'interstitial_max_num', defaultvalue=2,
                        modelpath='input-parameter.interstitial.max_num',
                        description=' '.join([
                            "Specifies the max number of interstitial sites to fill.",
                            "Note that the actual max may be less if the number of",
                            "atoms or the number of interstitial sites at the grain",
                            "boundary are smaller.  Default value is 2."]))
        self._add_value('float', 'interstitial_width', defaultvalue=1.5,
                        modelpath='input-parameter.interstitial.width',
                        description=' '.join([
                            "Only interstitial sites within the interstitial_width",
                            "distance from the grain boundary will be considered for",
                            "filling, and only atoms within 2*interstitial_width from",
                            "the grain boundary will be considered for moving."]))
        self._add_value('str', 'interstitial_sample_style',
                        modelpath='input-parameter.interstitial.sample_style',
                        defaultvalue='uniform',
                        allowedvalues=['uniform', 'volume'],
                        description=' '.join([
                            "The sampling weight to use for selection of the interstitial",
                            "sites to fill.  If 'uniform' (default), all sites will have",
                            "equal likelihood of being filled.  If 'volume', then the",
                            "likelihood of sites being filled is weighted based on the",
                            "site's volume."]))

        # Generated calculation settings and metadata
        self._add_value('int', 'randomseed',
                        modelpath='calculation-parameter.randomseed',
                        description="The random number generator seed.")
        self._add_value('float', 'shift1',
                        modelpath='calculation-parameter.shift1',
                        description="The relative shift along the first in-plane box vector.")
        self._add_value('float', 'shift2',
                        modelpath='calculation-parameter.shift2',
                        description="The relative shift along the second in-plane box vector.")
        self._add_value('int', 'sizemult1',
                        modelpath='calculation-parameter.sizemult1',
                        description="The size multiplier along the first in-plane box vector.")
        self._add_value('int', 'sizemult2',
                        modelpath='calculation-parameter.sizemult2',
                        description="The size multiplier along the second in-plane box vector.")
        self._add_value('int', 'runsteps',
                        modelpath='calculation-parameter.runsteps',
                        description="The number of MD relaxation steps to perform.")
        self._add_value('float', 'temperature',
                        modelpath='calculation-parameter.temperature',
                        description="The temperature at which to perform the MD relaxation.")
        self._add_value('float', 'density',
                        modelpath='calculation-parameter.density',
                        description="The atomic density of the grain boundary after atom deletion.")
        self._add_value('int', 'ninterstitials',
                        modelpath='calculation-parameter.ninterstitials',
                        description="The number of interstitial positions filled.")

    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return 'grip-settings'

    def boundary(self,
                 grainboundary: GrainBoundary,
                 randomseed: Optional[int] = None,
                 decimals: int = 6,
                 verbose: bool = False,
                 **kwargs):
        
        # Update any calculation input settings if provided
        if len(kwargs) > 0:
            self.set_values(**kwargs)

        # Set the random number seed and create a generator
        if randomseed is None:
            randomseed = secrets.randbits(32)
        self.randomseed = randomseed
        rng = np.random.default_rng(seed=randomseed)

        # Select random input values
        self.__set_shifts(rng, verbose)
        self.__set_sizemults(rng, verbose)
        self.__set_runsteps(rng, verbose)
        self.__set_temperature(rng, verbose)
        
        # Build and modify the grain boundary system
        system, natoms1 = self.__build_boundary(grainboundary, verbose)
        system, natoms1 = self.__delete_atoms(system, natoms1, rng,
                                              grainboundary, decimals, verbose)
        self.__perturb_atoms(system, natoms1, rng, grainboundary, verbose)
        self.__interstitial_atoms(system, rng, grainboundary, verbose)

        return system, natoms1

    def __set_shifts(self,
                     rng: np.random.Generator,
                     verbose: bool):
        """
        Selects the two in-plane vector shifts to apply to the system using
        random uniform samples.
        """
        # Build array of possible shifts to sample from 
        possibleshifts = np.arange(0, 1, self.shift_delta)

        # Select gb shifts uniformly from possible values
        self.shift1 = rng.choice(possibleshifts)
        self.shift2 = rng.choice(possibleshifts)
        
        if verbose:
            print('shift1:', self.shift1)
            print('shift2:', self.shift2)

    def __set_sizemults(self,
                        rng: np.random.Generator,
                        verbose: bool):
        """
        Selects random sizemults values for the two in-plane box vectors.
        """
        # Build array of all possible sizemults
        sizemults1 = np.arange(1, self.sizemults_max1 + 1)
        sizemults2 = np.arange(1, self.sizemults_max2 + 1)

        if self.sizemults_sample_style == 'uniform':
            # Uniformly sample from allowed values
            self.sizemult1 = rng.choice(sizemults1)
            self.sizemult2 = rng.choice(sizemults2)

        elif self.sizemults_sample_style == 'exponential':
            # Place higher weights on smaller systems
            sizemultsweights1 = np.exp(-sizemults1) / np.sum(np.exp(-sizemults1))
            sizemultsweights2 = np.exp(-sizemults2) / np.sum(np.exp(-sizemults2))
        
            # Select samples using the weights
            self.sizemult1 = rng.choice(sizemults1, p=sizemultsweights1)
            self.sizemult2 = rng.choice(sizemults2, p=sizemultsweights2)
    
        else:
            raise ValueError('Unsupported sizemults_sample_style value')
        
        if verbose:
            print('sizemult1:', self.sizemult1)
            print('sizemult2:', self.sizemult2, flush=True)

    def __set_runsteps(self,
                       rng: np.random.Generator,
                       verbose: bool):
        """
        Select a random number of runsteps and temperature for performing the
        MD relaxation.
        """
        
        if rng.random() > self.runsteps_chance:
            # Skip MD run if probability check fails 
            self.runsteps = 0
        
        elif self.runsteps_sample_style == 'uniform':
            # Uniformly sample the number of MD steps
            self.runsteps = int(np.round(rng.choice(np.arange(self.runsteps_min, self.runsteps_max + 1, self.runsteps_delta))))
        
        elif self.runsteps_sample_style == 'exponential':
            # Exponentially scale the number of MD steps
            C = np.log(self.runsteps_max / self.runsteps_min)
            self.runsteps = int(np.round(self.runsteps_min * np.exp(C * rng.random()) / self.runsteps_delta) * self.runsteps_delta)

        else:
            raise ValueError('Unsupported runsteps_sample_style value')
        
        if verbose:
            print('runsteps:', self.runsteps, flush=True)
    
    def __set_temperature(self,
                          rng: np.random.Generator,
                          verbose: bool):
        """
        Select a random number of runsteps and temperature for performing the
        MD relaxation.
        """
        if self.runsteps == 0:
            self.temperature = 0.0

        elif self.temperature_sample_style == 'uniform':
            # Select a random uniform sample for the temperature
            self.temperature = np.round(rng.choice(np.arange(self.temperature_min, self.temperature_max + 1, self.temperature_delta)))

        else:
            raise ValueError('Unsupported temperature_sample_style value')
        
        if verbose:
            print('temperature:', self.temperature, flush=True)

    def __build_boundary(self,
                         gb: GrainBoundary,
                         verbose: bool):
        """
        Builds the grain boundary system using the builder and the selected
        sizemults and shifts.
        """
        mults1, mults2, strain = gb.identifymults(maxmult=self.maxmult,
                                                  minwidth=self.minwidth,
                                                  setvalues=False)
        
        # Multiply the in-plane sizemults by the GRIP values
        i1 = 0
        i2 = 2
        if gb.cutboxvector == 'a':
            i1 = 1
        elif gb.cutboxvector == 'c':
            i2 = 1
        mults1[i1] *= self.sizemult1
        mults2[i1] *= self.sizemult1
        mults1[i2] *= self.sizemult2
        mults2[i2] *= self.sizemult2

        if verbose:
            print('mults1:', mults1)
            print('mults2:', mults2)

        # Call super to generate the boundary system
        system, natoms1 = gb.boundary(mults1=mults1, mults2=mults2,
                                      freesurface=True, straintype='top',
                                      shift1=self.shift1, shift2=self.shift2,
                                      deleter=0.0)
        
        if verbose:
            print('system width:', system.box.vects[gb.cutindex, gb.cutindex])
            print('natoms (initial):', system.natoms, flush=True)

        return system, natoms1
        
    def __delete_atoms(self,
                       system: System,
                       natoms1: int,
                       rng: np.random.Generator,
                       gb: GrainBoundary,
                       decimals: int,
                       verbose: bool):
        """
        Randomly selects atoms for deletion from the grain boundary.
        """        
        # GB plane region ranges from min coord in top grain to that plus dlat
        dlat = gb.dlat(decimals) - 10**-decimals
        mincoord = system.atoms.pos[:natoms1, gb.cutindex].min()
        boundaryz = mincoord + dlat

        # Identify and count gb plane atoms
        inplane = np.where(system.atoms.pos[:natoms1, gb.cutindex] < boundaryz)[0]
        natomsplane = len(inplane)
        
        # Pick a random number of atoms to delete
        ndel = int(rng.integers(np.floor(natomsplane * (1 - self.delete_max)),
                                np.ceil(natomsplane * (1 - self.delete_min)),
                                endpoint=True))
        
        # Randomly select ndel atoms for deletion
        todelete = rng.choice(inplane, size=ndel, replace=False)
        keepindex = [x for i, x in enumerate(range(system.natoms)) if i not in todelete]
        
        if verbose:
            print('# atoms being deleted:', ndel)
            for index in todelete:
                print('atom deleted at', system.atoms.pos[index])

        # Delete the selected atoms
        newsystem = system.atoms_ix[keepindex]
        newnatoms1 = natoms1 - ndel
        
        # Compute the grain boundary atomic density
        self.density = (natomsplane - ndel) / natomsplane

        if verbose:
            print('natoms (now):', newsystem.natoms)
            print('gb density:', self.density, flush=True)

        return newsystem, newnatoms1

    def __perturb_atoms(self,
                        system: System,
                        natoms1: int,
                        rng: np.random.Generator,
                        gb: GrainBoundary,
                        verbose: bool):
        """
        Randomly perturb atoms near the grain boundary.
        """
        # Split atomic positions by grain
        pos1 = system.atoms.pos[:natoms1]
        pos2 = system.atoms.pos[natoms1:]

        # Perturb grain 1 atoms
        mincoord = pos1[:, gb.cutindex].min()
        boundary = mincoord + self.perturb_width
        inboundary = pos1[:, gb.cutindex] < boundary
        pos1[inboundary, :] += self.perturb_max1 * rng.random([np.sum(inboundary), 3])

        # Perturb grain 2 atoms
        maxcoord = pos2[:, gb.cutindex].max()
        boundary = maxcoord - self.perturb_width
        inboundary = pos2[:, gb.cutindex] > boundary
        pos2[inboundary, :] += self.perturb_max2 * rng.random([np.sum(inboundary), 3])

        # Join pos and update in system
        system.pos = np.vstack([pos1, pos2])

        if verbose:
            print('atoms perturbed', flush=True)

    def __interstitial_atoms(self,
                             system: System,
                             rng: np.random.Generator,
                             gb: GrainBoundary,
                             verbose: bool):
        """
        Randomly moves atoms near the grain boundary into interstitial sites.
        This identifies both interstitial sites and atoms near the grain
        boundary, then randomly selects a random number of atoms to move into
        randomly selected interstitial sites.
        """
        # Quick return if max is 0
        if self.interstitial_max_num <= 0:
            return
        
        # Create slice of the atomic system around the grain boundary
        search_width = 10 * self.interstitial_width
        in_search = ((system.atoms.pos[:, gb.cutindex] > -search_width) &
                     (system.atoms.pos[:, gb.cutindex] < search_width))
        search_system = system.atoms_ix[in_search]

        # Find interstitial sites from the search system
        allsites = interstitial_site_finder(search_system)

        # Filter out sites away from the grain boundary
        sites = []
        for site in allsites:
            if (site.pos[gb.cutindex] > -self.interstitial_width and 
                site.pos[gb.cutindex] <  self.interstitial_width):
                sites.append(site)
        all_site_ids = [i for i in range(len(sites))]

        if verbose:
            print('# total interstitial sites:', len(sites), flush=True)

        # Find ids of atoms near the grain boundary
        in_boundary = ((system.atoms.pos[:, gb.cutindex] > -2 * self.interstitial_width) &
                       (system.atoms.pos[:, gb.cutindex] < 2 *  self.interstitial_width))
        atom_ids = np.where(in_boundary)[0]

        if verbose:
            print('# atoms for interstitial shifts:', len(atom_ids))

        # Select max number to move based on max_num, num atoms and num sites
        max_num = min([self.interstitial_max_num, len(sites), len(atom_ids)])

        # Pick a random number of atoms to move from 0 to max_num
        nmove = int(rng.integers(0, max_num, endpoint=True))

        if verbose:
            print('# atoms being shifted:', nmove)

        # Shuffle atom_ids and only select the first nmove
        np.random.shuffle(atom_ids)
        atom_ids = atom_ids[:nmove]

        if self.interstitial_sample_style == 'uniform':
            # Randomly select interstitial sites to fill with no weights
            site_ids = rng.choice(all_site_ids, size=nmove, replace=False)
        
        elif self.interstitial_sample_style == 'volume':
            # Randomly select interstitial sites to fill weighted towards large volumes
            volumes = np.array([site.volume for site in sites])
            weights = volumes / np.sum(volumes)
            site_ids = rng.choice(all_site_ids, size=nmove, replace=False, weights=weights)

        # Move atoms to the interstitial sites
        for atom_id, site_id in zip(atom_ids, site_ids):
            if verbose:
                print('atom at', system.atoms.pos[atom_id])
                print('moved to', sites[site_id].pos, flush=True)
            system.atoms.pos[atom_id] = sites[site_id].pos

        self.ninterstitials = nmove
        