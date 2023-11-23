'''
This is a module for the Ecohydrology model that can be imported from the main driver/coupler notebook.
'''

import numpy as np
from landlab import RasterModelGrid 
from landlab.components import (Radiation, PotentialEvapotranspiration)
from generate_uniform_precip import PrecipitationDistribution
from vegetation_dynamics import Vegetation
from soil_moisture_dynamics import SoilMoisture


class EcoHyd:
    def __init__(self):
        #self.i = i
        #self.j = j

        #TODO might want to make the config file something that is passed in from the outside.
        self.config = {
            # assume Canicula starts on 1st of January and lasts for 100 days
            'canicula_start':1,
            'canicula_end':100, 
            # I pulled the following out of my behind, need to check with local climatology if this is sensible. 
            # Times are all in hours.
            'mean_interstorm_wet':5*24,
            'mean_storm_wet':2*24,
            'mean_interstorm_dry':20*24,
            'mean_storm_dry':2*24
        }

        # Represent current time in years (N.B. this is a float, so 0.5 would be mid-June of the first model year)
        self.current_time = 1/365  # Start from first day of Jan
        self.n=365 # length of single model loop in days

        # declaring few variables that will be used in the storm loop
        self.WS = 0.0  # Buffer for Water Stress
        self.Tg = 270  # Growing season in days

        #initialise some of the timeseries 

        self.P = np.zeros(self.n) #empty array for rainfall

        self.Time = np.zeros(self.n) #empty array to record timestamps at which calculations are made in main loop 
                        #(non-uniform time step length)


        #set up grid of size 53*53. This will result in 51*51 cells plus a rim of nodes around them (hence 53*53).
        #the inputs and outputs we need to pass all live on cells, not nodes. 
        self.mg = RasterModelGrid((53, 53), 5.)

        #let's try to add an idealised elevation profile to this grid.
        valley = np.zeros((53,53))

        def valleyfunc(x, y):
            e = 0.02*(x-25)**2 - 0.02*(y-25)**2 + 60
            return e

        for x in np.arange(0, 53, 1):
            for y in np.arange(0, 53, 1):
                valley[y][x] = valleyfunc(x,y)

        self.mg.add_field("topographic__elevation", valley, at="node", units="m", copy=True, clobber=True) 
        # The elevation field needs to have exactly 
        # this name in order for components to be 
        # able to work with it.

        #############################
        #generate precipitation data#
        #############################

        #'''
        #Commenting this out because it assigns things to the grid. Current implementation needs just an array though.
        # N.B. to make this (grid-based implementation of precipitation) run, find the 'generate_uniform_precip.py' file in your 
        # LandLab distribution, and add a 'clobber=True' statement to the self.grid.add_field() call in L185 
        # (this becomes more obvious if you run the cell and look at the error message).
        self.PD_D = PrecipitationDistribution(self.mg, mean_storm_duration=self.config['mean_storm_dry'], 
                                              mean_interstorm_duration=self.config['mean_interstorm_dry'],
                                              mean_storm_depth=0.5, total_t=40)

        self.PD_W = PrecipitationDistribution(self.mg, mean_storm_duration=self.config['mean_storm_wet'], 
                                              mean_interstorm_duration=self.config['mean_interstorm_wet'],
                                              mean_storm_depth=0.5, total_t=325)

        #if we choose this way of assigning precipitation, there should then be a `rainfall__flux` field on our grid, although I am not sure 
        #how to access it as it is added to the entire grid rather than individual nodes.
        #'''

        '''
        PD_D = PrecipitationDistribution(mean_storm_duration=config['mean_storm_dry'], mean_interstorm_duration=config['mean_interstorm_dry'],
                                        mean_storm_depth=0.5, total_t=40)

        PD_W = PrecipitationDistribution(mean_storm_duration=config['mean_storm_wet'], mean_interstorm_duration=config['mean_interstorm_wet'],
                                        mean_storm_depth=0.5, total_t=325)
        '''

        #-------------------------------#
        #instantiate radiation component#
        self.rad = Radiation(self.mg, current_time=self.current_time)
        self.rad.update()

        #--------------------------------------------------#
        #instantiate Potential Evapotranspiration Component#
        #TODO set the PET method to Priestley-Taylor once we have temperature inputs (daily Tmax, Tmin, Tavg)
        self.PET = PotentialEvapotranspiration(self.mg, current_time=self.current_time)
        self.PET.update() #running this initialises the output fields on the grid (which the soil moisture component needs)

        #-----------------------------------#
        #instantiate Soil Moisture Component#

        # pick some initial values for the fields that the soil moisture component requires to run - if we don't supply all of these fields, 
        # the initialisation throws an error!
        self.mg.at_cell['soil_moisture__initial_saturation_fraction'] = 0.75*np.ones(self.mg.number_of_cells) #75%
        self.mg.at_cell['rainfall__daily_depth'] = np.ones(self.mg.number_of_cells) #is 1mm of daily rainfall a valid assumption? How does it 
                                                                        #use this, does it multiply this value with the legth of 
                                                                        #the time step expressed in days?
        self.mg.at_cell['vegetation__cover_fraction'] = np.ones(self.mg.number_of_cells) #full cover everywhere
        self.mg.at_cell['vegetation__live_leaf_area_index'] = 2*np.ones(self.mg.number_of_cells) #not sure what a sensible value is here.
                                                                                    # from doc:
                                                                                    # one-sided green leaf area per unit ground surface area.
        self.mg.at_cell['vegetation__plant_functional_type'] = 0*(np.ones(self.mg.number_of_cells)).astype(int) #from doc:
                                                            # "classification of plants (int), grass=0, shrub=1, tree=2, "
                                                            #"bare=3, shrub_seedling=4, tree_seedling=5" - i.e., just let everything be 'grass'.

        self.SM = SoilMoisture(self.mg)
        self.SM.initialize()


        #Instantiate Vegetation Component

        #again we need to initialise some fields to get this to run
        self.mg.at_cell['surface__potential_evapotranspiration_30day_mean'] = self.mg.at_cell['surface__potential_evapotranspiration_rate'] 
                        # set the mean equal to the initial value (is this sensible?)
        self.mg.at_cell['surface__WSA_soilhealth'] = np.ones(self.mg.number_of_cells) 
        self.VEG = Vegetation(self.mg)


        # finally, we define a lower bound for WSA_soilhealth (probs just 1) and an upper bound 
        # (say e.g. 1.3). If WSA=True on a field, compute upperbound-WSA_soilhealth and increment 
        # WSA_soilhealth by a fixed percentage of that value
        # (e.g., half.). Do the reverse on fields where WSA=False.
        self.WSA_sh_lower = 1.
        self.WSA_sh_upper = 1.3


    def stepper(self, WSA_array):
        '''
        Run a one-year loop of the Ecohydrology model at a daily time step.
        '''

        #biomass = np.zeros((365, 101**2))

        #assume that '1' is for using WSA and '0' is for not using WSA. We also assume that WSA is only about 
        #cover cropping for now. This means that in the dry season, we would set the 
        #'vegetation__plant_functional_type' field to 3 (bare) in non-WSA plots, while it is 6 (Cover Crop) in WSA
        #plots. In both cases, the biomass would have to be reset to zero at the beginning and end
        #of the dry season/whenever the crop is changed. 

        #set functional type to 0 everywhere in growing season
        functype_growing = 0*np.ones(self.mg.number_of_cells).astype(int)

        #set functional type according to mask in non-growing season
        functype_nongrowing = 6*np.ones(WSA_array.shape).astype(int)
        functype_nongrowing[WSA_array == 0] = int(3)
        functype_nongrowing = functype_nongrowing.flatten()
        if len(functype_nongrowing) != self.mg.number_of_cells:
            print('grid size: ', self.mg.number_of_cells, 'wsa size:', len(functype_nongrowing))
            raise Exception('sorry, WSA array provided has wrong shape for the grid')

        for i in range(0, 365):
            # Update objects

            # Calculate Day of Year
            Julian = int(np.floor((self.current_time - np.floor(self.current_time)) * 365.0))
            print(Julian)
            print(self.current_time)

            # Generate seasonal storms
            # for Dry season
            if Julian < self.config['canicula_end'] or Julian > self.config['canicula_start']: 
                self.PD_D.update()
                self.P[i] = self.PD_D.storm_depth
            # Wet Season 
            else:
                self.PD_W.update()
                self.P[i] = self.PD_W.storm_depth

            # At the start of the canicula, harvest all fields and change the PFT to bare soil on 
            # non-WSA fields
            if Julian == self.config['canicula_start']:
                self.mg.at_cell['vegetation__plant_functional_type'] = functype_nongrowing
                #record the last state of the biomass before we overwrite by re-initialising ('harvesting')
                biomass = self.mg.at_cell['vegetation__live_biomass']
                self.VEG.initialize()

            # At the end of the canicula, harvest WSA fields and set non-WSA PFT back to grass
            if Julian == self.config['canicula_end']:
                self.mg.at_cell['vegetation__plant_functional_type'] = functype_growing
                #record soil moisture at end of canicula to see if WSA makes a difference
                SM_canic_end = self.mg.at_cell['soil_moisture__saturation_fraction']
                self.VEG.initialize()
                
            # calculate radiation for each field based on day of the year
            #self.rad.current_time = self.current_time
            self.rad.update()

            # Spatially distribute PET and its 30-day-mean (analogous to degree day)
            #mg["cell"]["surface__potential_evapotranspiration_rate"] = PET_[Julian]
            #mg["cell"]["surface__potential_evapotranspiration_30day_mean"] = EP30[Julian]

            # calculate PET for each field based on day of the year
            #TODO read in temperature data for each day and update PET with those
            self.PET.update()

            # Assign spatial rainfall data
            self.mg.at_cell["rainfall__daily_depth"] = self.P[i] * np.ones(self.mg.number_of_cells)

            # Update soil moisture component
            #TODO need to fudge this so WSA makes a difference. Add a term that increases resistance to evaporation as 
            #well as reducing root zone leakage in the presence of vegetation. 
            self.current_time = self.SM.update()

            # Decide whether its growing season or not (comment this out, think it is irrelevant as the 
            # canicula kind of is this)
            # if Julian != 364:
            #    if (Julian > 100 or Julian < 300):
            #        PET_threshold = 1
            #        # 1 corresponds to ETThresholdup (begin growing season)
            #    else:
            #        PET_threshold = 0
            #        # 0 corresponds to ETThresholddown (end growing season)

            # Update vegetation component
            self.VEG.update()

            # Update yearly cumulative water stress data
            self.WS += (self.mg["cell"]["vegetation__water_stress"]) # need multiply this by time step in days if dt!=1day

            # Record time (optional)
            self.Time[i] = self.current_time

            #print some outputs
            print('soil moisture sat.:', self.mg.at_cell['soil_moisture__saturation_fraction'])
            print('live biomass: ', self.mg.at_cell['vegetation__live_biomass'])
            print('ET: ', self.mg.at_cell['surface__potential_evapotranspiration_rate'])
            print('ET30: ', self.mg.at_cell['surface__potential_evapotranspiration_30day_mean'])
            print('PFT: ', self.mg.at_cell['vegetation__plant_functional_type'])

            #write to biomass
            #biomass[i, :] = self.mg.at_cell['vegetation__live_biomass']

        # update soil health parameter
        WSA_sh_mask = np.ones(WSA_array.shape)
        WSA_sh_mask[WSA_array == 0] = self.WSA_sh_lower
        WSA_sh_mask[WSA_array == 1] = self.WSA_sh_upper
        WSA_sh_mask = WSA_sh_mask.flatten()
        self.mg.at_cell['surface__WSA_soilhealth'] += (WSA_sh_mask - self.mg.at_cell['surface__WSA_soilhealth'])*0.5 


        return biomass, SM_canic_end

