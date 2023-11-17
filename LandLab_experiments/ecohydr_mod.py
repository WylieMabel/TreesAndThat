'''
This is a module for the Ecohydrology model that can be imported from the main driver/coupler notebook.
'''

import numpy as np
from landlab import RasterModelGrid 
from landlab.components import (Vegetation, SoilMoisture, Radiation, PotentialEvapotranspiration,
                                PrecipitationDistribution)


class EcoHyd:
    def __init__(self, WSA_array):
        #self.i = i
        #self.j = j

        self.WSA_array = WSA_array 
        #TODO still need to write code that converts WSA array into inputs for components.

        #TODO might want to make the config file something that is passed in from the outside.
        self.config = {
            # assume Canicula starts on 1st of July and lasts for 40 days
            'canicula_start':182,
            'canicula_end':223, 
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


        #set up grid of size 103*103. This will result in 101*101 cells plus a rim of nodes around them (hence 103*103).
        #the inputs and outputs we need to pass all live on cells, not nodes. 
        self.mg = RasterModelGrid((103, 103), 5.)

        #let's try to add an idealised elevation profile to this grid.
        valley = np.zeros((103,103))

        def valleyfunc(x, y):
            e = 0.02*(x-51)**2 - 0.02*(y-51)**2 + 60
            return e

        for x in np.arange(0, 103, 1):
            for y in np.arange(0, 103, 1):
                valley[y][x] = valleyfunc(x,y)

        self.mg.add_field("topographic__elevation", valley, at="node", units="m", copy=True, clobber=True) # The elevation field needs to have exactly 
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
        self.PD_D = PrecipitationDistribution(self.mg, mean_storm_duration=self.config['mean_storm_dry'], mean_interstorm_duration=self.config['mean_interstorm_dry'],
                                        mean_storm_depth=0.5, total_t=40)

        self.PD_W = PrecipitationDistribution(self.mg, mean_storm_duration=self.config['mean_storm_wet'], mean_interstorm_duration=self.config['mean_interstorm_wet'],
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

        #################################
        #instantiate radiation component#
        self.rad = Radiation(self.mg, current_time=self.current_time)
        self.rad.update()

        ####################################################
        #instantiate Potential Evapotranspiration Component#
        self.PET = PotentialEvapotranspiration(self.mg, current_time=self.current_time)
        self.PET.update() #running this initialises the output fields on the grid (which the soil moisture component needs)

        #####################################
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
        self.VEG = Vegetation(self.mg)



    def stepper(self):
        '''
        Run a one-year loop of the Ecohydrology model at a daily time step.
        '''

        for i in range(0, 365):
            # Update objects

            # Calculate Day of Year (DOY)
            Julian = int(np.floor((current_time - np.floor(current_time)) * 365.0))

            # Generate seasonal storms
            # for Dry season
            if Julian < self.config['canicula_end'] or Julian > self.config['canicula_start']: 
                self.PD_D.update()
                self.P[i] = self.PD_D.storm_depth
            # Wet Season 
            else:
                self.PD_W.update()
                self.P[i] = self.PD_W.storm_depth

            # calculate radiation for each field based on day of the year
            self.rad.current_time = self.current_time
            self.rad.update()

            # Spatially distribute PET and its 30-day-mean (analogous to degree day)
            #mg["cell"]["surface__potential_evapotranspiration_rate"] = PET_[Julian]
            #mg["cell"]["surface__potential_evapotranspiration_30day_mean"] = EP30[Julian]

            # calculate PET for each field based on day of the year
            #PET.current_time = current_time
            self.PET.update()

            # Assign spatial rainfall data
            self.mg.at_cell["rainfall__daily_depth"] = self.P[i] * np.ones(self.mg.number_of_cells)

            # Update soil moisture component
            current_time = self.SM.update()

            # Decide whether its growing season or not
            if Julian != 364:
                if (Julian > 100 or Julian < 300):
                    PET_threshold = 1
                    # 1 corresponds to ETThresholdup (begin growing season)
                else:
                    PET_threshold = 0
                    # 0 corresponds to ETThresholddown (end growing season)

            # Update vegetation component
            VEG = Vegetation(self.mg, PETthreshold_switch=PET_threshold, Tb=self.Tb[i], Tr=self.Tr[i]) #need to re-instatiate in order to 
                                                                                        #be able to set threshold switch;
                                                                                        #need to make sure this does not overwrite 
                                                                                        #vegetation fields from previous time steps.
            VEG.update()

            # Update yearly cumulative water stress data
            WS += (self.mg["cell"]["vegetation__water_stress"]) # need multiply this by time step in days if dt!=1day

            # Record time (optional)
            self.Time[i] = current_time

