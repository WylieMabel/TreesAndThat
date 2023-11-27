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
    def __init__(self, config, init_min_T, init_max_T, init_avg_T):

        self.config = config

        self.canicula_length = self.config['canicula_end']-self.config['canicula_start']

        # Represent current time in years (N.B. this is a float, so 0.5 would be mid-June of the first model year)
        self.current_time = 0  # Start from first day of Jan
        self.n=365 # length of single model loop in days

        # declaring few variables that will be used in the storm loop
        self.WS = 0.0  # Buffer for Water Stress
        self.Tg = 270  # Growing season in days

        #initialise some of the timeseries 

        self.P = np.zeros(self.n) #empty array for rainfall

        self.Time = [] #empty list to record timestamps at which calculations are made in main loop 

        #set up grid of size 53*53. This will result in 51*51 cells plus a rim of nodes around them (hence 53*53).
        #the inputs and outputs we need to pass all live on cells, not nodes. 
        #We define the side length of grid cells to be 70m - this corresponds to an average farm being about 1.5 
        #hectares and consisting of 3.25 fields. 
        self.mg = RasterModelGrid((53, 53), 70.)

        #let's try to add an idealised elevation profile to this grid.
        valley = np.zeros((53,53))

        def valleyfunc(x, y):
            e = 0.001*(x-25)**2 - 0.001*(y-25)**2 + 60
            return e

        for x in np.arange(0, 53, 1):
            for y in np.arange(0, 53, 1):
                valley[y][x] = valleyfunc(x,y)

        #valley = np.zeros((53,53))

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
                                              mean_storm_depth=self.config['mean_raindpth_dry'], 
                                              total_t=self.canicula_length*24)

        self.PD_W = PrecipitationDistribution(self.mg, mean_storm_duration=self.config['mean_storm_wet'], 
                                              mean_interstorm_duration=self.config['mean_interstorm_wet'],
                                              mean_storm_depth=self.config['mean_raindpth_wet'], 
                                              total_t=(365-self.canicula_length)*24)

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
        self.Tmin = init_min_T
        self.Tmax = init_max_T
        self.Tavg = init_avg_T
        self.PET = PotentialEvapotranspiration(self.mg, method='PriestleyTaylor', current_time=self.current_time,
                                               Tmin = self.Tmin, Tmax = self.Tmax, Tavg = self.Tavg, latitude=10.)
        self.PET.update() 
        #running this initialises the output fields on the grid (which the soil moisture component needs)

        #-----------------------------------#
        #instantiate Soil Moisture Component#

        # pick some initial values for the fields that the soil moisture component requires to run - if we don't supply all of these fields, 
        # the initialisation throws an error!
        self.mg.at_cell['soil_moisture__initial_saturation_fraction'] = 0.75*np.ones(self.mg.number_of_cells) #75%
        self.mg.at_cell['rainfall__daily_depth'] = np.ones(self.mg.number_of_cells) #is 1mm of daily rainfall a valid assumption? How does it 
                                                                        #use this, does it multiply this value with the legth of 
                                                                        #the time step expressed in days?
        self.mg.at_cell['vegetation__cover_fraction'] = 0.5*np.ones(self.mg.number_of_cells) #full cover everywhere
        self.mg.at_cell['vegetation__live_leaf_area_index'] = np.ones(self.mg.number_of_cells) #not sure what a sensible value is here.
                                                                                    # from doc:
                                                                                    # one-sided green leaf area per unit ground surface area.
        self.mg.at_cell['vegetation__plant_functional_type'] = 0*(np.ones(self.mg.number_of_cells)).astype(int) #from doc:
                                                            # "classification of plants (int), grass=0, shrub=1, tree=2, "
                                                            #"bare=3, shrub_seedling=4, tree_seedling=5" - i.e., just let everything be 'grass'.

        self.SM = SoilMoisture(self.mg)
        self.SM.initialize()

        #--------------------------------#
        #Instantiate Vegetation Component#

        #again we need to initialise some fields to get this to run
        self.mg.at_cell['surface__potential_evapotranspiration_30day_mean'] = self.mg.at_cell['surface__potential_evapotranspiration_rate'] 
                        # set the mean equal to the initial value (is this sensible?)
        self.mg.at_cell['surface__WSA_soilhealth'] = np.ones(self.mg.number_of_cells) 
        # decrease ET threshold from the default value of 3.8 bc that meant farmers on N-facing slopes had huge losses
        self.VEG = Vegetation(self.mg, PETthreshold_switch=1, ETthreshold_up=3.)


        # finally, we define a lower bound for WSA_soilhealth (probs just 1) and an upper bound 
        # (say e.g. 1.3). If WSA=True on a field, compute upperbound-WSA_soilhealth and increment 
        # WSA_soilhealth by a fixed percentage of that value
        # (e.g., half.). Do the reverse on fields where WSA=False.
        self.WSA_sh_lower = 1.
        self.WSA_sh_upper = 1.3

        #-------------------------------#
        # initialise output time series #
        self.WSA_SM_tseries = [np.mean(self.mg.at_cell['soil_moisture__saturation_fraction'])]
        self.noWSA_SM_tseries = [np.mean(self.mg.at_cell['soil_moisture__saturation_fraction'])]
        self.WSA_biomass_tseries = [np.mean(self.mg.at_cell['vegetation__live_biomass'])]
        self.noWSA_biomass_tseries = [np.mean(self.mg.at_cell['vegetation__live_biomass'])]
        self.ET30_tseries = [np.mean(self.mg.at_cell['surface__potential_evapotranspiration_30day_mean'])]
        self.rain_tseries = [np.mean(self.mg.at_cell['rainfall__daily_depth'])]
        

    #--------------#
    # time stepper #
    #--------------#

    def stepper(self, WSA_array, avg_temp, maximum_temp, minimum_temp):
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
        
        #generate precipitation time series
        PD_raw = self.PD_D.get_storm_time_series()
        PW_raw = self.PD_W.get_storm_time_series()

        #put data into useful format
        self.P = np.zeros(365)
        # Iterate over each rainfall event and update the daily precipitation
        for i in range(len(PD_raw)):
            # Distribute the intensity over the corresponding days
            start_day = int(PD_raw[i][0]) // 24
            end_day = int(PD_raw[i][1]) // 24
            if end_day <= self.canicula_length:
                self.P[start_day:end_day+1] = PD_raw[i][2]*24

        for i in range(len(PW_raw)):
            # Distribute the intensity over the corresponding days
            start_day = int(PW_raw[i][0]) // 24 + self.canicula_length
            end_day = int(PW_raw[i][1]) // 24 + self.canicula_length
            if end_day <= 365:
                self.P[start_day:end_day+1] = PW_raw[i][2]*24  

        print(self.P)      


        for i in range(0, 365):
            # Update objects

            # Calculate Day of Year
            Julian = int(np.floor((self.current_time - np.floor(self.current_time)) * 365.0))
            #print(Julian)
            #print(self.current_time)

            # At the start of the canicula, harvest all fields and change the PFT to bare soil on 
            # non-WSA fields and cover crop on WSA fields
            if Julian == self.config['canicula_start_expected']:
                self.mg.at_cell['vegetation__plant_functional_type'] = functype_nongrowing
                #need to re-initialize SM component for it to recognise new PFT
                self.SM.initialize()
                self.VEG.initialize(Blive_init=10.0)

            # At the end of the canicula, harvest WSA fields and set all PFT back to grass
            if Julian == self.config['canicula_end_expected']:
                self.mg.at_cell['vegetation__plant_functional_type'] = functype_growing
                print(self.mg.at_cell['vegetation__plant_functional_type'])
                #record soil moisture at end of canicula to see if WSA makes a difference
                SM_canic_end = self.mg.at_cell['soil_moisture__saturation_fraction'].copy()
                self.SM.initialize()
                self.VEG.initialize(Blive_init=10.0)

            # do first harvest at the end of the first maize crop cycle (100 days)
            #if Julian == self.config['canicula_end_expected'] + 100:
            #    biomass = self.mg.at_cell['vegetation__live_biomass'].copy()
            #    self.VEG.initialize(Blive_init=10.0)
                
            # calculate radiation for each field based on day of the year
            self.rad.update()

            # calculate PET for each field based on day of the year
            self.PET.Tmin = minimum_temp[i]
            self.PET.Tmax = maximum_temp[i]
            self.PET.Tavg = avg_temp[i]
            self.PET.update()

            # Assign spatial rainfall data
            self.mg.at_cell["rainfall__daily_depth"] = self.P[i] * np.ones(self.mg.number_of_cells)

            # Update soil moisture component
            self.current_time = self.SM.update()

            # Update vegetation component
            self.VEG.update()

            # Update yearly cumulative water stress data
            self.WS += (self.mg["cell"]["vegetation__water_stress"]) # need multiply this by time step in days if dt!=1day

            # Record time 
            self.Time.append(self.current_time)

            #horrific hack to get around time stepping bug and still make sure PET and radiation know about current time
            if self.Time[i-1] < self.Time[i]:
                self.PET.current_time = self.current_time
                self.rad.current_time = self.current_time

            
            #write time series output for soil moisture and biomass
            self.WSA_SM_tseries.append(np.mean(
                self.mg.at_cell['soil_moisture__saturation_fraction'][WSA_array.flatten() == 1]))
            self.noWSA_SM_tseries.append(np.mean(
                self.mg.at_cell['soil_moisture__saturation_fraction'][WSA_array.flatten() == 0]))
            self.WSA_biomass_tseries.append(np.mean(
                self.mg.at_cell['vegetation__live_biomass'][WSA_array.flatten() == 1]))
            self.noWSA_biomass_tseries.append(np.mean(
                self.mg.at_cell['vegetation__live_biomass'][WSA_array.flatten() == 0]))
            self.ET30_tseries.append(np.mean(self.mg.at_cell['surface__potential_evapotranspiration_30day_mean']))
            self.rain_tseries.append(np.mean(self.mg.at_cell['rainfall__daily_depth']))
            
        # update soil health parameter at the end of the year
        WSA_sh_mask = np.ones(WSA_array.shape)
        WSA_sh_mask[WSA_array == 0] = self.WSA_sh_lower
        WSA_sh_mask[WSA_array == 1] = self.WSA_sh_upper
        WSA_sh_mask = WSA_sh_mask.flatten()
        self.mg.at_cell['surface__WSA_soilhealth'] += (WSA_sh_mask - self.mg.at_cell['surface__WSA_soilhealth'])*0.5 

        # assume second crop cycle just always ends at the end of the year 
        biomass = self.mg.at_cell['vegetation__live_biomass'].copy()


        return biomass, SM_canic_end

