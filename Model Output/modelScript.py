import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pynetlogo
import numpy as np
import sys
import datetime
sys.path.append('../')

from ecohydr_mod import EcoHyd

def get_yearly_temp(csv_path, num_years):
    df = pd.read_csv(csv_path)
    df.dt = pd.to_datetime(df.dt)
    avg_temp_per_year = []
    max_temp_per_year = []
    min_temp_per_year = []
    for year in range(0, num_years):
        earliest_yr = df.iloc[0, 0].year
        year = year + earliest_yr
        if len(df.loc[(df.dt.dt.year == year), :]) > 365:
            avg_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'AverageTemperature'].iloc[:-1].tolist())
            max_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'MaxTemperature'].iloc[:-1].tolist())
            min_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'MinTemperature'].iloc[:-1].tolist())
        else:
            avg_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'AverageTemperature'].tolist())
            max_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'MaxTemperature'].tolist())
            min_temp_per_year.append(df.loc[(df.dt.dt.year == year), 'MinTemperature'].tolist())
    
    return avg_temp_per_year, max_temp_per_year, min_temp_per_year

def setUpNetLogoModel(leadFarmers, desperation, jealousy, grace):
    # think this is for the GUI idk?
    sns.set_style("white")
    sns.set_context("talk")

    # starts a NetLogo link, point the netlogo home path to where you have it installed on your machine
    netlogo = pynetlogo.NetLogoLink(
        gui=False,
        netlogo_home="/Volumes/NetLogo 6.3.0/NetLogo 6.3.0"
    )

    # loads a .nlogo model from provided path
    netlogo.load_model("./modelv3.nlogo")

    # runs the model setup command
    netlogo.command("setup")

     # sets globals
    globals = "update-globals " + str(leadFarmers) + " " + str(desperation) + " " + str(jealousy) + " " + str(grace)
    netlogo.command(globals)
    

    return netlogo

def reportsToDataFrame(netlogo):
    # gets field attributes and puts it in a data frame
    fieldAttributes = netlogo.report("get-info")
    sorted_list = sorted(fieldAttributes, key=lambda x: (-x[1],x[0]))
    fieldData = pd.DataFrame(columns=["who", "xcor","ycor","owner-id","implements-WSA", "owner-knows-WSA", "yield"], data=sorted_list)
    return fieldData

def convertWSAToNPArray(data):
    # sets bool into correct format to pass to hydrology model
    emptynp = np.empty((51,51))
    for i in range(0,51):
        for j in range(0,51):
            ycorFromIndex = i * -1 + 25
            xcorFromIndex = j - 25
            emptynp[i][j]= data.loc[(data['ycor'] == ycorFromIndex) & (data['xcor'] == xcorFromIndex)]["implements-WSA"].iloc[0]
    return emptynp

def convertHydrologyToDF(hydrologyArray, data):
    # method converts the hydrology model output into a pandas dataframer
    # 2 is the placeholder to multiply the biomass to get yield
    hydrologyData = data.copy()
    hydrologyData["yield"] = hydrologyArray.reshape((2601,1))
    return hydrologyData

def fullModelRun(climate, leadFarmers, social, input_csv_path, no_of_years):
    # sets up model
    netlogo = setUpNetLogoModel(leadFarmers, social[0], social[1], social[2])

    # this will record all field attributes throughout the simulation - need to add year index to differentiate
    baseFieldData = reportsToDataFrame(netlogo)

    returnedData = baseFieldData.copy()
    
    baseFieldData["Year"] = 0

    WSA_records = []

    
    #get input temperature data
    avg, maxi, mini = get_yearly_temp(input_csv_path, no_of_years)
    avg = np.array(avg)
    maxi = np.array(maxi)
    mini = np.array(mini)

    Ecohyd_model = EcoHyd(climate, 20, 26, 23)

    #--------------------------------------------#
    # let hydrology model spin up for five years #
    for i in range(0,5):
        #just use same initial WSA array for each year
        WSA_array = convertWSAToNPArray(returnedData)
        _,_ = Ecohyd_model.stepper(WSA_array, avg[0]+climate['tempshift'], maxi[0]+climate['tempshift'], mini[0]+climate['tempshift'])

    #---------------------------------#
    # actual coupled model loop whooo #
    for year in range(0, no_of_years):
        print("year:", year)

        # converts the usingWSA bool for each field into an NP array
        WSA_array = convertWSAToNPArray(returnedData)
        WSA_records.append([WSA_array])

        biomass_harvest, SM_canic_end = Ecohyd_model.stepper(WSA_array, avg[year]+climate['tempshift'], maxi[year]+climate['tempshift'], mini[year]+climate['tempshift'])

        #record outputs for yearly rainfall
        cum_rainfall = np.cumsum(Ecohyd_model.rain_tseries[(year)*365:(year+1)*365])[-1]

        returnedData.sort_values(by=['ycor', 'xcor'], ascending=[False, True], inplace=True)

        # converts the updated yields to the dataframe
        hydrologyData = convertHydrologyToDF(biomass_harvest,returnedData)
        
        # writes this new yield information to the netlogo implementation
        netlogo.write_NetLogo_attriblist(hydrologyData, "field")

        # runs one step of social model
        netlogo.command("farming-year")

        # converts field data to a df
        returnedData = reportsToDataFrame(netlogo)

        dataToRecord = returnedData.copy()
        
        dataToRecord["Year"] = year + 1

        dataToRecord["TotalYearRainfall"] = cum_rainfall

        # adds this years results to the dataframe
        baseFieldData = pd.concat([baseFieldData, dataToRecord], ignore_index=True)

    summarisedData = baseFieldData.groupby(["owner-id","Year"]).agg({'xcor':'mean','ycor':'mean','implements-WSA':'mean','owner-knows-WSA':'mean', 'yield':'sum', 'who':'count'}).reset_index()
    summarisedData.rename(columns={'owner-id':'FarmerID', 'xcor':'MeanXCor', 'ycor':'MeanYCor', 'implements-WSA':'ImplementingWSA','owner-knows-WSA':'KnowsWSA','yield':'TotalYield','who':'NumberofFields'})

    summarisedData["LeadFarmers"] = leadFarmers
    summarisedData["SocialScenario"] = social[3]
    if climate['mean_raindpth_wet'] == 10:
        summarisedData["ClimateScenario"] = "Current Climate"
    else:
        summarisedData["ClimateScenario"] = "Warm Climate"
    summarisedData["UniqueID"] = str(datetime.datetime.now())
    
    # this writes to a csv
    summarisedData.to_csv(path_or_buf="sampleModelOutput", mode = "a", index=False, header = True)