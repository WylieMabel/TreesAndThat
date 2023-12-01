The contents of this folder are:
 - modelDriver.ipynb, which is the main driver notebook from which the full coupled model should be run.
 - this notebook calls modelScript.py, which sets up the coupling between netlogo and python and then calls the social
   system and ecohydrology models in order to timestep the model and write outputs
 - the social system model that is called from modelScript.py is coded up in modelv3.nlogo, while the ecohydrological model 
   time stepper that is called from modelScript.py lives in ecohydr_mod.py.
 - ecohydr_mod.py itself has dependencies, namely the landlab components we modified. These are soil_moisture_dynamics.py,
   vegetation_dynamics.py and generate_uniform_precip.py (in the last one we just had to fix a bug, no actual science here).
 - new_temp_data.csv is the temperature data the ecohydrological model needs as an input. It is read in from modelScript.py.
 - sampleModelOutput is an example output from the coupled model, but does not contain all of our runs as that would have 
   been too much data
 - test_ecohydro_module.ipynb is a notebook we used to run and explore the ecohydrological model on its own. Note 
   that to get output for multiple years, you have to run the cell that calls Ecohyd_model.stepper() multiple times.
   (running this once triggers a 365 daily time steps). 

To run the model from the driver, you need to be in a Python environment that has the Landlab, Pynetlogo and multiprocessing
libraries installed (as well as all the default stuff such as numpy, time etc.).
The full repository can be found on https://github.com/WylieMabel/TreesAndThat/tree/main 
