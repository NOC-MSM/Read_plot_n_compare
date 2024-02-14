# Read_plot_n_compare

a set of python functions to read NEMO-MEDUSA outputs, plot and compare them.

## How to call these ?

call them in a python function or jupyter notebook, the following way :

``` python

## If you want to test in Jupyter notebook,
## you might need the first lot of packages to work with :
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
#import h5py

## Import the functions defined in the `read_plot_n_compare.py` file :
from read_plot_n_compare import runclass
from read_plot_n_compare import subplot_diff_list, subplot_diff, subplot_no_diff, subplot_timeseries
from read_plot_n_compare import subplot_comp_list, subplot_comp
from read_plot_n_compare import prep_cube_transect, prep_cube_extract_depth, prep_cube_vert_inv
from read_plot_n_compare import prep_cube_Atl_Pac, prep_cube_Atl_Pac_hov_from_monthly, regrid_ORCA, fix_nemo_cube
from read_plot_n_compare import prep_cube_surf, prep_cube_Obs, prep_cube, prep_cube_diff, prep_Masked_cube, read_cube
from read_plot_n_compare import subplot_proj_regional, subplot_proj_Orcagrid, plot_proj_Orcagrid
from read_plot_n_compare import bbox_extract_2Dcoords, nc_obs2D_varname, mat_obs2D_varname, mat_obs3D_varname

### !!! To Be Adapted !!!
### define the experiments outputs to use in your analyses :
### in `set_runs_input_to_plot_surf.py` :
###    -- You need to define the path and output file names to be used in the script --
###===================================
from set_runs_input_to_plot_surf import rundict_file ## , GOSI9dict_file ## you can def and add more if needed.
###===================================

### dont want the huge unit var list :
import warnings
warnings.filterwarnings('ignore')



yr = "YY"

### Ancillary files -- Needed by some functions :
reg_grid_file = "../MESH/Data_reg_grid/woa13_all_i00_01.nc"
mesh_mask     = "../MESH/eORCA1/NEMO42/mesh_mask.nc"
maskfile      = "../MESH/eORCA1/NEMO42/eORCA100_masks.nc"
#maskfile = "/noc/users/axy/Matlab/Macros/woa_mask100.mat"



#### Experiment name (used for plot legends and as
#### dictionary key to register the outputs path and names):
runlist = [ "cste",
            "clim_chl",
            "1_way",
            "2_ways",
            "coupled",
            "2_ways_ifr",
            "2_ways_ifr_dc"]

#### Title associated to experiments
ttl_list = ["cste",
            "clim_chl",
            "1_way",
            "2_ways",
            "coupled",
            "2_ways_ifr",
            "2_ways_ifr_dc"]

#### Path to experiment output files
rundict = {
   runlist[0] : "MED_CHL/CHL_cst",
   runlist[1] : "MED_CHL/CHL_clim",
   runlist[2] : "MED_CHL/CHL_1way",
   runlist[3] : "MED_CHL/CHL_2ways",
   runlist[4] : "MED_CHL/CHL_cpl",
   runlist[5] : "MED_CHL/CHL_2ways_ifr",
   runlist[6] : "MED_CHL/CHL_2ways_dc_ifr"
}


[rundict_ptrc,rundict_diad,rundict_grid] = rundict_file(yr, runlist, rundict)
#[GOSI9dict_ptrc,GOSI9dict_diad,GOSI9dict_grid] = GOSI9dict_file(yr, runlist, rundict)
#[rundict_sumptrc,rundict_sumdiad,rundict_sumgrid] = rundict_sum_file(yr)
### for example : get ptrc output of 1st exp with file_ptrc = rundict_ptrc[runlist[0]]


###############################
###############################
# init done

---------------------
## choose the expiriments you want to compare:
## Here we compare  the experiment 3 to 6 and 5
## the numbers are the one set in the dictionnary and lists above.

#run_exp_list=[1,3,2,4,0]
run_exp_list=[3,6,5]

## Some plot examples :
##------------------------------------

## Surface plots
subplot_diff("CHL", rundict_ptrc, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("PP", rundict_diad, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("XPAR3D", rundict_diad, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("ML_PP", rundict_diad, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("SUB_ML_PP", rundict_diad, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("mldr10_1", rundict_grid, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("qsr_oce", rundict_grid, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("tos", rundict_grid, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)
subplot_diff("DIN", rundict_ptrc, runlist, run_exp_list, Meth="Surface",region = "NorthOrtho", MASK_ZERO=True, SINGLE_PLOTS=False)

## Atl-pac transect
subplot_diff("CHL", rundict_ptrc, runlist, run_exp_list, Meth="Atl_Pac_Lat",NEW_regrid=False, EPIPEL=True, SINGLE_PLOTS=False)
subplot_diff("TPP3", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_Lat",NEW_regrid=False, EPIPEL=True, SINGLE_PLOTS=False)
subplot_diff("XPAR3D", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_Lat",NEW_regrid=False, EPIPEL=True, SINGLE_PLOTS=False)
subplot_diff("DIN", rundict_ptrc, runlist, run_exp_list, Meth="Atl_Pac_Lat",NEW_regrid=False, EPIPEL=True, SINGLE_PLOTS=False)
subplot_diff("thetao", rundict_grid, runlist, run_exp_list, Meth="Atl_Pac_Lat",NEW_regrid=False, EPIPEL=False, SINGLE_PLOTS=False)


## Hovmoeller plots
## on the Hovmoeller script
##======================================
#subplot_diff("CHL", rundict_ptrc, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("PP", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("XPAR3D", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("ML_PP", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("SUB_ML_PP", rundict_diad, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("mldr10_1", rundict_grid, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("qsr_oce", rundict_grid, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("tos", rundict_grid, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)
#subplot_diff("DIN", rundict_ptrc, runlist, run_exp_list, Meth="Atl_Pac_hov",NEW_regrid=False, SINGLE_PLOTS=False)

```

## More informations: 
the most used functions (like `subplot_diff` ) have a documentation of their own. 
You can se it with `sub_plot_diff.__doc__` ; Or see the head of the function in `read_plot_n_compare.py`
