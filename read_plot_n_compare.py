
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma


### dont want the huge unit var list :
#import warnings
#warnings.filterwarnings('ignore')
#
#%matplotlib inline


'''
def def_var_n_dic(YY) :
	################
	################
	## def file and var : (default have to be updated)
    yr = "YY"
    ### Ancillary files -- Needed by some functions : 
    reg_grid_file = "../MESH/Data_reg_grid/woa13_all_i00_01.nc" 
    mesh_mask     = "../MESH/eORCA1/NEMO42/mesh_mask.nc"
    try:
        ## running local
        meshfile  = "../MESH/eORCA1/NEMO42/mesh_mask.nc"
    except:
        ## on NOC system 
        meshfile  = "/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
    maskfile      = "../MESH/eORCA1/NEMO42/eORCA100_masks.nc"
    #maskfile = "/noc/users/axy/Matlab/Macros/woa_mask100.mat"


    #### #experiment name :
    runlist = [ "Init-cond",
                "No-olivine",
                "surf",
                "100m_max",
                "200m_max"]

    #### Title associated to experiments
    ttl_list = ["Init-cond",
                "No-olivine",
                "surf",
                "100m_max",
                "200m_max"]

    #### Path to experiment output files
    rundict = {
        runlist[0] : "Out",
        runlist[1] : "Out",
        runlist[2] : "Out",
        runlist[3] : "Out",
        runlist[4] : "Out"
    }

    #### Somme Var lists :
    ##-----------------------------
    #Inv_2d_var = [ "RIV_ALK",
    #            "TOT_SHALK",
    #            "CO2FLUX"]

    Inv_2d_var = ["CO2FLUX"]

    Inv_2d_phys = [ "somxl010",
            "tos",
            "sos"]

    Inv_2d_ZMP_var = ["FRAGDF1",
               "FRAGDF2",
               "FRAGD",
               "GMPDF1",
               "GMPDF2",
               "GMPD"]

    Surf_3d_var = ["ZME",
               "CHL",
              "ZMI",
              "PHN",
              "PHD",
              "ALK",
              "DIC",
              "DIN",
              "FER",
              "SIL"]

    Surf_3d_zmp_var = ["ZMP"]

    Inv_3d_var =  ["ZME_E3T",
              "ZMI_E3T",
              "PHN_E3T",
              "PHD_E3T",
              "PHYTO",
              "ZOO",
              "PLKT",
              "ALK_E3T",
              "DIC_E3T",
              "DIN_E3T",
              "FER_E3T",
              "SIL_E3T"]

    Inv_3d_zmp_var = ["ZMP_E3T"]

    Transect_3d_var = ["TPP3",
                  "DETFLUX3"]


    Obs_surf_compare_var = ["CHL",
                        "ALK",
                        "DIC",
                        "DIN",
                        "SIL"]

    Obs_surf_diad_var = ["PP",
                     "OCN_PCO2",
                     "CO2FLUX"]

    Obs_3D_compare_var = ["OXY",
                      "ALK",
                      "DIC",
                      "DIN",
                      "SIL"]

    monthnm = {
         "01" : "Jan",
         "02" : "Feb",
         "03" : "Mar",
         "04" : "Apr",
         "05" : "May",
         "06" : "Jun",
         "07" : "Jul",
         "08" : "Aug",
         "09" : "Sep",
         "10" : "Oct",
         "11" : "Nov",
         "12" : "Dec" }

'''


def nc_obs2D_varname(var) :
#dimensions:
#        x = 362 ;
#        y = 332 ;
#variables:
#        double nav_lon(y, x) ;
#        double nav_lat(y, x) ;
#        double obs_din(y, x) ;
#        double obs_sil(y, x) ;
#        double obs_oxy(y, x) ;
#        double obs_dic(y, x) ;
#        double obs_alk(y, x) ;
#        double obs_chl(y, x) ;
#        double obs_npp(y, x) ;
#
    if var == 'ALK' :
        obsvar = 'obs_alk'
    elif var == 'DIC' :
        obsvar = 'obs_dic'
    elif var == 'DIN' :
        obsvar = 'obs_din'
    elif var == 'SIL' :
        obsvar = 'obs_sil'
    elif var == 'PP' :
        obsvar = 'obs_npp'
    elif var == 'CHL' :
        obsvar = 'obs_chl'
    else :
        print("E R R O R -- No obs for this var -- ")
    #
    return obsvar








def mat_obs2D_varname(var) :
    #'rg_glodap1_surf_alk', 'rg_glodap1_surf_dic', 'rg_glodap1_surf_pidic',
    #'rg_glodap2_surf_alk', 'rg_glodap2_surf_dic', 'rg_glodap2_surf_nit',
    #'rg_glodap2_surf_pidic', 'rg_glodap2_surf_sil', 'rg_modis_npp',
    #'rg_roden_clim_flux', 'rg_roden_clim_pco2', 'rg_seawifs_surf_chl',
    #'rg_woa13_surf_nit', 'rg_woa13_surf_oxy', 'rg_woa13_surf_sil'
    if var == 'ALK' :
        obsvar = 'rg_glodap1_surf_alk'
    elif var == 'DIC' :
        obsvar = 'rg_glodap1_surf_dic'
    elif var == 'DIN' :
        obsvar = 'rg_woa13_surf_nit'
    elif var == 'SIL' :
        obsvar = 'rg_woa13_surf_sil'
    elif var == 'PP' :
        obsvar = 'rg_modis_npp'
    elif var == 'CHL' :
        obsvar = 'rg_seawifs_surf_chl'
    elif var == 'OCN_PCO2' :
        obsvar = 'rg_roden_clim_pco2'
    elif var == 'CO2FLUX' :
        obsvar = 'rg_roden_clim_flux'

    return obsvar

def mat_obs3D_varname(var) :
    #'woa75_alk', 'woa75_cfc11', 'woa75_dic', 'woa75_din', 'woa75_oxy', 'woa75_sil'
    if var == 'ALK' :
        obsvar = 'woa75_alk'
    elif var == 'DIC' :
        obsvar = 'woa75_dic'
    elif var == 'DIN' :
        obsvar = 'woa75_din'
    elif var == 'SIL' :
        obsvar = 'woa75_sil'
    elif var == 'CFC11' :
        obsvar = 'woa75_cfc11'
    elif var == 'OXY' :
        obsvar = 'woa75_oxy'

    return obsvar



#################
#################
### Read and plot functions - if need to update - update on jupyter, then update back here






class runclass:

    def __init__(self, modlist):

        #
        print(modlist)
        #
        if modlist == 0 :
            name = "Half sink"
            exp_name = adddict
            dict_ptrc = adddict_ptrc
            dict_diad = adddict_diad
            dict_grid = adddict_grid
            #outfreq   = "1M"
            outfreq   = "5D"
        elif modlist == 1 :
            name = "Dble sink"
            exp_name = rundict
            dict_ptrc = rundict_ptrc
            dict_diad = rundict_diad
            dict_grid = rundict_grid
            outfreq   = "5D"
        elif modlist == 2 :
            name = "MEDUSA"
            exp_name = rundictold
            dict_ptrc = olddict_ptrc
            dict_diad = olddict_diad
            dict_grid = olddict_grid
            outfreq   = "5D"
        else :
            print("modlist out of list")

        self.name = name
        self.exp_name = exp_name
        self.dict_ptrc = dict_ptrc
        self.dict_diad = dict_diad
        self.dict_grid = dict_grid
        self.outfreq   = outfreq


        
def diad_replace(file_name) :
    '''
    if not the diad file, change file name to make it the diad file
    '''
    if "diad" not in file_name :
        ## We need a way to know the diad file (use it for ORCA1 2D cube template)
        try :
            file_name = file_name.replace("ptrc", "diad")
        except :
            file_name = file_name.replace("grid", "diad")
    return file_name



def ptrc_replace(file_name) :
    '''
    if not the ptrc file, change file name to make it the ptrc file
    '''
    if "ptrc" not in file_name :
        ## We need a way to know the diad file (use it for ORCA1 2D cube template)
        try :
            file_name = file_name.replace("diad", "ptrc")
        except :
            file_name = file_name.replace("grid", "ptrc")
    return file_name



def grid_replace(file_name) :
    '''
    if not the grid file, change file name to make it the grid file
    '''
    if "grid" not in file_name :
        ## We need a way to know the diad file (use it for ORCA1 2D cube template)
        try :
            file_name = file_name.replace("ptrc", "grid")
        except :
            file_name = file_name.replace("diad", "grid")
    return file_name



def make_big_average_timeseries (var, runlist, run_exp_list, rundict, out_file="ptrc", MASK=None, region=None, Regrid=False, nemogrid="ORCA1_NEMO42", machine="NOC_linux", year=1850 ):
    '''
    call as : 
    make_big_average_timeseries (var, runlist, run_exp_list, out_file="ptrc", year=1850)
    where 
    -- var is the variable to read
    -- runlist is the list of possible experiments(used in the dict)
    -- run_exp_list is the list of exp to go through 
    -- out_file can be "ptrc" "diad" or "grid" 
    -- year the year to treat
    '''

    from os.path import exists
    #from set_runs_input_for_timeseries import rundict_file
    from set_runs_input_to_plot import rundict_file

    # def file_name for given year and run
    # check in saved registered array if year is already here for this run
    # if not --> read the file and avg and save it.
    #------------------------
    ### prep regions : 
    
    ## regional plot :
    if region == "Amazon" :
        jj_min = 175
        jj_max = 220
        ii_min = 220
        ii_max = 260
    elif region == "Bengal" :
        jj_min = 200
        jj_max = 230
        ii_min = 0
        ii_max = 30
    elif region == "UKshelv" :
        jj_min = 250
        jj_max = 280
        ii_min = 275
        ii_max = 295
        #ii_min=-7
        #ii_max=2
        #jj_min=45
        #jj_max=55
    elif region == "NorthOrtho" or region == "NorthPolarStereo":
        jj_min = 275
        jj_max = 332
        ii_min = 0
        ii_max = 361
        MASK="arc_msk"
    
    ### get the output files name for this year - all experiments
    
    [rundict_ptrc, rundict_diad, rundict_grid] = rundict_file(year, runlist, rundict)
    print(rundict_ptrc)
    for rundd in [rundict_ptrc, rundict_diad, rundict_grid] :
        for ii in run_exp_list :
            if out_file == "ptrc" :
            	bibi = rundict_ptrc[runlist[ii]]
            elif out_file == "diad" :
            	bibi = rundict_diad[runlist[ii]]
            elif out_file == "grid" :
            	bibi = rundict_grid[runlist[ii]]
            
            if type(bibi) is list : 
                try : 
                    if out_file == "ptrc" :
                        rundict_ptrc.update({runlist[ii] : bibi[0]})
                    if out_file == "diad" :
                        rundict_diad.update({runlist[ii] : bibi[0]})
                    if out_file == "grid" :
                        rundict_grid.update({runlist[ii] : bibi[0]})
                except : 
                    print('!! list problem -- probably no file or file corrupted !!')
                    EXIT = True
    ##
    ### run through all runs of the list
    for exp in run_exp_list :
        EXIT=False
        ###
        ##1 - read saved time_serie if exist
        if region is None :
            tm_name=var+"_"+runlist[exp]+"_time_serie"
        else : 
            tm_name=var+"_"+runlist[exp]+"_"+region+"_time_serie"
        print(tm_name)    
        if exists(tm_name+".npy") : 
            tm = np.load(tm_name+".npy")
            print(tm.shape)
            #print(tm)
            NEW = False
            ### Check if year already filled for this run : 
            if year in tm[0,:] : 
                EXIT=True
        else :
            NEW = True
            
        ###
        if EXIT is not True :
            ##2 - read the output 
            if out_file == "ptrc" :
                ref_file_name = rundict_ptrc[runlist[exp]]
            elif out_file == "diad" :
                ref_file_name = rundict_diad[runlist[exp]]
            elif out_file == "grid" :
                ref_file_name = rundict_grid[runlist[exp]]
            #
            ref_ttl       = runlist[exp]
            #
            try : 
                cube = prep_cube(ref_file_name, var,Obs=False,MASK_ZERO=True )
            except : 
                print(' !! SOMETHING WRONG IN ', ref_file_name )
                EXIT = True
            ### we have the field --> now average :
            # need  surface only :
 
        if EXIT is not True :        
            try :
                cube = prep_cube_surf(cube)
                if cube is None :
                    print("cube ref surface gives None")
            except :
                print("already 2D -- ", cube.data.shape ," no need to extract ")
        ##
            ## mask if needed :
            if MASK is not None :
                cube = prep_Masked_cube(cube,MASK,nemogrid=nemogrid, machine=machine)
                if cube is None :
                    print("cube ref Masking gives None")
            ##
            ## regrid to get a real average

            if Regrid :
                print("Regridding cube_ref")
                cube = regrid_ORCA(cube, runlist, Meth="Surface", run_nb=exp, MASK=MASK, NEW=True)
                if cube is None :
                    print("cube ref Regrid gives None")
                print("Regridded cube : ", cube.shape)
            #print("cube_ref : ", cube_ref.shape)
    ##

            ##
            ## extract region if needed     
            #if region is not None :
            #    cubi = cube[...,jj_min:jj_max,ii_min:ii_max].copy()
            #else :
            #    cubi = cube.copy()
            ### check 
            #print("cube shape after extractions :", cubi.shape )
        
            ### now average : 
            ## can be done proprely, but try first the easy way : 
            #### we want to regrid the surface, and then cube collapse lat/lon or nanmean. \o/
            #### no need for region extraction --
            avg = np.nanmean(cube.data)
            #print(year, avg)
            ytm=np.empty([2,1])
            ytm[0,0]=year
            ytm[1,0] = avg
            if NEW :
                ntm=ytm
            else : 
                ntm=np.empty([2,tm.shape[1]+1])
                ntm[:,:-1]=tm
                ntm[:,-1]=ytm[:,0]
            #print("new array's shape :",ntm.shape)
            #print(ntm)
            np.save(tm_name, ntm)
            ## 
            ## should be it
        
        


        
def fix_nemo_cube(src_cubi, mesh_mask="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/MEDUSA_well_defined.nc") :
    ## cube badly written - have to fix it...
    ## there is no x and y and z coord in the 4.2 NEMO files
    ## get a file with the right coord :
    tmp_src= prep_cube(mesh_mask,"PRN")
    #print(tmp_src)
    # and create new, well defined, cube
                    ###
    print(" NO coord... fix NEMO file")
    # guess 2D or 3D
    regrid2d = False
    regrid3d = False
    if src_cubi.data.ndim == 4 :
        print("-- Fixing 3D --")
        regrid3d = True
    elif src_cubi.data.ndim == 3 :
        ## check the first dim is depth (more than one layer)
        ## if not, no need to extract
        if src_cubi.data.shape[0] <= 12 :
            print("-- Fixing 2D --")
            regrid2d = True
        else :
            # real 3D
            print("-- Fixing 3D --")
            regrid3d = True
    elif src_cubi.data.ndim == 2 :
        # already 2D
        print("-- Fixing 2D --")
        regrid2d = True
        #area = src_cubi.copy(e1t.data[0,:,:] * e2t.data[0,:,:])
    else :
        print("data shape not understood : ", src_cubi.data.shape)

    # create Cube :
    #print("area shape     : ", area.shape)
    xx = tmp_src.coord(axis='x')
    #print(xx)
    yy = tmp_src.coord(axis='y')
    #print(yy)
    #tt = src_cubi.coords(axis='t')[0]
    #print(tt)
    ##
    #ttdim = src_cubi.shape[0]
    #yydim = tmp_src.shape[1]
    #xxdim = tmp_src.shape[2]
    ##
    if regrid3d :
        src_cubi.add_aux_coord(yy, data_dims=[2, 3])
        src_cubi.add_aux_coord(xx, data_dims=[2, 3])
        #zz = src_cubi.coord(axis='z')
        #print(zz)
        #zzdim = src_cubi.shape[1]
        #reg_cubi = iris.cube.Cube(np.zeros((ttdim, zzdim, yydim, xxdim)),
        #                dim_coords_and_dims=[(tt, 0),
        #                                    (zz, 1),
        #                                    (yy, 2),
        #                                    (xx, 3)])
    else:
        #reg_cubi = iris.cube.Cube(np.zeros((ttdim, yydim, xxdim)))
        #reg_cubi.add_dim_coord(tt, 0)
        src_cubi.add_aux_coord(yy, data_dims=[1, 2])
        src_cubi.add_aux_coord(xx, data_dims=[1, 2])

    #reg_cubi.data = src_cubi.data
    #reg_cubi.name = src_cubi.name
    #reg_cubi.var_name = src_cubi.var_name
    #reg_cubi.long_name = src_cubi.long_name
    #reg_cubi.units = src_cubi.units

    return src_cubi





def regrid_ORCA(src_cubi, runlist, Meth=None, NEW = False, run_nb=None, MASK=None, reg_grid_file = "/noc/users/jpp1m13/WORKING/UKESM/MESH/Data_reg_grid/woa13_all_i00_01.nc") :
    '''
    regrid_ORCA(src_cubi, Meth=None, NEW = False, run_nb=run_nb) :
    regrids on WOA 1x1 grid (if not already here ) and saves it 
    -- src_cubi : source cube to regrid
    -- runlist : same as in subplot_diff -- list of available run names, used to name the regrid file 
    -- if Meth is Surface : regrid only surface field
    -- NEW to reproduce the file -- to regrid again
    -- run_nb : run number to get the run name for the file to save.
    -- MASK   : if you want to mask the file to regrid (mask is the same than in subplot_diff )
    -- reg_grid_file : a regrid file used a ref model for the grid -- default WOA 1x1 grid.
    '''
    import os

    dirr="SAVE_NC_FILE"  ## Directory with all saved regrided nc files
    ##
    ##first check dir exists and create it if needed:
    if os.path.isdir(dirr) != True :
        os.mkdir(dirr)
    ##
    ## regrid :
    ## Target grid : WOA 
    reg_cube = prep_cube(reg_grid_file,"i_mn")
    #if Meth=="Surface":
    reg_cube = prep_cube_surf(reg_cube)
    ##
    if Meth=="Surface":
        print("Regrid -- extract surf -- cube shape : ",src_cubi.shape)
        src_cubi = prep_cube_surf(src_cubi)
        if src_cubi is None : 
            print("Regrid -- cube Surface gives None")
    ##
    ## extra precautions -- mask problematic values
    ## target -- should not be needed
    reg_cube.data = np.ma.masked_invalid(reg_cube.data)
    reg_cube.data = np.ma.masked_where(reg_cube.data==0.0, reg_cube.data)
    reg_cube_mask = reg_cube.data.mask
    ## sources 
    src_cubi.data = np.ma.masked_invalid(src_cubi.data)
    src_cubi.data = np.ma.masked_where(src_cubi.data==0.0, src_cubi.data)

        #print(reg_cube)
    reg_cube.coord(axis='x').coord_system = iris.coord_systems.GeogCS(6371229.0)
    reg_cube.coord(axis='y').coord_system = iris.coord_systems.GeogCS(6371229.0)
    #reg_grid = reg_cube.coords()
    #print(reg_grid)
    ## Source grid to regrid: NEMO
    #src_cubi = prep_cube("Out/medusa_cp799o_sum-1m_20200101-21001231_diad-T.nc","TOT_SHALK")
    regrid2d = False
    regrid3d = False
    print("src cubi shape : ", src_cubi.shape)
    print("reg cube shape : ", reg_cube.shape)
    src_mask = src_cubi.copy()
    src_mask.data[...] = 1
    src_mask.data[src_cubi.data.mask] = 0
    
    if src_cubi.data.ndim == 4 :
        print("-- 3D Regrid --")
        regrid3d = True
    elif src_cubi.data.ndim == 3 :
        ## check the first dim is depth (more than one layer)
        ## if not, no need to extract
        if src_cubi.data.shape[0] <= 12 :
            print("-- 2D Regrid --")
            regrid2d = True
        else :
            # real 3D
            print("-- 3D Regrid --")
            regrid3d = True
    elif src_cubi.data.ndim == 2 :
        # already 2D
        print("-- 2D Regrid --")
        regrid2d = True
        #area = src_cubi.copy(e1t.data[0,:,:] * e2t.data[0,:,:])
    else : 
        print("data shape not understood : ", src_cubi.data.shape)
   ###
   ### need src_cubi coord system (but might need to fix the cubes...)
    try:
        src_cubi.coord(axis='x').coord_system = reg_cube.coord(axis='x').coord_system
        src_cubi.coord(axis='y').coord_system = reg_cube.coord(axis='y').coord_system
    except :
        src_cubi = fix_nemo_cube(src_cubi)
        src_cubi.coord(axis='x').coord_system = reg_cube.coord(axis='x').coord_system
        src_cubi.coord(axis='y').coord_system = reg_cube.coord(axis='y').coord_system
    ##
    ## Scheme --
    #scheme = iris.analysis.AreaWeighted(mdtol=0.5)
    #scheme = iris.analysis.PointInCell()
    scheme = iris.analysis.UnstructuredNearest()
    #scheme = iris.analysis.Nearest()
    #scheme = iris.analysis.Linear()
    ##
    isfile = False
    if MASK==None : 
        MASK="no"
    ##
    if regrid2d : 
        fname = dirr + "/regridded_" + src_cubi.var_name + "_" + runlist[run_nb] + "_" + MASK + "_mask_2D.nc"
        print("check if 2D file exists : ", fname)
        if os.path.isfile(fname) and not NEW: 
             # file exists : read it 
            print("read regridded file : ",fname)
            reg_cubi = iris.load_cube(fname)
            isfile = True
        else :     
            # create it 
            print("No file to read : regridding")
            reg_cubi=src_cubi.regrid(reg_cube, scheme)
            #reg_mask=src_mask.regrid(reg_cube, scheme)
    elif regrid3d : 
        fname = dirr + "/regridded_" + src_cubi.var_name + "_" + runlist[run_nb] + "_" + MASK + "_mask_3D.nc"
        print("check if 3D file exists : ", fname)
        if os.path.isfile(fname) and not NEW: 
            # file exists : read it 
            print("read regridded file : ",fname)
            reg_cubi = iris.load_cube(fname)
            isfile = True
        else :     
            ### 
            print("No file to read : regridding")
            # create Cube :
            #print("area shape     : ", area.shape)
            xx = reg_cube.coord(axis='x')
            yy = reg_cube.coord(axis='y')
            #print(reg_cube)
            tt = src_cubi.coords(axis='t')[0]
            zz = src_cubi.coord(axis='z')
            #if zz==None :
            #    zz = src_cubi.coord('Vertical T levels')
            #print(cc)
            ttdim = src_cubi.shape[0]
            zzdim = src_cubi.shape[1]
            yydim = reg_cube.shape[1]
            xxdim = reg_cube.shape[2]

            reg_cubi = iris.cube.Cube(np.zeros((ttdim, zzdim, yydim, xxdim)),
                        dim_coords_and_dims=[(tt, 0),
                                             (zz, 1),
                                             (yy, 2), 
                                             (xx, 3)])
            reg_cubi.name = src_cubi.name
            reg_cubi.var_name = src_cubi.var_name
            reg_cubi.long_name = src_cubi.long_name
            reg_cubi.units = src_cubi.units
            #reg_mask = reg_cubi.copy()
            for ZZZ in np.arange(zzdim) : 
                print("Regrid Z level --",ZZZ)
                print("src cubi slice shape : ", src_cubi[:,ZZZ,:,:].shape)
                tcube=src_cubi[:,ZZZ,:,:].regrid(reg_cube, scheme)
                reg_cubi.data[:,ZZZ,:,:] = tcube.data
                #tmask = src_mask[:,ZZZ,:,:].regrid(reg_cube, scheme)
                #reg_mask.data[:,ZZZ,:,:] = tmask.data
            
    
    if not isfile : 
        ## need a masked array
        #bibi = np.ma.masked_where(reg_mask.data == 0, reg_cubi.data)
        bibi = np.ma.masked_where(reg_cubi.data > 1e17, reg_cubi.data)
        reg_cubi.data = bibi
        ## Save file : 
        print("saving file : ",fname)          
        iris.save(reg_cubi, fname)
    else : 
        print(reg_cubi)
        
    return reg_cubi






def load_meta(datapath, fxpath=None):
    """Load metadata of the netCDF file.

    Parameters
    ----------
    datapath: str
        path to the netCDF file with data
    fxpath: str
        path to the netCDF file with fx files

    Returns
    -------
    datafile: instance of netCDF4 Dataset
        points to the file
    lon2d: numpy array
        Two dimentional longitude information
    lat2d: numpy array
        Two dimentional latitude information
    lev: numpy array
        depths of the model levels
    time: numpy array
        dates are converted to datetime objects
    areacello: numpy array
        values of areacello
    """
    from netCDF4 import Dataset
    
    datafile = Dataset(datapath)

    if fxpath:
        datafile_area = Dataset(fxpath)
        areacello = datafile_area.variables['areacello'][:]
    else:
        areacello = None

    lon = datafile.variables['lon'][:]
    lat = datafile.variables['lat'][:]
    lev = datafile.variables['lev'][:]
    time = num2date(datafile.variables['time'][:],
                    datafile.variables['time'].units)
    # hack for HadGEM2-ES
    lat[lat > 90] = 90

    if lon.ndim == 2:
        lon2d, lat2d = lon, lat
    elif lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)

    metadata = {}
    metadata['datafile'] = datafile
    metadata['lon2d'] = lon2d
    metadata['lat2d'] = lat2d
    metadata['lev'] = lev
    metadata['time'] = time
    metadata['areacello'] = areacello
    return metadata






def plot_budget_bar(var, xx, globi, exp_names, yr = None, Convert=1.0, ConvUnit = "None", op = "sum", log_scale = False, centered=False, Titre = None, REF_COMP=None, RETURN=False, **attrs) :

    print(var)

    if op == "sum" :
        f_name = "Glub_sum_" + var
    elif op == "avg" :
        f_name = "Glov_avg_" + var


    comp_nb = np.size(globi)

    fig = plt.figure(figsize=(7.5,7.5))
    ax = fig.add_subplot(1, 1, 1)
    #fig = plt.figure(figsize=(5*6,5*6))
    ## no pac


    ##
    ## plot :
    print(xx)
    print(exp_names)
    print(globi)
    p1 = ax.bar(xx, globi*Convert, edgecolor="white", linewidth=0.8)
    if ConvUnit == "None" :
        ax.set_ylabel("-")
    else :
        ax.set_ylabel(ConvUnit)
    ax.bar_label(p1, label_type='center', size=12)
    ax.set_xticks(xx, labels= exp_names, size=12 )
    #ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
    #   ylim=(0, 8), yticks=np.arange(1, 8))

    if Titre == None :
        plt.suptitle( "No title.... add one" , size=16)
    else :
        plt.suptitle( Titre , size=16)
    plt.rcParams.update({'font.size': 12})
    plt.savefig(f_name)
    plt.show()

    if RETURN :
        return xx, globi, exp_names



def comp_budget_bar(var, rundict, run_exp_list, Meth="None", yr="None", Convert=1.0, ConvUnit = "None", op = "sum", log_scale = False, centered=False, Diff = False, Titre = None, REF_COMP=None, RETURN=False, **attrs) :

    print(var)

    if op == "sum" :
        f_name = "Glob_sum_" + var
    elif op == "avg" :
        f_name = "Glob_avg_" + var
    if yr != "None" :
        f_name = f_name + "_" + yr

    if REF_COMP is not None :
        runcomp = REF_COMP
        ref_comp_fname = rundict[runcomp]
        comp_ttl       = runcomp
    else :
        ref_comp_fname = None

    comp_nb = np.size(run_exp_list)

    fig = plt.figure(figsize=(7.5,7.5))
    ax = fig.add_subplot(1, 1, 1)
    #fig = plt.figure(figsize=(5*6,5*6))
    ## no pac


    loop = 0
    xx = np.arange(comp_nb)
    for compref in run_exp_list :
        comp_file_name = rundict[runlist[compref]]
        comp_ttl       = runlist[compref]

        cube = prep_cube(comp_file_name, var, REF_COMP=ref_comp_fname)
        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Surface":
            try :
                cube = prep_cube_surf(cube)
            except :
                print("already 2D -- no need to extract ")
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ")
        ##
        if Diff == False :
            if loop == 0 :
                globi = np.array([glob_budget_2D(cube, op = op)])
                exp_names = np.array([comp_ttl])
            else :
                globi = np.append(globi,[glob_budget_2D(cube, op = op)])
                exp_names = np.append(exp_names,[comp_ttl])
        else :
            if loop == 0 :
                globref = np.array([glob_budget_2D(cube, op = op)])
            elif loop == 1 :
                globi = np.array([glob_budget_2D(cube, op = op) - globref])
                exp_names = np.array([comp_ttl])
            else :
                globi = np.append(globi,[glob_budget_2D(cube, op = op) - globref])
                exp_names = np.append(exp_names,[comp_ttl])
        loop = loop + 1
    ##
    ## plot :
    if Diff == True :
        xx = xx[1:]
    print(xx)
    print(exp_names)
    print(globi)
    p1 = ax.bar(xx, globi*Convert, edgecolor="white", linewidth=0.8)
    if ConvUnit == "None" :
        ax.set_ylabel(cube.units)
    else :
        ax.set_ylabel(ConvUnit)
    ax.bar_label(p1, label_type='center', size=12)
    ax.set_xticks(xx, labels= exp_names, size=12 )
    #ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
    #   ylim=(0, 8), yticks=np.arange(1, 8))
    if Diff == False :
        if op == "sum" :
            tittle = "Global sum of " + cube.long_name
        elif op == "avg" :
            tittle = "Global avg of " + cube.long_name
    else :
        if op == "sum" :
            tittle = "Global sum of " + cube.long_name + " -- diff with " + runlist[0]
        elif op == "avg" :
            tittle = "Global avg of " + cube.long_name + " -- diff with " + runlist[0]
    if Titre == None :
        plt.suptitle( tittle , size=16)
    else :
        plt.suptitle( Titre , size=16)
    plt.rcParams.update({'font.size': 12})
    plt.savefig(f_name)
    plt.show()

    if RETURN :
        return xx, globi, exp_names





def glob_budget_2D(cube, op = "sum",meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc") :

    def glob_sum(cube, area, op = "sum") :
        #try :
        #cube = read_cube(file_name, var)
        area.mask = cube.data.mask
        ppp = cube.data.copy()
        ppp1 = cube.data.copy()
        if op == "sum" :
            ppp1 = ppp * area ## mmol--2D
            #
            ppp2 = ppp1.sum(axis = 1)
            ppp3 = ppp2.sum(axis = 1)
            gls = ppp3[0]  # glob sum Tmol-N
        elif op == "avg" :
            print("averaging")
            ppp2 = ma.masked_where(ppp1 == np.nan, ppp1 )
            marea = ma.masked_where(ppp1 == np.nan, area )
            ppp3 = ppp2.compressed()
            carea = marea.compressed()
            gls = np.average(ppp3, weights=carea)
            #gls = ppp4[0]
        #except :
        #    ## var not in the file (no pac probably)
        #    gls = 0.0
            print("gls = ",gls)
        #    #
        return gls

    ## all vars are 2D inventory, only needs 2D area and sum-up
    #meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/mesh_mask.nc"
    e1t = prep_cube(meshfile, "e1t")
    e2t = prep_cube(meshfile, "e2t")
    area = e1t.data.copy()
    area = e1t.data * e2t.data

    global_sum     = glob_sum(cube, area, op)
    return global_sum






def prod_graz_budget(file_name, run_ttl,meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc") :

    def glob_sum(file_name, var, area) :
        try :
            cube = read_cube(file_name, var)
            area.mask = cube.data.mask
            ppp = cube.data.copy()
            ppp1 = cube.data.copy()
            ppp1 = ppp * area ## mmol-N/d
            #
            ppp2 = ppp1.sum(axis = 1)
            ppp3 = ppp2.sum(axis = 1)
            gls = ppp3[0] * 360 /1E15 # glob sum Tmol-N/y
        except :
            ## var not in the file (no pac probably)
            gls = 0.0
            #
        return gls

    ## all vars are 2D inventory, only needs 2D area and sum-up
    try:
        ## running local
        #meshfile="../MESH/eORCA1/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
    except:
        ## on NOC system 
        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
    e1t = prep_cube(meshfile, "e1t")
    e2t = prep_cube(meshfile, "e2t")
        ##
    area = e1t.data.copy()
    area = e1t.data * e2t.data

    glob_PRN     = glob_sum(file_name, "PRN", area)
    glob_PRD     = glob_sum(file_name, "PRD", area)
    glob_ZIGRO   = glob_sum(file_name, "ZI_GROW", area)
    glob_GMIPN   = glob_sum(file_name, "GMIPN", area)
    glob_GMID    = glob_sum(file_name, "GMID", area)
    glob_ZEGRO   = glob_sum(file_name, "ZE_GROW", area)
    glob_GMEPN   = glob_sum(file_name, "GMEPN", area)
    glob_GMEPD   = glob_sum(file_name, "GMEPD", area)
    glob_GMEZMI  = glob_sum(file_name, "GMEZMI", area)
    glob_GMEZMP  = glob_sum(file_name, "GMEZMP", area)
    glob_GMED    = glob_sum(file_name, "GMED", area)
    glob_GMEDF1  = glob_sum(file_name, "GMEDF1", area)
    glob_GMEDF2  = glob_sum(file_name, "GMEDF2", area)
    glob_ZPGRO   = glob_sum(file_name, "ZP_GROW", area)
    glob_GMPPN   = glob_sum(file_name, "GMPPN", area)
    glob_GMPPD   = glob_sum(file_name, "GMPPD", area)
    glob_GMPZMI  = glob_sum(file_name, "GMPZMI", area)
    glob_GMPD    = glob_sum(file_name, "GMPD", area)
    glob_GMPDF1  = glob_sum(file_name, "GMPDF1", area)
    glob_GMPDF2  = glob_sum(file_name, "GMPDF2", area)

    tgrazZMI = glob_GMIPN + glob_GMID
    tgrazZME = glob_GMEPN + glob_GMEPD + glob_GMEZMI + glob_GMEZMP + glob_GMED + glob_GMEDF1 + glob_GMEDF2
    tgrazZMP = glob_GMPPN + glob_GMPPD + glob_GMPZMI + glob_GMPD + glob_GMPDF1 + glob_GMPDF2

    f= open(run_ttl + "_graz_prod_budget.csv","w")

    f.writelines(" , , " +  run_ttl  + " in T mol-N/Y \n")
    f.writelines(" , PRN, PRD, ZMI growth, ZME gowth, PAC growth \n")
    f.writelines(" , " + str(glob_PRN) + "," + str(glob_PRD) + "," + str(glob_ZIGRO) + "," + str(glob_ZEGRO) + "," + str(glob_ZPGRO) + " \n")
    f.writelines(" , loss to ZMI, loss to ZME, loss to PAC  \n")
    f.writelines("PN ," + str(glob_GMIPN) + "," + str(glob_GMEPN) + "," + str(glob_GMPPN) + " \n")
    f.writelines("PD ,0.0," + str(glob_GMEPD) + "," + str(glob_GMPPD) + " \n")
    f.writelines("slow det ," + str(glob_GMID) + "," + str(glob_GMED) + "," + str(glob_GMPD) + " \n")
    f.writelines("Fast det 1 , 0.0 ," + str(glob_GMEDF1) + "," + str(glob_GMPDF1) + " \n")
    f.writelines("Fast det 2 , 0.0 ," + str(glob_GMEDF2) + "," + str(glob_GMPDF2) + " \n")
    f.writelines("ZMI , 0.0 ," + str(glob_GMEZMI) + "," + str(glob_GMPZMI) + " \n")
    f.writelines("PAC , 0.0 ," + str(glob_GMEZMP) + ", 0.0 \n")
    f.writelines("Total grazing , "+ str(tgrazZMI) + "," + str(tgrazZME) + "," + str(tgrazZMP) + " \n")
    f.close()






class subplot_POC_flx_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, out_path , var = "CONC",
                 fig=None, iii=1, jjj=1, rect=1, **attrs):

        """
    Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']


        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            speed_pref_diad = out_path + "/" + diadf(yr)
            speed_pref_dext = out_path + "/" + diadf(yr)
            speed_pref_ptrc = out_path + "/" + ptrcf(yr)

            try:
                FD2 = prep_cube(speed_pref_dext, "FD2_CAR3")
                #grazpn.standard_name = "grazpn"
                print(FD2)
                FD1 = prep_cube(speed_pref_dext, "FD1_CAR3")
                #grazdn.standard_name = "grazdn"
                print(FD1)
                FDS = prep_cube(speed_pref_ptrc, "DTC")
                FDS.data = FDS.data * 0.5
                #grazzmi.standard_name = "grazzmi"
                print(FDS)

                if var == "CONC" :
                    FDS.data = FDS.data / 0.5
                    FD1.data = FD1.data / 35.0
                    FD2.data = FD2.data / 115.0
                    FDS.long_name="Detritus concentration"
                    FDS.units = "mmol-C/m3"
                    tot =  FDS.copy()
                    tot.data = FDS.data + FD1.data + FD2.data
                else :
                    tot =  FDS.copy()
                    tot.data = FDS.data  + FD1.data  + FD2.data
                    tot.long_name="Detritus flux"

                cubelist = {FDS:"slow det", FD1:"Fast det", FD2:"Fast det 2", tot:"tot det"}

            except :
                FD1 = prep_cube(speed_pref_dext, "FD1_CAR3")
                #grazdn.standard_name = "grazdn"
                print(FD1)
                FDS = prep_cube(speed_pref_ptrc, "DTC")
                FDS.data = FDS.data * 0.5
                #grazzmi.standard_name = "grazzmi"
                print(FDS)

                if var == "CONC" :
                    tot =  FDS.copy()
                    FDS.long_name="Detritus concentration"
                    tot.data = FDS.data / 0.5 + FD1.data / 35.0
                else :
                    tot =  FDS.copy()
                    tot.long_name="Detritus flux"
                    tot.data = FDS.data  + FD1.data

                cubelist = {FDS:"slow det", FD1:"Fast det", tot:"tot det"}


            try :
                depth = FDS.coord('Vertical T levels').points
            except :
                depth = FDS.coord('depth').points


            #print(depth)
            for cube in cubelist :
                if loc == "PAP" :
                    jj = 258
                    ii = 272
                    pocprof = cube.data[:,:,jj,ii].copy()
                elif loc == "BATS" :
                    jj = 258
                    ii = 272
                    pocprof = cube.data[:,:,jj,ii].copy()
                elif loc == "GLOB" :
                    area = cube.data.copy()
                    for zz in np.arange(75) :
                        area[:,zz,:,:] = area2D[:,:,:]

                    area.mask = cube.data.mask
                    ppp = cube.data.copy()
                    ppp1 = cube.data.copy()
                    ppp1 = ppp * area

                    ppp2 = ppp1.sum(axis = 2)
                    ppp3 = ppp2.sum(axis = 2)

                    area1 = area.sum(axis = 2)
                    area2 = area1.sum(axis = 2)

                    pocprof = ppp3 / area2

                ax.plot(pocprof[0,:],depth[:],label = cubelist[cube] + " - mm " + yr )


            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)

        ###
        ## manage plot
        ax.legend()
        ax.set_ylim([0.0, 5000])
        ax.invert_yaxis()
        if var == "CONC" :
            ax.set_xlim([0.0, 3])
        else :
            ax.set_xlim([0.0, 7])
        ax.set_xlabel(cube.units)
        ax.set_ylabel("Depth (m)")
        if 'fig_ttl' in attrs:
            fig_ttl = attrs['fig_ttl']
            ax.set_title(cube.long_name + " - " + fig_ttl)
        else :
            ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)





class subplot_PAC_grazing_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, out_path , var = "TOT",
                 fig=None, iii=1, jjj=1, rect=1, LOGSCALE=True,**attrs):

        """
    Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']


        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            speed_pref_diad = out_path + "/" + diadf(yr)
            speed_pref_dext = out_path + "/" + diadf(yr)
            speed_pref_ptrc = out_path + "/" + ptrcf(yr)

            grazpn = prep_cube(speed_pref_dext, "GMPPN3")
            #grazpn.standard_name = "grazpn"
            print(grazpn)
            grazdn = prep_cube(speed_pref_dext, "GMPPD3")
            #grazdn.standard_name = "grazdn"
            print(grazdn)
            grazzmi = prep_cube(speed_pref_dext, "GMPZMI3")
            #grazzmi.standard_name = "grazzmi"
            print(grazzmi)
            grazd = prep_cube(speed_pref_dext, "GMPD3")
            #grazd.standard_name = "grazd"
            print(grazd)
            grazdf1 = prep_cube(speed_pref_dext, "GMPDF13")
            #grazdf1.standard_name = "grazdf1"
            print(grazdf1)
            grazdf2 = prep_cube(speed_pref_dext, "GMPDF23")
            #grazdf2.standard_name = "grazdf2"
            print(grazdf2)
            ZMP = prep_cube(speed_pref_ptrc, "ZMP")
            #ZMP.standard_name = "ZMP"
            print(ZMP)
            ## adapt mask
            grazpn.data.mask = ZMP.data.mask
            grazdn.data.mask = ZMP.data.mask
            grazzmi.data.mask = ZMP.data.mask
            grazd.data.mask = ZMP.data.mask
            grazdf1.data.mask = ZMP.data.mask
            grazdf2.data.mask = ZMP.data.mask

            ## Tot graz
            totgraz = grazpn.copy()
            totgraz.data = grazpn.data + grazdn.data + grazzmi.data + grazd.data + grazdf1.data + grazdf2.data
            totgraz.long_name = "PAC total grazing"
            #totgraz.standard_name = "totgraz"
            #Specif graz
            spegraz = ZMP.copy()
            print(totgraz)
            print(ZMP)
            spegraz.data =  totgraz.data / ZMP.data
            spegraz.long_name = "PAC specific grazing"
            #spegraz.standard_name = "spegraz"

            try :
                depth = ZMP.coord('Vertical T levels').points
            except :
                depth = ZMP.coord('depth').points

            if var == "TOT" :
                cubelist = {grazpn:"grazpn", grazdn:"grazdn", grazzmi:"grazzmi", grazd:"grazd", grazdf1:"grazdf1", grazdf2:"grazdf2", totgraz:"totgraz"}
            else :
                cubelist = {spegraz:"spegraz"}

            #print(depth)
            for cube in cubelist :
                if loc == "PAP" :
                    jj = 258
                    ii = 272
                    pocprof = cube.data[:,:,jj,ii].copy()
                elif loc == "BATS" :
                    jj = 258
                    ii = 272
                    pocprof = cube.data[:,:,jj,ii].copy()
                elif loc == "GLOB" :
                    area = cube.data.copy()
                    for zz in np.arange(75) :
                        area[:,zz,:,:] = area2D[:,:,:]

                    area.mask = cube.data.mask
                    ppp = cube.data.copy()
                    ppp1 = cube.data.copy()
                    ppp1 = ppp * area

                    ppp2 = ppp1.sum(axis = 2)
                    ppp3 = ppp2.sum(axis = 2)

                    area1 = area.sum(axis = 2)
                    area2 = area1.sum(axis = 2)

                    pocprof = ppp3 / area2

                ax.plot(pocprof[0,:],depth[:],label = cubelist[cube] + " - mm " + yr )


            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)

        ###
        ## manage plot
        ax.legend()
        if LOGSCALE :
                    ax.set_yscale('log')
        ax.set_ylim([0.0, 5000])
        ax.invert_yaxis()
        if var == "TOT" :
            ax.set_xlim([0.0, 0.05])
        else :
            ax.set_xlim([0.0, 0.4])
        ax.set_xlabel(cube.units)
        ax.set_ylabel("Depth (m)")
        if 'fig_ttl' in attrs:
            fig_ttl = attrs['fig_ttl']
            ax.set_title(cube.long_name + " - " + fig_ttl)
        else :
            ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)







class subplot_3Dvar_profiles_diff(object):
    """
    Plot projected maps
    """

    def __init__(self, var, ref_run_nm, run_nm,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, LOGSCALE=True, Manage=True,
                 YLIM=[0.0, 5000], **attrs):

        """
        Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))
        if ax is None :
            ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']
        else :
            loc = "GLOB"

        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            #speed_pref_diad = out_path + "/eORCA1_MED_UKESM_1m_18540101_18541230_diad_T_1854"+yr+"-1854"+yr+".nc"
            #speed_pref_file = out_path + "/eORCA1_MED_UKESM_1m_0101_18591230_"+ fkind +"_T_1859"+yr+"-1859"+yr+".nc"
            try :
                refcube = prep_cube(rundict_ptrc[ref_run_nm], var)
            except :
                refcube = prep_cube(rundict_diad[ref_run_nm], var)
            #
            try :
                cube = prep_cube(rundict_ptrc[run_nm], var)
            except :
                cube = prep_cube(rundict_diad[run_nm], var)
            #
            #SDDT = prep_cube(speed_pref_ptrc, "DTC")
            #SDDT.data = SDDT.data * 0.5 ## 0.5 m/d
            #SDDT.units = "mmol-C/m2/d"
            #cube =  FDDT.copy()
            #cube.data = FDDT.data + SDDT.data
            #cube.long_name = "sinking detrical carbon flux"
            try :
                depth = cube.coord('Vertical T levels').points
            except :
                depth = cube.coord('depth').points

            if loc == "PAP" :
                jj = 258
                ii = 272
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "BATS" :
                jj = 237
                ii = 224
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "ATL1" :
                jj = 85
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "ATL2" :
                jj = 110
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "ATL3" :
                jj = 150
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "ATL4" :
                jj = 200
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "ATL5" :
                jj = 260
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
                refprof = refcube.data[:,:,jj,ii].copy()
            elif loc == "GLOB" :
                area = cube.data.copy()
                for zz in np.arange(75) :
                    area[:,zz,:,:] = area2D[:,:,:]

                area.mask = cube.data.mask
                ppp = cube.data.copy()
                ppp1 = cube.data.copy()
                ppp1 = ppp * area

                ppp2 = ppp1.sum(axis = 2)
                ppp3 = ppp2.sum(axis = 2)

                rppp = refcube.data.copy()
                rppp1 = refcube.data.copy()
                rppp1 = rppp * area

                rppp2 = rppp1.sum(axis = 2)
                rppp3 = rppp2.sum(axis = 2)

                area1 = area.sum(axis = 2)
                area2 = area1.sum(axis = 2)

                pocprof = ppp3 / area2
                refprof = rppp3 / area2

            diffi =  pocprof - refprof
            print(pocprof.size)
            if 'fig_lab' in attrs:
                lab = attrs['fig_lab']
            else :
                lab = None

            if (yr == "NO") and (Manage==True):
                ax.plot(diffi[0,:],depth[:],'k--', label = lab )
            else :
                ax.plot(diffi[0,:],depth[:],label = lab )


            aaa_l = diffi.min()
            bbb_l = diffi.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)
        ###
        ## manage plot
        ax.legend()
        if Manage==True :
            if LOGSCALE :
                ax.set_yscale('log')
            ax.set_ylim(YLIM)
            ax.invert_yaxis()
            if 'xlim' in attrs:
                xlim = attrs['xlim']
                ax.set_xlim(xlim)
            ax.set_xlabel(cube.units)
            ax.set_ylabel("Depth (m)")
            if 'fig_ttl' in attrs:
                fig_ttl = attrs['fig_ttl']
                ax.set_title(cube.long_name + " - " + fig_ttl)
            else :
                ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)









class subplot_3Dvar_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, var, run_nm,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, LOGSCALE=True, Manage=True,
                 YLIM=[0.0, 5000], **attrs):

        """
       Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))
        if ax is None :
            ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']
        else :
            loc = "GLOB"

        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            #speed_pref_diad = out_path + "/eORCA1_MED_UKESM_1m_18540101_18541230_diad_T_1854"+yr+"-1854"+yr+".nc"
            #speed_pref_file = out_path + "/eORCA1_MED_UKESM_1m_0101_18591230_"+ fkind +"_T_1859"+yr+"-1859"+yr+".nc"
            try :
                cube = prep_cube(rundict_ptrc[run_nm], var)
            except :
                cube = prep_cube(rundict_diad[run_nm], var)
            #
            #SDDT = prep_cube(speed_pref_ptrc, "DTC")
            #SDDT.data = SDDT.data * 0.5 ## 0.5 m/d
            #SDDT.units = "mmol-C/m2/d"
            #cube =  FDDT.copy()
            #cube.data = FDDT.data + SDDT.data
            #cube.long_name = "sinking detrical carbon flux"
            try :
                depth = cube.coord('Vertical T levels').points
            except :
                depth = cube.coord('depth').points

            if loc == "PAP" :
                jj = 258
                ii = 272
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "BATS" :
                jj = 237
                ii = 224
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL1" :
                jj = 85
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL2" :
                jj = 110
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL3" :
                jj = 150
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL4" :
                jj = 200
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL5" :
                jj = 260
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "GLOB" :
                area = cube.data.copy()
                for zz in np.arange(75) :
                    area[:,zz,:,:] = area2D[:,:,:]

                area.mask = cube.data.mask
                ppp = cube.data.copy()
                ppp1 = cube.data.copy()
                ppp1 = ppp * area

                ppp2 = ppp1.sum(axis = 2)
                ppp3 = ppp2.sum(axis = 2)

                area1 = area.sum(axis = 2)
                area2 = area1.sum(axis = 2)

                pocprof = ppp3 / area2
            print(pocprof.size)
            if 'fig_lab' in attrs:
                lab = attrs['fig_lab']
            else :
                lab = None

            if (yr == "NO") and (Manage==True):
                ax.plot(pocprof[0,:],depth[:],'k--', label = lab )
            else :
                ax.plot(pocprof[0,:],depth[:],label = lab )


            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)
        ###
        ## manage plot
        ax.legend()
        if Manage==True :
            if LOGSCALE :
                ax.set_yscale('log')
            ax.set_ylim(YLIM)
            ax.invert_yaxis()
            if 'xlim' in attrs:
                xlim = attrs['xlim']
                ax.set_xlim(xlim)
            ax.set_xlabel(cube.units)
            ax.set_ylabel("Depth (m)")
            if 'fig_ttl' in attrs:
                fig_ttl = attrs['fig_ttl']
                ax.set_title(cube.long_name + " - " + fig_ttl)
            else :
                ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)










class subplot_NormPOCup_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, out_path ,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, LOGSCALE=True,
                 Martin=True, Manage=True, **attrs):

        """
        Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))
        if ax is None:
            ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']

        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            speed_pref_diad = out_path + "/" + diadf(yr)
            speed_pref_dext = out_path + "/" + diadf(yr)
            speed_pref_ptrc = out_path + "/" + ptrcf(yr)

            try:
                FD2 = prep_cube(speed_pref_dext, "FD2_CAR3")
                #grazpn.standard_name = "grazpn"
                print(FD2)
                FD1 = prep_cube(speed_pref_dext, "FD1_CAR3")
                #grazdn.standard_name = "grazdn"
                print(FD1)
                FDS = prep_cube(speed_pref_ptrc, "DTC")
                FDS.data = FDS.data * 0.5
                #grazzmi.standard_name = "grazzmi"
                print(FDS)

                cube =  FDS.copy()
                cube.data = FDS.data + FD1.data  + FD2.data
                #Sexp = FDS.data[0,23,:,:] + FD1.data[0,23,:,:]  + FD2.data[0,23,:,:]

                #for zz in np.arange(75) :
                #    cube.data[0,zz,:,:] = cube.data[0,zz,:,:] / cube.data[0,23,:,:]

                #cube.long_name="Detritus flux"

                #cubelist = {FDS:"slow det", FD1:"Fast det", FD2:"Fast det 2", tot:"tot det"}

            except :
                FD1 = prep_cube(speed_pref_dext, "FD1_CAR3")
                #grazdn.standard_name = "grazdn"
                print(FD1)
                FDS = prep_cube(speed_pref_ptrc, "DTC")
                FDS.data = FDS.data * 0.5
                #grazzmi.standard_name = "grazzmi"
                print(FDS)

                cube =  FDS.copy()
                #cube.long_name="Detritus flux"
                cube.data = FDS.data  + FD1.data

                #Sexp = FDS.data[0,23,:,:] + FD1.data[0,23,:,:]


            try :
                depth = cube.coord('Vertical T levels').points
            except :
                depth = cube.coord('depth').points
            for zz in np.arange(75) :
                cube.data[0,zz,:,:] = cube.data[0,zz,:,:] / cube.data[0,23,:,:]

            cube.long_name = "sinking detrical carbon flux"
            #depth = FDDT.coord('Vertical T levels').points
            #print(depth)
            if loc == "PAP" :
                jj = 258
                ii = 272
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "BATS" :
                jj = 237
                ii = 224
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL1" :
                jj = 85
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL2" :
                jj = 110
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL3" :
                jj = 150
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL4" :
                jj = 200
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL5" :
                jj = 260
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "GLOB" :
                area = cube.data.copy()
                for zz in np.arange(75) :
                    area[:,zz,:,:] = area2D[:,:,:]

                area.mask = cube.data.mask
                ppp = cube.data.copy()
                ppp1 = cube.data.copy()
                ppp1 = ppp * area

                ppp2 = ppp1.sum(axis = 2)
                ppp3 = ppp2.sum(axis = 2)

                area1 = area.sum(axis = 2)
                area2 = area1.sum(axis = 2)

                pocprof = ppp3 / area2

            if 'fig_ttl' in attrs:
                yr = attrs['fig_ttl']
            if (yr == "YY") and (Manage==True):
                ax.plot(pocprof[0,:],-depth[:],'k--', label = yr )
            else :
                ax.plot(pocprof[0,:],-depth[:],label = yr )


            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)
        ###
        # plot martin curve :
        if Martin==True :
            dd = depth.copy()
            FF = 1.53
            bb = -0.858
            M100 = FF * (100/100)**bb
            Mart = FF * (dd/100)**bb
            Mnorm = Mart / M100
            ax.plot(Mnorm,-depth[:],'k:',label = 'Martin' )
        ###
        ## manage plot
        ax.legend()
        if Manage==True:
            if LOGSCALE :
                ax.set_yscale('log')
            ax.set_ylim([-600, 0.0])
            #ax.set_xlim([0.0, 2.0])
            ax.set_xlabel(cube.units)
            ax.set_ylabel("Depth (m)")
            if 'fig_ttl' in attrs:
                fig_ttl = attrs['fig_ttl']
                ax.set_title(cube.long_name + " - " + fig_ttl)
            else :
                ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)










class subplot_NormPOC_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, run_nm ,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, LOGSCALE=True,
                 Martin=True, Manage=True, Norm=True, YLIM=[0.0, 5000], **attrs):

        """
        Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA

        meshfile="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc"
        e1t = prep_cube(meshfile, "e1t")
        e2t = prep_cube(meshfile, "e2t")
        area2D = e1t.data.copy()
        area2D = e1t.data * e2t.data
        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))
        if ax is None:
            ax = fig.add_subplot(iii, jjj, rect)
        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        if 'loc' in attrs:
            loc = attrs['loc']

        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            #speed_pref_diad = out_path + "/" + diadf(yr)
            #speed_pref_ptrc = out_path + "/" + ptrcf(yr)

            FDDT = prep_cube(rundict_diad[run_nm], "FD_CAR3")
            SDDT = prep_cube(rundict_ptrc[run_nm], "DTC")
            SDDT.data = SDDT.data * 0.5 ## 0.5 m/d
            SDDT.units = "mmol-C/m2/d"

            cube =  FDDT.copy()
            #expf100 = FDDT.extract(iris.Constraint(depth=100))
            #exps100 = SDDT.extract(iris.Constraint(depth=100))
            cube.data = (FDDT.data + SDDT.data)
            if Norm==True:
                Sexp = FDDT.data[0,23,:,:] + SDDT.data[0,23,:,:]
                #Stemp = cube.copy()
                #Sexp = Stemp.extract(iris.Constraint(depth=100))
                ##
                for zz in np.arange(75) :
                    cube.data[0,zz,:,:] = cube.data[0,zz,:,:] / Sexp.data[:,:]

            cube.long_name = "sinking detrical carbon flux"
            try :
                depth = FDDT.coord('Vertical T levels').points
            except :
                depth = FDDT.coord('depth').points
            #print(depth)
            if loc == "PAP" :
                jj = 258
                ii = 272
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "BATS" :
                jj = 237
                ii = 224
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL1" :
                jj = 85
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL2" :
                jj = 110
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL3" :
                jj = 150
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL4" :
                jj = 200
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "ATL5" :
                jj = 260
                ii = 262
                pocprof = cube.data[:,:,jj,ii].copy()
            elif loc == "GLOB" :
                area = cube.data.copy()
                for zz in np.arange(75) :
                    area[:,zz,:,:] = area2D[:,:,:]

                area.mask = cube.data.mask
                ppp = cube.data.copy()
                ppp1 = cube.data.copy()
                ppp1 = ppp * area

                ppp2 = ppp1.sum(axis = 2)
                ppp3 = ppp2.sum(axis = 2)

                area1 = area.sum(axis = 2)
                area2 = area1.sum(axis = 2)

                pocprof = ppp3 / area2

            if 'fig_lab' in attrs:
                lab = attrs['fig_lab']
            else :
                lab= None

            if (yr == "NO") and (Manage==True):
                ax.plot(pocprof[0,:],depth[:],'k--', label = lab )
            else :
                ax.plot(pocprof[0,:],depth[:],label = lab )


            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)
        ###
        # plot martin curve :
        if Martin==True:
            dd = depth.copy()
            FF = 1.53
            bb = -0.858
            M100 = FF * (100/100)**bb
            Mart = FF * (dd/100)**bb
            if Norm==True:
                Mnorm = Mart / M100
                ax.plot(Mnorm,depth[:],'k:',label = 'Martin' )
            else:
                ax.plot(Mart,depth[:],'k:',label = 'Martin' )
        ###
        ## manage plot
        ax.legend()
        if Manage==True:
            if LOGSCALE :
                ax.set_yscale('log')
            ax.set_ylim(YLIM)
            ax.invert_yaxis()
            if Norm==True:
                ax.set_xlim([0.0, 2.0])
                ax.set_xlabel("- normalized - ")
            elif 'xlim' in attrs:
                xlim = attrs['xlim']
                ax.set_xlim(xlim)
                ax.set_xlabel(FDDT.units)
            else:
                ax.set_xlabel(FDDT.units)

            ax.set_ylabel("Depth (m)")
            if 'fig_ttl' in attrs:
                fig_ttl = attrs['fig_ttl']
                ax.set_title(cube.long_name + " - " + fig_ttl)
            else :
                ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_xlabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)










class subplot_POC_monthprofiles(object):
    """
    Plot projected maps
    """

    def __init__(self, run_nm ,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, **attrs):

        """
        Make the monthly profile subplots
        """

        import matplotlib.pyplot as plt

        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA


        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))
        if ax is None:
            ax = fig.add_subplot(iii, jjj, rect)

        #fig.add_subplot(ax)
        #### Subplot defined - now plot :
        aaa = 9000.0
        bbb = 0.0

        #for yr in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"] :
        #for yr in ["01", "04", "07", "10", "YY"] :
        for yr in ["YY"] :

            #speed_pref_diad = out_path + "/" + diadf(yr)
            #speed_pref_ptrc = out_path + "/" + ptrcf(yr)

            FDDT = prep_cube(rundict_diad[run_nm], "FD_CAR3")
            SDDT = prep_cube(rundict_ptrc[run_nm], "DTC")
            SDDT.data = SDDT.data * 0.5 ## 0.5 m/d
            SDDT.units = "mmol-C/m2/d"
            cube =  FDDT.copy()
            cube.data = FDDT.data + SDDT.data
            cube.long_name = "sinking detrical carbon flux"
            try :
                depth = FDDT.coord('Vertical T levels').points
            except :
                depth = FDDT.coord('depth').points

            if 'fig_lab' in attrs:
                lab = attrs['fig_lab']
            else :
                lab= None

            if yr == "YY":
                ax.plot(cube.data[0,:,258, 272],-depth[:],'k--', label = lab )
            else :
                ax.plot(cube.data[0,:,258, 272],-depth[:],label = lab )

            aaa_l = cube.data.min()
            bbb_l = cube.data.max()
            aaa = np.min([aaa, aaa_l])
            bbb = np.max([bbb, bbb_l])

        print('min var = ',aaa,'; max var = ',bbb)
        ###
        # plot martin curve :
        dd = depth.copy()
        FF = 1.53
        bb = -0.858
        Mart = FF * (dd/100)**bb
        ax.plot(Mart,-depth[:],'k:',label = 'Martin' )
        ###
        ## manage plot
        ax.legend()
        ax.set_ylim([-600, 0.0])
        ax.set_xlim([0.0, 18.0])
        ax.set_xlabel(FDDT.units)
        ax.set_ylabel("Depth (m)")
        if 'fig_ttl' in attrs:
            fig_ttl = attrs['fig_ttl']
            ax.set_title(cube.long_name + " - " + fig_ttl)
        else :
            ax.set_title(cube.long_name)



        self.ax = ax                   # Graphical axes



    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_ylabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)









class subplot_to_put(object):
    """
    Plot (presumably) Obs vs model to put in already/outside defined subplot
    -- only plot the 2 first experiment of the list - not more --
    """

    def __init__(self, var, rundict, run_exp_list,
                 fig=None, ax = None, iii = 1, jjj = 2, rect = 1,
                 Meth="None", transect_name = None,
                 yr="None", log_scale = False, centered=False,
                 colorbar=True, location=None,
                 xrev = False, ytick_loc = "left",
                 title=None, **attrs) :

        print(var)

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
        runref = run_exp_list[0]
        ref_file_name = rundict[runlist[runref]]


        if var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL", "TPP3"] :
            Epipel_plot = True
        else :
            Epipel_plot = False

        #plt.rcParams.update({'font.size': 12})
        cmap0 = plt.cm.get_cmap('viridis', 15)
        cmap0 = plt.cm.get_cmap('turbo', 15)
        #cmap1 = plt.cm.get_cmap('RdBu', 15)
        cmap1 = plt.cm.get_cmap('RdBu_r', 15)
        cmap3 = plt.cm.get_cmap('RdBu', 15)
        cmap2 = plt.cm.get_cmap('viridis_r', 15)

        cmap0.set_bad(color='white')
        cmap1.set_bad(color='gray')
        cmap2.set_bad(color='gray')
        cmap3.set_bad(color='gray')

        if Meth == "None" :
            f_name = "Surf_plot_" + var
        else :
            f_name = "Surf_plot_" + var + "_" + Meth
        if yr != "None" :
            f_name = f_name + "_" + yr

        comp_nb = np.size(run_exp_list) - 1


        #fig = plt.figure(figsize=(5*jjj,5*iii))
        ## no pac
        cube = prep_cube(ref_file_name, var)
        print(cube)
        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Surface":
            try :
                cube = prep_cube_surf(cube)
            except :
                print("already 2D -- no need to extract ")
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ")
        elif Meth =="Transect" :
            if transect_name == None :
                print("---- Need to specify the transect name ----")
            else :
                if transect_name == "atlantic" :
                    tlat = 262
                elif transect_name == "pacific" :
                    tlat = 115
        ##
        #aaa = cube.data.min()
        #bbb = cube.data.max()
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            cube.data = ma.masked_where(cube.data == np.nan, cube.data )
            aaa = np.percentile(cube.data.compressed(), 2)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            cube.data = ma.masked_where(cube.data == np.nan, cube.data )
            bbb = np.percentile(cube.data.compressed(), 98)
        if centered :
            #ccc = (-aaa + bbb)/2
            ccc = np.max([-aaa, bbb])
            aaa = -ccc
            bbb = ccc
            cmap = cmap1
        else :
            cmap = cmap0
        ### plot :
        if Meth=="Transect":
            if 'Min' in attrs:
                aaa = attrs['Min']
            else :
                aaa = np.percentile(cube.data[0,:,:,tlat].compressed(), 2)
            if 'Max' in attrs:
                bbb = attrs['Max']
            else :
                bbb = np.percentile(cube.data[0,:,:,tlat].compressed(), 98)
            #
            plot = subplot_proj_Orcagrid(cube[0,:,:,tlat], fig = fig, ax=ax, iii=iii, jjj=jjj, rect=rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = colorbar, location=location, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, xrev=xrev, ytick_loc=ytick_loc )
        else:
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax, iii=iii, jjj=jjj, rect=rect, colorbar = colorbar, location=location, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale )

        if title is not None :
            ref_ttl=title
        else :
            #ref_ttl       = runlist[runref]
            ref_ttl       = cube.long_name

        plot.subtitle(ref_ttl)








class subplot_ObsMod(object):
    """
    Plot (presumably) Obs vs model to put in already/outside defined subplot
    -- only plot the 2 first experiment of the list - not more --
    """

    def __init__(self, var, rundict, run_exp_list,
                 fig=None, ax = None, ax1 = None, iii = 1, jjj = 2, rect = 1,
                 Meth="None", yr="None", log_scale = False, centered=False,
                 title=None, **attrs) :

        print(var)

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
        if title is not None :
            ref_ttl=title
        else :
            runref = run_exp_list[0]
            ref_file_name = rundict[runlist[runref]]
            ref_ttl       = runlist[runref]

        if var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL", "TPP3"] :
            Epipel_plot = True
        else :
            Epipel_plot = False

        #plt.rcParams.update({'font.size': 12})
        cmap0 = plt.cm.get_cmap('viridis', 15)
        cmap0 = plt.cm.get_cmap('turbo', 15)
        #cmap1 = plt.cm.get_cmap('RdBu', 15)
        cmap1 = plt.cm.get_cmap('RdBu_r', 15)
        cmap3 = plt.cm.get_cmap('RdBu', 15)
        cmap2 = plt.cm.get_cmap('viridis_r', 15)

        cmap0.set_bad(color='white')
        cmap1.set_bad(color='gray')
        cmap2.set_bad(color='gray')
        cmap3.set_bad(color='gray')

        if Meth == "None" :
            f_name = "Surf_plot_" + var
        else :
            f_name = "Surf_plot_" + var + "_" + Meth
        if yr != "None" :
            f_name = f_name + "_" + yr

        comp_nb = np.size(run_exp_list) - 1


        #fig = plt.figure(figsize=(5*jjj,5*iii))
        ## no pac
        cube = prep_cube(ref_file_name, var)
        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Surface":
            try :
                cube = prep_cube_surf(cube)
            except :
                print("already 2D -- no need to extract ")
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ")
        #aaa = cube.data.min()
        #bbb = cube.data.max()
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            cube.data = ma.masked_where(cube.data == np.nan, cube.data )
            aaa = np.percentile(cube.data.compressed(), 2)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            cube.data = ma.masked_where(cube.data == np.nan, cube.data )
            bbb = np.percentile(cube.data.compressed(), 98)
        if centered :
            #ccc = (-aaa + bbb)/2
            ccc = np.max([-aaa, bbb])
            aaa = -ccc
            bbb = ccc
            cmap = cmap1
        else :
            cmap = cmap0
        ### plot :
        if Meth=="Transect":
            if 'Min' in attrs:
                aaa = attrs['Min']
            else :
                aaa = np.percentile(cube.data[0,:,:,262].compressed(), 2)
            if 'Max' in attrs:
                bbb = attrs['Max']
            else :
                bbb = np.percentile(cube.data[0,:,:,262].compressed(), 98)
            #
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, ax=ax, iii=iii, jjj=jjj, rect=rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = False, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, location="left", remove_cb = 'whtvr' )
        else:
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax, iii=iii, jjj=jjj, rect=rect, colorbar = False, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, location="left" )
        if rect==1 :
            plot.subtitle("Observation")
        plot.y_label(cube.long_name)

        loop = 1
        compref = run_exp_list[1]
        comp_file_name = rundict[runlist[compref]]
        comp_ttl       = runlist[compref]

        # no_speed_pref
        cube = prep_cube(comp_file_name, var)
        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Surface":
            try :
                cube = prep_cube_surf(cube)
            except :
                print("already 2D -- no need to extract ")
                #cube = prep_cube_surf(cube)
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ")
        ### plot :
        #
        rect = rect+1
        ##
        if Meth=="Transect":
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect=rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = False, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, location="right")
        else:
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect=rect, colorbar = False, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, location="right")
        if rect==2 :
            plot.subtitle("PAC-model")









def subplot_short(var, rundict, run_exp_list, Meth="None", yr="None", **attrs) :

    print(var)

    # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    runref = run_exp_list[0]
    ref_file_name = rundict[runlist[runref]]
    ref_ttl       = runlist[runref]

    if var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL"] :
        Epipel_plot = True
    else :
        Epipel_plot = False

    plt.rcParams.update({'font.size': 12})
    cmap0 = plt.cm.get_cmap('viridis', 15)
    #cmap1 = plt.cm.get_cmap('RdBu', 15)
    cmap1 = plt.cm.get_cmap('RdBu_r', 15)
    cmap3 = plt.cm.get_cmap('RdBu', 15)
    cmap2 = plt.cm.get_cmap('viridis_r', 15)

    cmap0.set_bad(color='white')
    cmap1.set_bad(color='gray')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='gray')

    if Meth == "None" :
        f_name = "Surf_plot_" + var
    else :
        f_name = "Surf_plot_" + var + "_" + Meth
    if yr != "None" :
        f_name = f_name + "_" + yr

    comp_nb = np.size(run_exp_list) - 1
    ## prep subplot display
    if comp_nb == 1 :
        iii=2
        jjj=2
    elif comp_nb == 2 :
        iii=2
        jjj=3
    elif comp_nb == 3 :
        iii=2
        jjj=4
    elif 4 <= comp_nb <= 5 :
        iii=3
        jjj=4
    elif 6 <= comp_nb <= 7 :
        iii=4
        jjj=4
    elif 8 <= comp_nb <= 11 :
        iii=4
        jjj=6


    fig = plt.figure(figsize=(5*jjj,5*iii))
    ## no pac
    cube = prep_cube(ref_file_name, var)
    if Meth=="Surface":
        try :
            cube = prep_cube_surf(cube)
        except :
            print("already 2D -- no need to extract ")
    elif Meth=="Vert_Inv":
        cube = prep_cube_vert_inv(cube)
    #aaa = cube.data.min()
    #bbb = cube.data.max()
    if 'Min' in attrs:
        aaa = attrs['Min']
    else :
        cube.data = ma.masked_where(cube.data == np.nan, cube.data )
        aaa = np.percentile(cube.data.compressed(), 1)
    if 'Max' in attrs:
        bbb = attrs['Max']
    else :
        cube.data = ma.masked_where(cube.data == np.nan, cube.data )
        bbb = np.percentile(cube.data.compressed(), 99)
    ### plot :
    if Meth=="Transect":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[0,:,:,262].compressed(), 1)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[0,:,:,262].compressed(), 99)
        #
        print("ploting transect between min = ", Min," and max = ",Max)
        plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb )
    else:
        print("ploting flat between min = ", Min," and max = ",Max)
        plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = 1,colorbar = True, v_min= aaa,v_max=bbb )
    plot.subtitle(ref_ttl)


    plt.suptitle(cube.long_name, size=16)
    plt.savefig(f_name)
    plt.show()








def subplot_diff_list(var, rundict_list, loc_list, filetype="ptrc", Meth="None", yr="None", log_scale = False, outfreq = "1M", centered=False, region=None, Epipel=False, SHOW=False, **attrs) :
               #(var, rundict, run_exp_list, Meth="None", yr="None", log_scale = False, centered=False, region=None, **attrs) :

    print(var)

    # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    runref = loc_list[0]
    print("loc_list : ",loc_list)

    runnb = rundict_list[0]
            #
            #ax1 = fig.add_subplot(iii, jjj, rect)
            #
    modclass = runclass(runnb)
            #
    print(filetype)
    if filetype=="ptrc" :
        rundict        = modclass.dict_ptrc
    elif filetype=="diad" :
        rundict        = modclass.dict_diad
    elif filetype=="grid" :
        rundict        = modclass.dict_grid
    print(rundict)
    ref_file_name = rundict[runlist[runref]]
    ref_ttl       = modclass.name


    if var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL", "TPP3"] :
        Epipel_plot = True
    else :
        Epipel_plot = Epipel

    if var in ["CHL", "PP", "CHN", "CHD", "PHN", "PHD", "DIN"] and Meth=="Surface":
        log_scale = True

    plt.rcParams.update({'font.size': 12})
    cmap0 = plt.cm.get_cmap('viridis', 15)
    cmap0 = plt.cm.get_cmap('turbo', 15)
    #cmap1 = plt.cm.get_cmap('RdBu', 15)
    cmap1 = plt.cm.get_cmap('RdBu_r', 15)
    cmap3 = plt.cm.get_cmap('RdBu', 15)
    cmap2 = plt.cm.get_cmap('viridis_r', 15)

    cmap0.set_bad(color='white')
    cmap1.set_bad(color='gray')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='gray')

    if Meth == "None" :
        f_name = "Surf_plot_" + var
    else :
        f_name = "Surf_plot_" + var + "_" + Meth
    if yr != "None" :
        f_name = f_name + "_" + yr
    if region is not None :
        f_name = f_name + "_" + region

    ## regional plot :
    if region == "Amazon" :
        jj_min = 175
        jj_max = 220
        ii_min = 220
        ii_max = 260
    elif region == "Bengal" :
        jj_min = 200
        jj_max = 230
        ii_min = 0
        ii_max = 30
    elif region == "UKshelv" :
        jj_min = 250
        jj_max = 280
        ii_min = 275
        ii_max = 295
        #ii_min=-7
        #ii_max=2
        #jj_min=45
        #jj_max=55

    comp_nb = (np.size(loc_list)*np.size(rundict_list)) - 1
    ## prep subplot display
    if comp_nb == 0 :
        iii=1
        jjj=1
    elif comp_nb == 1 :
        iii=2
        jjj=2
    elif comp_nb == 2 :
        iii=2
        jjj=3
    elif comp_nb == 3 :
        iii=2
        jjj=4
    elif 4 <= comp_nb <= 5 :
        iii=3
        jjj=4
    elif 6 <= comp_nb <= 7 :
        iii=4
        jjj=4
    elif 8 <= comp_nb <= 11 :
        iii=4
        jjj=6


    fig = plt.figure(figsize=(5*jjj,5*iii))
    #fig = plt.figure(figsize=(5*6,5*6))
    ## no pac

    rect=1

    if Meth=="OneD" or Meth=="surface_OneD" :
        oneD=True
    else :
        oneD=False

      ## no pac
    cube = prep_cube(ref_file_name, var, oneD=oneD)
                #
    if region is not None :
        #c_coord = prep_cube("../VESTA/Out/medusa_cn857o_1y_20301201-20311201_ptrc-T.nc","SIL")
        #longit = c_coord.coord('longitude').points
        #latit = c_coord.coord('latitude').points
        longitt = read_cube("/noc/users/jpp1m13/WORKING/UKESM/","nav_lon", mask = False)
        latitt  = read_cube("/noc/users/jpp1m13/WORKING/UKESM/","nav_lat", mask = False)
        cube.coord('latitude').points = latitt.data
        cube.coord('longitude').points = longitt.data
        longit = longitt[jj_min:jj_max,ii_min:ii_max]
        latit  = latitt[jj_min:jj_max,ii_min:ii_max]

    if Meth=="Surface" or Meth=="surface_OneD" :
        cube = prep_cube_surf(cube, oneD=oneD)
    elif Meth=="Vert_Inv":
        cube = prep_cube_vert_inv(cube)
        #aaa = cube.data.min()
        #bbb = cube.data.max()
    if 'Min' in attrs:
        aaa = attrs['Min']
    else :
        aaa = np.percentile(cube.data.compressed(), 5)
    if 'Max' in attrs:
        bbb = attrs['Max']
    else :
        bbb = np.percentile(cube.data.compressed(), 95)
    if centered :
        #ccc = (-aaa + bbb)/2
        ccc = np.max([-aaa, bbb])
        aaa = -ccc
        bbb = ccc
        cmap = cmap1
    else :
        cmap = cmap0
    ### plot :
    if Meth=="Transect":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[0,:,:,262].compressed(), 5)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[0,:,:,262].compressed(), 95)
        #
        plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect = rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    elif Meth=="OneD" or Meth=="surface_OneD":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[...,0,0].compressed(), 5)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[...,0,0].compressed(), 95)
        #
        plot = subplot_timeseries(cube[...,0,0], fig = fig, iii=iii, jjj=jjj, rect = rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=modclass.outfreq, log_scale=log_scale )
        if cube.ndim == 4 :
            plot.y_label("Depth (m)")
        elif cube.ndim < 4 :
            plot.y_label(cube.units)
            plot.x_label("time in Year")
    else:
        plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = rect, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    plot.subtitle(modclass.name)

    loop = 1

    for runnb in rundict_list[1:] :

            #runnb = rundict_list[0]
            #
            #ax1 = fig.add_subplot(iii, jjj, rect)
            #
        modclass = runclass(runnb)
            #
        print(filetype)
        if filetype=="ptrc" :
            rundict        = modclass.dict_ptrc
        elif filetype=="diad" :
            rundict        = modclass.dict_diad
        elif filetype=="grid" :
            rundict        = modclass.dict_grid
        print(rundict)
        comp_file_name = rundict[runlist[runref]]
        comp_ttl       = modclass.name


        # no_speed_pref
        cube = prep_cube(comp_file_name, var, oneD=oneD)

        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Surface" or Meth=="surface_OneD" :
            cube = prep_cube_surf(cube, oneD=oneD)
        elif Meth=="Vert_Inv":
            cube = prep_cube_vert_inv(cube)
        ### plot :
        if comp_nb < 3 :
            rect = loop+1
        else :
            rect = (2 * loop)+1
        ##
        if Meth=="Transect":
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
        elif Meth=="OneD" or Meth=="surface_OneD":
            plot = subplot_timeseries(cube[...,0,0], fig = fig, iii=iii, jjj=jjj, rect = rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=modclass.outfreq, log_scale=log_scale )
            if cube.ndim == 4 :
                plot.y_label("Depth (m)")
            elif cube.ndim < 4 :
                plot.y_label(cube.units)
                plot.x_label("time in Year")
        elif region is None :
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale )
        else :
            plot = subplot_proj_regional(cube[0,jj_min:jj_max,ii_min:ii_max], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, proj=region, longit=longit.data, latit=latit.data )
        plot.subtitle(modclass.name)
        #plot.subtitle(comp_ttl)
        #####
        ##
        # Diff :
        cube = prep_cube_diff(ref_file_name, comp_file_name, var, oneD=oneD)
        #if region is None :
        #    cube = prep_cube_diff(ref_file_name, comp_file_name, var)
        #else :
        #    cube_tmp = prep_cube_diff(ref_file_name, comp_file_name, var)
        #    cube = cube_tmp[...,jj_min:jj_max,ii_min:ii_max].copy()

        if Meth=="Surface"or Meth=="surface_OneD" :
            cube = prep_cube_surf(cube, oneD=oneD)
        elif Meth=="Vert_Inv":
            cube = prep_cube_vert_inv(cube)
        ### get diff min max
        if loop == 1 :
            #aaa_d = cube.data.min()
            #bbb_d = cube.data.max()
            #ccc = np.max([- aaa_d, bbb_d])
            if 'Min_diff' in attrs:
                aaa_d = attrs['Min_diff']
            elif Meth=="OneD" or Meth=="surface_OneD":
                    aaa_d = np.percentile(cube.data[...,0,0].compressed(), 5)
            else :

                cube.data = ma.masked_where(cube.data == np.nan, cube.data )
                if region is None :
                    aaa_d = np.percentile(cube.data.compressed(), 2)
                else :
                    aaa_d = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 2)

            if 'Max_diff' in attrs:
                bbb_d = attrs['Max_diff']
            elif Meth=="OneD" or Meth=="surface_OneD":
                bbb_d = np.percentile(cube.data[...,0,0].compressed(), 95)
            else :
                cube.data = ma.masked_where(cube.data == np.nan, cube.data )
                if region is None :
                    bbb_d = np.percentile(cube.data.compressed(), 98)
                else :
                    bbb_d = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 98)
            #
            #ccc = (-aaa_d + bbb_d)/2
            ccc = np.max([-aaa_d, bbb_d])
        ### plot :
        if comp_nb < 3 :
            rect = loop + comp_nb + 2
        else :
            rect = (2 * loop) + 2
        ##

        if Meth=="Transect":
            if loop == 1 :
                if 'Min_diff' in attrs:
                    aaa_d = attrs['Min_diff']
                else :
                    aaa_d = np.percentile(cube.data[0,:,:,262].compressed(), 2)
                if 'Max_diff' in attrs:
                    bbb_d = attrs['Max_diff']
                else :
                    bbb_d = np.percentile(cube.data[0,:,:,262].compressed(), 98)
                #
                ccc = np.max([- aaa_d, bbb_d])
                #ccc = (-aaa_d + bbb_d)/2
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        elif Meth=="OneD" or Meth=="surface_OneD":
            plot = subplot_timeseries(cube[...,0,0], fig = fig, iii=iii, jjj=jjj, rect = rect, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc, outfreq=modclass.outfreq, log_scale=log_scale )
            if cube.ndim == 4 :
                plot.y_label("Depth (m)")
            elif cube.ndim < 4 :
                plot.y_label(cube.units)
                plot.x_label("time in Year")
        elif region is None :
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        else :
            plot = subplot_proj_regional(cube[0,jj_min:jj_max,ii_min:ii_max], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc, proj=region, longit=longit.data, latit=latit.data)

        plot.subtitle(comp_ttl + " minus " + ref_ttl)
        # ++ loop
        loop = loop + 1
    ## End loop


    plt.suptitle(cube.long_name, size=16)
    plt.rcParams.update({'font.size': 28})
    plt.savefig(f_name)
    if SHOW :
        plt.show()











def subplot_diff(var, rundict, runlist, run_exp_list, Meth="None", depth_extract = 100, yr="None", Regrid = False, NEW_regrid=False, MASK=None, MASK_ZERO=False, SINGLE_PLOTS=False, COAST_LINE=True, log_scale = False, EPIPEL=None, centered=False, region=None, proj=None, Obs= False, meshmask="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/eORCA100_masks.nc", Convert=1.0, ConvUnit = "None", LOC="Atl", Titre=None, REF_COMP=None, SHOW=False, **attrs) :
    '''
    Subplot Diff :
    called as subplot_diff(var,
                           rundict,
                           run_exp_list,
                           Meth="None",
                           depth_extract = 100,
                           yr="None",
                           Regrid = False,
                           MASK=None,
                           MASK_ZERO=False,
                           SINGLE_PLOTS=False,
                           COAST_LINE=True,
                           log_scale = False,
                           EPIPEL = None,
                           centered=False,
                           region=None, 
                           proj=None, 
                           Obs= False, 
                           meshmask="../MESH/eORCA1/eORCA100_masks.nc",
                           Convert=1.0,
                           ConvUnit = "None",
                           LOC="Atl",
                           Titre=None,
                           REF_COMP=None,
                           **attrs)                           
    With :
    -- Var the variables to plot as in the netcdf files
    -- rundict : dictionary with all files we could possibly compare
    -- run_exp_list : number of the experiments from rundict to compare. the first one being the reference.
    -- Meth could be :
            "Surface"    : extract surface 2D field if not already 2D
            "Vert_Inv"   : perform a vertical sum
            "2D_Extract" : extract a specific depth -- triggers depth_extract (1 ; 100; 200; 500 and 1000 m depth)
            "Long_Avg"   : perform avg along longitude (results on depth * lat array)
            "Lat_Avg"    : perform avg along latitude (results on depth * long array)
            "Atl_Pac_Lat": Perform Long avg in Atl and Pac ocean and stick them in one plot.
                           overide Mask.
            "Atl_Pac_hov": Perform surface Atl_Pac_Lat average for hovmoller plot.
    -- depth_extract : in m ; used if Meth is "2D_extract; could be 1, 100, 200, 4500, 1000
    -- Regrid : set to True to regrid on 1x1 regular (WOA like) grid -- regridded files are masked with MASK, and saved.
                reggrided files are read from previous saved file if it exists.
                Renew the calculation with NEW_regrid=True
    -- NEW_regrid : don't read the saved regrided file, re-do the regrid and overwrite the saved nc for future use.
    -- MASK : apply a Mask if not None. could be set as :
          wor_msk  ;
          atl_msk  ;
          pac_msk  ;
          ind_msk  ;
          sou_msk  ;
          arc_msk  ;
          all_atl  ;
          all_pac  ;
          all_ind  ;
          nor_atl  ;
          sou_atl  ;
          nor_pac  ;
          sou_pac  ;
          atl_sou  ;
          pac_sou  ;
          ind_sou  ;
          msk_100  ;
          msk_200  ;
    -- MASK_ZERO : mask 0 values -- Masks are not based on 0 values otherwise, which is good, but one might want that sometimes
    -- SINGLE_PLOTS : produce all plots individually and not as joined subplots.
    -- COAST_LINE   : Add the coastline on the map.
    -- EPIPEL : set to True to plot sections or profiles down to 500m depth.
                if set to False will plot the whole vertical depth
                if left as None, will automatically set as True for var in "ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL", "TPP3"
    -- log_scale : Set to True to plot with log colorbar
    -- centered  : Set to True to plot with centered colorbar (True for the difference plots)
    -- region : for specific regional plots or analysis. can be :
               "Amazon", "Bengal", "UKshelv", "NorthOrtho" or "NorthPolarStereo" 
    -- proj : possible plot projections; to be chosen between: 
            'Mollweide','Robinson','Mollweide_pac','Mollweide_ind','PlateCarree',
            'NorthPolarStereo','SouthPolarStereo','Orthographic',
            SouthOrtho','NorthOrtho','OrthoAtl','OrthoPac','OrthoInd',
            'Amazon','Bengal,'UKshelv','PlateCarree'
    -- Obs : to be set as True if the file to read is not standard NEMO outputs
    -- meshmask : path and name of the mesh mask nc file
    -- Convert : convertion Factor (default is 1.0)
    -- ConvUnit : Resulting unit once converted (used if not None)
    -- LOC : used if Meth is "transect". could be :
        "Atl" for an Atlantic section
        "Pac" for a Pacific section
    -- Titre : to put a specific (Non default) title onto the final plot
    -- REF_COMP : Only used for Alk efficiency calculation.
    -- SHOW    :  print the plot on screen at run-time if True. (stop the script from running)
    -- Extra Optional Arguments to put at the end :
        Min (Minimum value for Colorbar)
        Max (Maximum value for Colorbar)
        Min_diff (Minimum value for comparison plot Colorbar)
        Max_diff (Minimum value for comparison plot Colorbar)
    '''
    ## depth_extract :  choices are :
    ##            1 ; 100; 200; 500 and 1000 m depth
    ##
    print(var)

    # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    runref = run_exp_list[0]
    ref_file_name = rundict[runlist[runref]]
    ref_ttl       = runlist[runref]
    ##
    if REF_COMP is not None :
        runcomp = REF_COMP
        ref_comp_fname = rundict[runcomp]
        comp_ttl       = runcomp
    else :
        ref_comp_fname = None
    ##
    if EPIPEL is not None:
        Epipel_plot = EPIPEL
    elif var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL", "TPP3"] :
        Epipel_plot = True
    else :
        Epipel_plot = False

    #if var in ["CHL", "PP", "CHN", "CHD", "PHN", "PHD", "DIN"] and Meth=="Surface":
    #    log_scale = True

    plt.rcParams.update({'font.size': 12})
    #cmap0 = plt.cm.get_cmap('viridis', 15)
    cmap0 = plt.cm.get_cmap('turbo', 30)
    #cmap1 = plt.cm.get_cmap('RdBu', 15)
    cmap1 = plt.cm.get_cmap('RdBu_r', 30)
    cmap3 = plt.cm.get_cmap('RdBu', 30)
    cmap2 = plt.cm.get_cmap('viridis_r', 30)

    cmap0.set_bad(color='white')
    cmap1.set_bad(color='white')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='white')

    if Meth == "None" :
        f_name = "plot_" + var
    elif Meth == "Surface" :
        f_name = "plot_" + var + "_Surf_"
    elif Meth == "2D_Extract" :
        f_name = "plot_" + var + "_" + str(depth_extract) +"m"
    else :
        f_name = "plot_" + var + "_" + Meth
    if yr != "None" :
        f_name = f_name + "_" + yr
    if Regrid:
        f_name = f_name + "_Regridded"
    if MASK is not None :
        f_name = f_name + "_" + MASK
    if region is not None :
        f_name = f_name + "_" + region
    if proj is not None :
        f_name = f_name + "_" + proj

    ## regional plot :
    if region == "Amazon" :
        jj_min = 175
        jj_max = 220
        ii_min = 220
        ii_max = 260
    elif region == "Bengal" :
        jj_min = 200
        jj_max = 230
        ii_min = 0
        ii_max = 30
    elif region == "UKshelv" :
        jj_min = 250
        jj_max = 280
        ii_min = 275
        ii_max = 295
        #ii_min=-7
        #ii_max=2
        #jj_min=45
        #jj_max=55
    elif region == "NorthOrtho" or region == "NorthPolarStereo":
        jj_min = 275
        jj_max = 332
        ii_min = 0
        ii_max = 361

    comp_nb = np.size(run_exp_list) - 1
    ## prep subplot display
    if (comp_nb == 0) or SINGLE_PLOTS :
        iii=1
        jjj=1
    elif comp_nb == 1 :
        iii=2
        jjj=2
    elif comp_nb == 2 :
        iii=2
        jjj=3
    elif comp_nb == 3 :
        iii=2
        jjj=4
    elif 4 <= comp_nb <= 5 :
        iii=3
        jjj=4
    elif 6 <= comp_nb <= 7 :
        iii=4
        jjj=4
    elif 8 <= comp_nb <= 11 :
        iii=4
        jjj=6


    ##
    ## define fig with size
    if SINGLE_PLOTS :
        fig = None
    elif Meth == "Atl_Pac_Lat" :
        fig = plt.figure(figsize=(20*jjj,10*iii))
    elif Meth == "Atl_Pac_hov" :
        fig = plt.figure(figsize=(10*jjj,20*iii))
    else :
        fig = plt.figure(figsize=(10*jjj,10*iii))
    #fig = plt.figure(figsize=(5*6,5*6))
    ##
    ## Want the colorbar on the right side for singlr plots :
    if SINGLE_PLOTS :
        location = "right"
    else :
        location = None

    ## Force Mask zero to be false for reggriding (masks everything otherwise-- need to check why): 
    if "Atl_Pac" in Meth or Regrid == True:
        MASK_ZERO=False
    ##
    ##
    ## Start :
    ## Load the first (ref) var :
    if Obs is True : 
        ref_comp_fname = rundict[runlist[1]]
        #
        ## need the diad file :
        ref_comp_fname = diad_replace(ref_comp_fname)
        ##
    if Meth=="Atl_Pac_hov":
        cube = prep_cube_Atl_Pac_hov_from_monthly(ref_file_name, var, runlist, run_nb=runref, REF_COMP=ref_comp_fname, Obs=Obs,MASK_ZERO=MASK_ZERO ,NEW_regrid=NEW_regrid)
        ##
    else:
        cube = prep_cube(ref_file_name, var, REF_COMP=ref_comp_fname, Obs=Obs,MASK_ZERO=MASK_ZERO )
    
    if Obs is True : 
        # need to re-init ref_comp_fname or it will break the rest of the code.
        ref_comp_fname = None
    ##
    if region is not None :
        try:
            longitt = read_cube(meshfile,"nav_lon", mask = False)
            latitt  = read_cube(meshfile,"nav_lat", mask = False)
            cube.coord('latitude').points = latitt.data
            cube.coord('longitude').points = longitt.data
            if region == "NorthOrtho" or region == "NorthPolarStereo":
                longit = longitt[jj_min:jj_max,ii_min:ii_max]
                latit  = latitt[jj_min:jj_max,ii_min:ii_max]
            else :
                longit = longitt[jj_min+1:jj_max+1,ii_min+1:ii_max+1]
                latit  = latitt[jj_min:jj_max,ii_min:ii_max]
        except :
            print(" !!! no lat-long attribute in the file... try without")
            longitt = read_cube(meshmask,"nav_lon", mask = False)
            latitt  = read_cube(meshmask,"nav_lat", mask = False)
            #cube.coord('latitude').points = latitt.data
            #cube.coord('longitude').points = longitt.data
            if region == "NorthOrtho" or region == "NorthPolarStereo":
                longit = longitt[jj_min:jj_max,ii_min:ii_max]
                latit  = latitt[jj_min:jj_max,ii_min:ii_max]
            else :
                longit = longitt[jj_min+1:jj_max+1,ii_min+1:ii_max+1]
                latit  = latitt[jj_min:jj_max,ii_min:ii_max]

    if "Atl_Pac" not in Meth:
        ## check here to fix the NEMO4.2 log/lat not being proprely defined problem...
        try :
            llat = cube.coord('latitude').points
            ## if no error raised, file is OK
        except :
            ## if not OK, file needs to be fixed :
            cube = fix_nemo_cube(cube)
            ## pas completement sur que ca marche...
            ##llat = cube.coord('latitude').points
    ##
    ##
    if Meth=="Surface":
        try :
            cube = prep_cube_surf(cube)
            if cube is None :
                print("cube ref surface gives None")
        except :
            print("already 2D -- ", cube.data.shape ," no need to extract ")
    ##
    if Meth != "Atl_Pac_Lat" or Meth != "Atl_Pac_hov":
        if MASK is not None :
            cube = prep_Masked_cube(cube,MASK)
            if cube is None :
                print("cube ref Masking gives None")
        ##
        if Regrid :
            print("Regridding cube_ref")
            cube = regrid_ORCA(cube, runlist, Meth=Meth, run_nb=runref, MASK=MASK, NEW=NEW_regrid)
            if cube is None :
                print("cube ref Regrid gives None")
            print("Regridded cube : ", cube.shape)
            #print("cube_ref : ", cube_ref.shape)
    ##
    ### if needed convertion :
    cube.data = cube.data * Convert
    if ConvUnit != "None" :
        cube.units = ConvUnit
    ###
    if log_scale :
        ## Make sure there are no neg values
        cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
    if Meth=="Atl_Pac_Lat":
        ## Does mask, regrid and stick atl pand PAC
        try :
            cube = prep_cube_Atl_Pac(cube, runlist, run_nb=runref, NEW=NEW_regrid)
            if cube is None :
                print("cube ref Atl_Pac_Lat gives None")
        except :
            print("Atl-Pac average not working... -- ", cube.data.shape )
    elif Meth=="Vert_Inv":
        try :
            cube = prep_cube_vert_inv(cube)
        except :
            print(" already 2D var -- ", cube.data.shape )
    elif Meth=="2D_Extract":
        try :
            cube = prep_cube_extract_depth(cube, depth_extract)
        except :
            print(" already 2D var -- ", cube.data.shape )
    elif Meth=="Long_Avg" :
        try :
            temp = cube.copy()
            cube = temp.collapsed('longitude', iris.analysis.MEAN)
            print("Average along longitude :")
            print("Long avg cube : ", cube.shape)
        except :
            print("Failed Average along longitude :")
            print(cube)
    elif Meth=="Lat_Avg" :
        try :
            temp = cube.copy()
            cube = temp.collapsed('latitude', iris.analysis.MEAN)
            print("Average along latitude :")
            print("Latg avg cube : ", cube.shape)
        except :
            print("Failed Average along latitude :")
            print(cube)
    #aaa = cube.data.min()
    #bbb = cube.data.max()
    if 'Min' in attrs:
        aaa = attrs['Min']
    else :
        cube.data = ma.masked_where(cube.data == np.nan, cube.data )
        if region is None :
            try :
                aaa = np.percentile(cube.data.compressed(), 2)
            except :
                aaa = np.min(cube.data)
        else :
            aaa = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 2)
    if 'Max' in attrs:
        bbb = attrs['Max']
    else :
        cube.data = ma.masked_where(cube.data == np.nan, cube.data )
        if region is None :
            try :
                bbb = np.percentile(cube.data.compressed(), 98)
            except :
                bbb = np.max(cube.data)
        else :
            bbb = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 98)
    if centered :
        #ccc = (-aaa + bbb)/2
        ccc = np.max([-aaa, bbb])
        aaa = -ccc
        bbb = ccc
        cmap = cmap1
    else :
        cmap = cmap0
    ### plot :
    if Meth=="Transect":
        if LOC == "Atl" :
            loci = 262
        elif LOC == "Pac" :
            loci = 112

        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[0,:,:,loci].compressed(), 2)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[0,:,:,loci].compressed(), 98)
        #
        plot = subplot_proj_Orcagrid(cube[0,:,:,loci], fig = fig, iii=iii, jjj=jjj, rect = 1,Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale )
    elif ( Meth=="Long_Avg") or (Meth=="Lat_Avg")  :
        plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=1, Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
    elif (Meth=="Atl_Pac_Lat")  :
        plot = subplot_proj_Orcagrid(cube[0,:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=1, Proj = False, location=location, Epipel_plot = Epipel_plot, Atl_Pac_section = True, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
    elif (Meth=="Atl_Pac_hov")  :
        plot = subplot_proj_Orcagrid(cube[:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=1, Proj = False, location=location, Epipel_plot = False, Atl_Pac_section = True, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
    else:
        if region is None :
            if proj is None :
                plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = 1,colorbar = True, COAST_LINE=COAST_LINE, location=location, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale )
            else :
                plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = 1,colorbar = True, COAST_LINE=COAST_LINE, location=location, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, proj=proj )
        else :
            plot = subplot_proj_regional(cube[0,jj_min:jj_max,ii_min:ii_max], fig = fig, iii=iii, jjj=jjj, rect = 1,colorbar = True, COAST_LINE=COAST_LINE, location=location, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, proj=region, longit=longit.data, latit=latit.data )
    if SINGLE_PLOTS :
        if (Meth=="Transect") or (Meth=="Long_Avg") or (Meth=="Atl_Pac_Lat") :
            plot.y_label("Depth [m]")
            plot.x_label("Latitude North [Degree]")
        elif (Meth=="Atl_Pac_hov") :
            plot.y_label("Latitude North [Degree]")
            #plot.x_label("Time (month)")
        elif  Meth=="Lat_Avg" :
            plot.y_label("Depth [m]")
            plot.x_label("Longitude East [Degree]")
        ##
        plt.savefig(f_name + ref_ttl, bbox_inches='tight', dpi=300)
        if SHOW :
            plt.show()
    else :
        plot.subtitle(ref_ttl)

    ### keep ref_cube
    cube_ref = cube.copy()
    ##############################################

    loop = 1
    for compref in run_exp_list[1:] :
        comp_file_name = rundict[runlist[compref]]
        comp_ttl       = runlist[compref]

        # no_speed_pref
        if Meth=="Atl_Pac_hov":
            cube = prep_cube_Atl_Pac_hov_from_monthly(comp_file_name, var, runlist, run_nb=compref, REF_COMP=ref_comp_fname, MASK_ZERO=MASK_ZERO ,NEW_regrid=NEW_regrid)
            ##
        else :
            cube = prep_cube(comp_file_name, var, REF_COMP=ref_comp_fname,MASK_ZERO=MASK_ZERO)
        ##
        #### Check anf fix NEMO4.2 nc file nc coord problem.N
        if "Atl_Pac" not in Meth:
            try :
                llat = cube.coord('latitude').points
                ## if no error raised, file is OK
            except :
                ## if not OK, file needs to be fixed :
                cube = fix_nemo_cube(cube)
                ## pas completement sur que ca marche...
                ##llat = cube.coord('latitude').points
        ##
        if Meth=="Surface":
            try :
                cube = prep_cube_surf(cube)
                if cube is None :
                    print("cube comp Surface gives None")
            except :
                print("already 2D -", cube.data.shape ,"- no need to extract ")
                #cube = prep_cube_surf(cube)
        #
        if Meth != "Atl_Pac_Lat" or Meth != "Atl_Pac_hov":
            if MASK is not None :
                cube = prep_Masked_cube(cube,MASK)
                if cube is None :
                    print("cube comp Masking gives None")
            ##
            if Regrid :
                print("Regridding cube_ref")
                cube = regrid_ORCA(cube, runlist, Meth=Meth, run_nb=compref, MASK=MASK, NEW=NEW_regrid)
                if cube is None :
                    print("cube comp Regrid gives None")
                print("cube : ", cube.shape)
                #print("cube_ref : ", cube_ref.shape)
        ##
        ### if needed convertion :
        cube.data = cube.data * Convert
        if ConvUnit != "None" :
            cube.units = ConvUnit
        ###
        if log_scale :
            ## Make sure there are no neg values
            cube.data = ma.masked_where(cube.data <= 0.0, cube.data )
        if Meth=="Atl_Pac_Lat":
            try :
                cube = prep_cube_Atl_Pac(cube, runlist, run_nb=compref, NEW=NEW_regrid)
                if cube is None :
                    print("cube comp Atl_Pac_Lat gives None")
            except :
                print("Atl-Pac average not working... -- ", cube.data.shape )
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ", cube.data.shape )
        elif Meth=="2D_Extract":
            try :
                cube = prep_cube_extract_depth(cube, depth_extract)
            except :
                print(" already 2D var -- ", cube.data.shape )
        elif Meth=="Long_Avg" :
            try :
                temp = cube.copy()
                cube = temp.collapsed('longitude', iris.analysis.MEAN)
                print("Average along longitude :")
                print(cube)
            except :
                print(" Failed Average along longitude :")
                print(cube)
        elif Meth=="Lat_Avg" :
            try :
                temp = cube.copy()
                cube = temp.collapsed('latitude', iris.analysis.MEAN)
                print("Average along latitude :")
                print(cube)
            except :
                print(" Failed Average along latitude :")
                print(cube)
        ### plot :
        if SINGLE_PLOTS :
            rect = 1
        elif comp_nb < 3 :
            rect = loop+1
        else :
            rect = (2 * loop)+1
        ##
        if Meth=="Transect":
            plot = subplot_proj_Orcagrid(cube[0,:,:,loci], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
        elif ( Meth=="Long_Avg") or (Meth=="Lat_Avg") :
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
        elif (Meth=="Atl_Pac_Lat") :
            plot = subplot_proj_Orcagrid(cube[0,:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, Atl_Pac_section = True, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
        elif (Meth=="Atl_Pac_hov") :
            plot = subplot_proj_Orcagrid(cube[:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = False, Atl_Pac_section = True, colorbar = True, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale)
        else:
            if region is None :
                if proj is None :
                    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location, COAST_LINE=COAST_LINE, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale )
                else :
                    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location, COAST_LINE=COAST_LINE, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, proj=proj )
            else :
                plot = subplot_proj_regional(cube[0,jj_min:jj_max,ii_min:ii_max], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location, COAST_LINE=COAST_LINE, v_min= aaa,v_max=bbb, cmap=cmap, log_scale=log_scale, proj=region, longit=longit.data, latit=latit.data )
        ##
        if SINGLE_PLOTS :
            if (Meth=="Transect") or (Meth=="Long_Avg") or (Meth=="Atl_Pac_Lat") :
                plot.y_label("Depth [m]")
                plot.x_label("Latitude North [Degree]")
            elif (Meth=="Atl_Pac_hov") :
                plot.y_label("Latitude North [Degree]")
                #plot.x_label("Time (month)")
            elif  Meth=="Lat_Avg" :
                plot.y_label("Depth [m]")
                plot.x_label("Longitude East [Degree]")
            ##
            plt.savefig(f_name + comp_ttl, bbox_inches='tight', dpi=300)
            if SHOW :
                plt.show()
        else :
            plot.subtitle(comp_ttl)
        # Diff :
        ## keep comp cube
        cube_comp = cube.copy()
        #####################################

        # Diff :
        #cube = prep_cube_diff(ref_file_name, comp_file_name, var)
        cube.data = cube_comp.data - cube_ref.data
        ### if needed convertion :
        #cube.data = cube.data * Convert
        #if ConvUnit != "None" :
        #    cube.units = ConvUnit
        ###
        if Meth=="Surface":
            cube = prep_cube_surf(cube)
        elif Meth=="Vert_Inv":
            try :
                cube = prep_cube_vert_inv(cube)
            except :
                print(" already 2D var -- ")
        elif Meth=="2D_Extract":
            try :
                cube = prep_cube_extract_depth(cube, depth_extract)
            except :
                print(" already 2D var -- ", FDDT.data.shape )
        ### get diff min max
        if loop == 1 :
            #aaa_d = cube.data.min()
            #bbb_d = cube.data.max()
            #ccc = np.max([- aaa_d, bbb_d])
            if 'Min_diff' in attrs:
                aaa_d = attrs['Min_diff']
            else :
                cube.data = ma.masked_where(cube.data == np.nan, cube.data )
                try :
                    if region is None :
                        aaa_d = np.percentile(cube.data.compressed(), 2)
                    else :
                        aaa_d = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 2)
                except :
                   print("can't get diff min-max for ", var)
                   aaa_d = cube.data.min()
                   # if region is None :
                   #     aaa_d = np.percentile(cube.data.compressed(), 2)
                   # else :
                   #     aaa_d = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 2) 
            if 'Max_diff' in attrs:
                bbb_d = attrs['Max_diff']
            else :
                cube.data = ma.masked_where(cube.data == np.nan, cube.data )
                try :
                    if region is None :
                        bbb_d = np.percentile(cube.data.compressed(), 98)
                    else :
                        bbb_d = np.percentile(cube.data[...,jj_min:jj_max,ii_min:ii_max].compressed(), 98)
                except :
                   print("can't get diff min-max for ", var)
                   bbb_d = cube.data.max()
            #
            #ccc = (-aaa_d + bbb_d)/2
            ccc = np.max([-aaa_d, bbb_d])
        ### plot :
        if SINGLE_PLOTS :
            rect = 1
        elif comp_nb < 3 :
            rect = loop + comp_nb + 2
        else :
            rect = (2 * loop) + 2
        ##

        if Meth=="Transect":
            if loop == 1 :
                if 'Min_diff' in attrs:
                    aaa_d = attrs['Min_diff']
                else :
                    aaa_d = np.percentile(cube.data[0,:,:,loci].compressed(), 2)
                if 'Max_diff' in attrs:
                    bbb_d = attrs['Max_diff']
                else :
                    bbb_d = np.percentile(cube.data[0,:,:,loci].compressed(), 98)
                #
                ccc = np.max([- aaa_d, bbb_d])
                #ccc = (-aaa_d + bbb_d)/2
            plot = subplot_proj_Orcagrid(cube[0,:,:,loci], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        elif (Meth=="Long_Avg") or (Meth=="Lat_Avg") :
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        elif (Meth=="Atl_Pac_Lat") :
            plot = subplot_proj_Orcagrid(cube[0,:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = Epipel_plot, Atl_Pac_section = True, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        elif (Meth=="Atl_Pac_hov") :
            plot = subplot_proj_Orcagrid(cube[:,0:-27], fig = fig, iii=iii, jjj=jjj, rect=rect, Proj = False, location=location, Epipel_plot = False, Atl_Pac_section = True, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        else:
            if region is None :
                if proj is None :
                    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location,cmap = cmap1, v_min= -ccc,v_max=ccc)
                else :
                    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location,cmap = cmap1, v_min= -ccc,v_max=ccc, proj=proj)
            else :
                plot = subplot_proj_regional(cube[0,jj_min:jj_max,ii_min:ii_max], fig = fig, iii=iii, jjj=jjj, rect=rect, colorbar = True, location=location,cmap = cmap1, v_min= -ccc,v_max=ccc, proj=region, longit=longit.data, latit=latit.data)
        ##
        if SINGLE_PLOTS :
            if (Meth=="Transect") or (Meth=="Long_Avg") or (Meth=="Atl_Pac_Lat"):
                plot.y_label("Depth [m]")
                plot.x_label("Latitude North [Degree]")
            elif (Meth=="Atl_Pac_hov") :
                plot.y_label("Latitude North [Degree]")
                #plot.x_label("Time (month)")
            elif  Meth=="Lat_Avg" :
                plot.y_label("Depth [m]")
                plot.x_label("Longitude East [Degree]")
            ##
            plt.savefig(f_name +"_"+ comp_ttl + " minus " + ref_ttl, bbox_inches='tight', dpi=300)
            if SHOW :
                plt.show()
        else :
            plot.subtitle(comp_ttl + " minus " + ref_ttl)
        # ++ loop
        loop = loop + 1
    ## End loop

    if not SINGLE_PLOTS :
        if Meth=="Transect":
            if LOC == "Atl" :
                subttle = "Atl. section -- "
            elif LOC == "Pac" :
                subttle = "Pacif. section -- "
            else :
                subttle = ""
        elif Meth == "2D_Extract" :
            subttle = str(depth_extract) +" m depth layer of --"
        elif Meth == "Surface" :
            subttle = "Surface values of --"
        elif Meth == "Vert_Inv" :
            subttle = "Vertically integrated --"
        elif Meth == "Long_Avg" :
            subttle = "Longitudinal average of --"
        elif Meth == "Loat_Avg" :
            subttle = "Latitudinal average of --"
        else :
            subttle = ""
        if Titre==None :
            plt.suptitle(subttle + cube.long_name, size=16)
        else :
            plt.suptitle(Titre, size=16)
        #plt.rcParams.update({'font.size': 28})
        plt.savefig(f_name)
        if SHOW :
            plt.show()










def subplot_no_diff(var, rundict, run_exp_list, Meth="None", depth_extract=100, yr="None", log_scale = False) :

    print(var)

    # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    runref = run_exp_list[0]
    ref_file_name = rundict[runlist[runref]]
    ref_ttl       = runlist[runref]

    if 'Epipel' in attrs:
        Epipel_plot = attrs['Epipel']
    elif var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL"] :
        Epipel_plot = True
    else :
        Epipel_plot = False

    plt.rcParams.update({'font.size': 12})
    cmap0 = plt.cm.get_cmap('viridis', 15)
    #cmap1 = plt.cm.get_cmap('RdBu', 15)
    cmap1 = plt.cm.get_cmap('RdBu_r', 15)
    cmap3 = plt.cm.get_cmap('RdBu', 15)
    cmap2 = plt.cm.get_cmap('viridis_r', 15)

    cmap0.set_bad(color='white')
    cmap1.set_bad(color='gray')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='gray')


    if Meth == "None" :
        f_name = "Surf_plot_" + var
    else :
        f_name = "Surf_plot_" + var + "_" + Meth
    if yr != "None" :
        f_name = f_name + "_" + yr

    comp_nb = np.size(run_exp_list) - 1
    ## prep subplot display
    if comp_nb == 1 :
        iii=1
        jjj=2
    elif comp_nb == 2 :
        iii=2
        jjj=2
    elif comp_nb == 3 :
        iii=2
        jjj=2
    elif 4 <= comp_nb <= 5 :
        iii=3
        jjj=2
    elif 6 <= comp_nb <= 8 :
        iii=3
        jjj=3
    elif 9 <= comp_nb <= 11 :
        iii=4
        jjj=3


    fig = plt.figure(figsize=(5*jjj,5*iii))
    ## no pac
    if Meth=="OneD" or Meth=="surface_OneD" :
        cube = prep_cube(ref_file_name, var, oneD=True)
    else:
        cube = prep_cube(ref_file_name, var)

    if Meth=="Surface" :
        cube = prep_cube_surf(cube)
    elif Meth=="surface_OneD" :
        cube = prep_cube_surf(cube, oneD=True)
    elif Meth=="Vert_Inv":
        cube = prep_cube_vert_inv(cube)
    elif Meth=="2D_Extract":
        try :
            cube = prep_cube_extract_depth(cube, depth_extract)
        except :
            print(" already 2D var -- ", cube.data.shape )
    #aaa = cube.data.min()
    #bbb = cube.data.max()
    if 'Min' in attrs:
        aaa = attrs['Min']
    else :
        aaa = np.percentile(cube.data.compressed(), 5)
    if 'Max' in attrs:
        bbb = attrs['Max']
    else :
        bbb = np.percentile(cube.data.compressed(), 95)
    ### plot :
    if Meth=="Transect":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[0,:,:,262].compressed(), 5)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[0,:,:,262].compressed(), 95)
        #
        plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    elif Meth=="OneD" or Meth=="surface_OneD":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :
            aaa = np.percentile(cube.data[...,0,0].compressed(), 5)
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :
            bbb = np.percentile(cube.data[...,0,0].compressed(), 95)
        #
        plot = subplot_timeseries(cube[...,0,0], fig = fig, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=outfreq, log_scale=log_scale )
        if cube.ndim == 4 :
            plot.y_label("Depth (m)")
        elif cube.ndim < 4 :
            plot.y_label(cube.units)
        plot.x_label("time in Year")
    else:
        plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = 1,colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    plot.subtitle(ref_ttl)

    loop = 1
    for compref in run_exp_list[1:] :
        comp_file_name = rundict[runlist[compref]]
        comp_ttl       = runlist[compref]

        # no_speed_pref
        if Meth == "OneD" or Meth=="surface_OneD" :
            cube = prep_cube(comp_file_name, var, oneD=True)
        else :
            cube = prep_cube(comp_file_name, var)

        if Meth=="Surface" :
            cube = prep_cube_surf(cube)
        elif Meth=="surface_OneD" :
            cube = prep_cube_surf(cube, oneD=True)
        elif Meth=="Vert_Inv":
            cube = prep_cube_vert_inv(cube)
        elif Meth=="2D_Extract":
            try :
                cube = prep_cube_extract_depth(cube, depth_extract)
            except :
                print(" already 2D var -- ", cube.data.shape )
        ### plot :
        if Meth=="Transect":
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
        elif Meth == "OneD" or Meth=="surface_OneD" :
            plot = subplot_timeseries(cube[...,0,0], fig = fig, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=outfreq, log_scale=log_scale)
            if cube.ndim == 4 :
                plot.y_label("Depth (m)")
            elif cube.ndim < 4 :
                plot.y_label(cube.units)
            plot.x_label("time in Year")
        else:
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect = (loop)+1,colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
        plot.subtitle(comp_ttl)
        # Diff :
        cube = prep_cube_diff(ref_file_name, comp_file_name, var)
        if Meth=="Surface":
            cube = prep_cube_surf(cube)
        elif Meth=="Vert_Inv":
            cube = prep_cube_vert_inv(cube)
        elif Meth=="2D_Extract":
            try :
                cube = prep_cube_extract_depth(cube, depth_extract)
            except :
                print(" already 2D var -- ", cube.data.shape )
        ### get diff min max
        #if loop == 1 :
        #    #aaa_d = cube.data.min()
        #    #bbb_d = cube.data.max()
        #    #ccc = np.max([- aaa_d, bbb_d])
        #    aaa_d = np.percentile(cube.data, 5)
        #    bbb_d = np.percentile(cube.data, 95)
        #    ccc = (-aaa_d + bbb_d)/2
        ### plot :
        #if Meth=="Transect":
        #    if loop == 1 :
        #        aaa_d = np.percentile(cube.data[0,:,:,262], 5)
        #        bbb_d = np.percentile(cube.data[0,:,:,262], 95)
        #        #ccc = np.max([- aaa_d, bbb_d])
        #        ccc = (-aaa_d + bbb_d)/2
        #    plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect=(2 * loop)+2, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        #else:
        #    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=(2 * loop)+2, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        #plot.subtitle(comp_ttl + " minus " + ref_ttl)
        # ++ loop
        loop = loop + 1
    ## End loop


    plt.suptitle(cube.long_name, size=16)
    plt.savefig(f_name)
    plt.show()







def subplot_comp_list(var, rundict_list, loc_list, filetype="ptrc", Meth="None", yr="None", log_scale = False, **attrs) :

    print(var)
        
    if 'Epipel' in attrs:
        Epipel_plot = attrs['Epipel']
    elif var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL"] :
        Epipel_plot = True
    else : 
        Epipel_plot = False
        
    plt.rcParams.update({'font.size': 12})
    cmap0 = plt.cm.get_cmap('viridis', 15) 
    #cmap1 = plt.cm.get_cmap('RdBu', 15) 
    cmap1 = plt.cm.get_cmap('RdBu_r', 15) 
    cmap3 = plt.cm.get_cmap('RdBu', 15) 
    cmap2 = plt.cm.get_cmap('viridis_r', 15)
        
    cmap0.set_bad(color='white')
    cmap1.set_bad(color='gray')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='gray')
    
    
    if Meth == "None" : 
        f_name = "Surf_plot_" + var
    else :
        f_name = "Surf_plot_" + var + "_" + Meth
    if yr != "None" : 
        f_name = f_name + "_" + yr
        
    comp_nb = np.size(loc_list) - 1
    ## prep subplot display
    if comp_nb == 0 : 
        iii=1
        jjj=1
    if comp_nb == 1 : 
        iii=1
        jjj=2
    elif comp_nb == 2 : 
        iii=2
        jjj=2
    elif comp_nb == 3 : 
        iii=2
        jjj=2    
    elif 4 <= comp_nb <= 5 : 
        iii=3
        jjj=2    
    elif 6 <= comp_nb <= 8 : 
        iii=3
        jjj=3    
    elif 9 <= comp_nb <= 11 : 
        iii=4
        jjj=3   
    elif 12 <= comp_nb <= 15 : 
        iii=4
        jjj=4       
    
    
    fig = plt.figure(figsize=(5*jjj,5*iii))
    
        # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    #runref = run_exp_list[0]
    #ref_file_name = rundict[runlist[runref]]
    #ref_ttl       = runlist[runref]
    ##
    #compref = run_exp_list[0]
    #add_file_name = adddict[runlist[runref]]
    #add_ttl       = runlist[runref]
    #
    #ax1 = fig.add_subplot(iii, jjj, 1)
    ## no pac
    
    loop = 0
    for compref in loc_list :
        print("loc_list : ",loc_list)
        rect = (loop)+1
        ax1 = fig.add_subplot(iii, jjj, rect)
        
        for runnb in rundict_list :
            #
            modclass = runclass(runnb)
            #
            print(filetype)
            if filetype=="ptrc" :
                rundict        = modclass.dict_ptrc
            elif filetype=="diad" :
                rundict        = modclass.dict_diad
            elif filetype=="grid" :
                rundict        = modclass.dict_grid  
            print(rundict)    
            comp_file_name = rundict[runlist[compref]]
            comp_ttl       = ttl_list[compref]
            #
            try: 
            # no_speed_pref
                if Meth == "OneD" or Meth=="surface_OneD" :
                    cube = prep_cube(comp_file_name, var, oneD=True)
                else :     
                    cube = prep_cube(comp_file_name, var)
                if Meth=="Surface" :
                    cube = prep_cube_surf(cube,rundict_diad=modclass.dict_diad)
                elif Meth=="surface_OneD" :
                    cube = prep_cube_surf(cube, oneD=True,rundict_diad=modclass.dict_diad) 
                elif Meth=="Vert_Inv":
                    cube = prep_cube_vert_inv(cube,rundict_diad=modclass.dict_diad)
                if 'Min' in attrs:
                    aaa = attrs['Min']
                else :     
                    aaa = np.percentile(cube.data.compressed(), 5)
                if 'Max' in attrs:
                    bbb = attrs['Max']
                else :      
                    bbb = np.percentile(cube.data.compressed(), 95)
                ### plot : 
                if Meth=="Transect":
                    if 'Min' in attrs:
                        aaa = attrs['Min']
                    else :     
                        aaa = np.percentile(cube.data[0,:,:,262].compressed(), 5) 
                    if 'Max' in attrs:
                        bbb = attrs['Max']
                    else :      
                        bbb = np.percentile(cube.data[0,:,:,262].compressed(), 95)
                    #        
                    plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
                elif Meth == "OneD" or Meth=="surface_OneD" : 
                    if 'Min' in attrs:
                        aaa = attrs['Min']
                    else :     
                        aaa = np.percentile(cube.data[...,0,0].compressed(), 5)
                        #aaa2 = np.percentile(addcube.data[...,0,0].compressed(), 5)
                        #aaa=np.min([aaa1,aaa2])
                    if 'Max' in attrs:
                        bbb = attrs['Max']
                    else :      
                        bbb = np.percentile(cube.data[...,0,0].compressed(), 95)
                        #bbb2 = np.percentile(addcube.data[...,0,0].compressed(), 95)
                        #bbb=np.min([bbb1,bbb2])
                    #    
                    plot = subplot_timeseries(cube[...,0,0], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=modclass.outfreq, log_scale=log_scale,label=modclass.name)
            
                    if cube.ndim == 4 : 
                        plot.y_label("Depth (m)")
                    elif cube.ndim < 4 :
                        plot.y_label(cube.units)
                    plot.x_label("time in Year")
                else:
                    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
                plot.subtitle(comp_ttl)
            except :
                print("Could not plot ",cube.long_name," from ", comp_file_name )
        loop = loop + 1
    ## End loop    
        
   
    plt.suptitle(cube.long_name, size=16)
    plt.savefig(f_name)
    plt.show()
    
    



def subplot_comp(var, rundict, adddict, run_exp_list, Meth="None", yr="None", log_scale = False, outfreq = "1M", **attrs) :

    print(var)
    
    # ref_file_name, interm_file_name, final_filename, ref_ttl, int_ttl, fin_ttl,
    runref = run_exp_list[0]
    ref_file_name = rundict[runlist[runref]]
    ref_ttl       = runlist[runref]
    ##
    #compref = run_exp_list[0]
    add_file_name = adddict[runlist[runref]]
    add_ttl       = runlist[runref]
    
    if 'Epipel' in attrs:
        Epipel_plot = attrs['Epipel']
    elif var in ["ZME", "ZMI", "ZMP", "PHN", "PHD", "CHL"] :
        Epipel_plot = True
    else : 
        Epipel_plot = False
        
    plt.rcParams.update({'font.size': 12})
    cmap0 = plt.cm.get_cmap('viridis', 15) 
    #cmap1 = plt.cm.get_cmap('RdBu', 15) 
    cmap1 = plt.cm.get_cmap('RdBu_r', 15) 
    cmap3 = plt.cm.get_cmap('RdBu', 15) 
    cmap2 = plt.cm.get_cmap('viridis_r', 15)
        
    cmap0.set_bad(color='white')
    cmap1.set_bad(color='gray')
    cmap2.set_bad(color='gray')
    cmap3.set_bad(color='gray')
    
    
    if Meth == "None" : 
        f_name = "Surf_plot_" + var
    else :
        f_name = "Surf_plot_" + var + "_" + Meth
    if yr != "None" : 
        f_name = f_name + "_" + yr
        
    comp_nb = np.size(run_exp_list) - 1
    ## prep subplot display
    if comp_nb == 1 : 
        iii=1
        jjj=2
    elif comp_nb == 2 : 
        iii=2
        jjj=2
    elif comp_nb == 3 : 
        iii=2
        jjj=2    
    elif 4 <= comp_nb <= 5 : 
        iii=3
        jjj=2    
    elif 6 <= comp_nb <= 8 : 
        iii=3
        jjj=3    
    elif 9 <= comp_nb <= 11 : 
        iii=4
        jjj=3    
    
    
    fig = plt.figure(figsize=(5*jjj,5*iii))
    ax1 = fig.add_subplot(iii, jjj, 1)
    ## no pac
    if Meth=="OneD" or Meth=="surface_OneD" :
        cube = prep_cube(ref_file_name, var, oneD=True)
        addcube = prep_cube(add_file_name, var, oneD=True)
    else:    
        cube = prep_cube(ref_file_name, var)
        addcube = prep_cube(add_file_name, var)
    if Meth=="Surface" :
        cube = prep_cube_surf(cube)
        addcube = prep_cube_surf(addcube,rundict_diad=adddict_diad)
    elif Meth=="surface_OneD" :
        cube = prep_cube_surf(cube, oneD=True) 
        addcube = prep_cube_surf(addcube, oneD=True,rundict_diad=adddict_diad) 
    elif Meth=="Vert_Inv":
        cube = prep_cube_vert_inv(cube)
        addcube = prep_cube_vert_inv(addcube,rundict_diad=adddict_diad)
    #aaa = cube.data.min() 
    #bbb = cube.data.max()
    if 'Min' in attrs:
        aaa = attrs['Min']
    else :     
        aaa = np.percentile(cube.data.compressed(), 5)
    if 'Max' in attrs:
        bbb = attrs['Max']
    else :      
        bbb = np.percentile(cube.data.compressed(), 95)
    ### plot : 
    if Meth=="Transect":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :     
            aaa = np.percentile(cube.data[0,:,:,262].compressed(), 5) 
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :      
            bbb = np.percentile(cube.data[0,:,:,262].compressed(), 95)
        #    
        plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    elif Meth=="OneD" or Meth=="surface_OneD":
        if 'Min' in attrs:
            aaa = attrs['Min']
        else :     
            aaa1 = np.percentile(cube.data[...,0,0].compressed(), 5)
            aaa2 = np.percentile(addcube.data[...,0,0].compressed(), 5)
            aaa=np.min([aaa1,aaa2])
        if 'Max' in attrs:
            bbb = attrs['Max']
        else :      
            bbb1 = np.percentile(cube.data[...,0,0].compressed(), 95)
            bbb2 = np.percentile(addcube.data[...,0,0].compressed(), 95)
            bbb=np.min([bbb1,bbb2])
        #    
        plot = subplot_timeseries(cube[...,0,0], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=outfreq, log_scale=log_scale,label="MEDUSA_1D" )
        plot = subplot_timeseries(addcube[...,0,0], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = 1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale,label="MEDUSA_3D")
        if cube.ndim == 4 : 
            plot.y_label("Depth (m)")
        elif cube.ndim < 4 :
            plot.y_label(cube.units)
        plot.x_label("time in Year")
    else:
        plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = 1,colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale )
    plot.subtitle(ref_ttl)
    
    loop = 1
    for compref in run_exp_list[1:] : 
        comp_file_name = rundict[runlist[compref]]
        comp_ttl       = runlist[compref]
        #
        add2_file_name = adddict[runlist[compref]]
        add2_ttl       = runlist[compref]
        ##
        rect = (loop)+1
        ax1 = fig.add_subplot(iii, jjj, rect)
        # no_speed_pref
        if Meth == "OneD" or Meth=="surface_OneD" :
            cube = prep_cube(comp_file_name, var, oneD=True)
            addcube = prep_cube(add2_file_name, var, oneD=True)
        else :     
            cube = prep_cube(comp_file_name, var)
            addcube = prep_cube(add2_file_name, var)
            
        if Meth=="Surface" :
            cube = prep_cube_surf(cube)
            addcube = prep_cube_surf(addcube,rundict_diad=adddict_diad)
        elif Meth=="surface_OneD" :
            cube = prep_cube_surf(cube, oneD=True) 
            addcube = prep_cube_surf(addcube, oneD=True,rundict_diad=adddict_diad) 
        elif Meth=="Vert_Inv":
            cube = prep_cube_vert_inv(cube)
            addcube = prep_cube_vert_inv(addcube)
        ### plot :    
        if Meth=="Transect":
            plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
        elif Meth == "OneD" or Meth=="surface_OneD" : 
            plot = subplot_timeseries(cube[...,0,0], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, outfreq=outfreq, log_scale=log_scale,label="MEDUSA_1D")
            plot = subplot_timeseries(addcube[...,0,0], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,Proj = False, Epipel_plot = Epipel_plot, colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale,label="MEDUSA_3D")

            if cube.ndim == 4 : 
                plot.y_label("Depth (m)")
            elif cube.ndim < 4 :
                plot.y_label(cube.units)
            plot.x_label("time in Year")
        else:
            plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, ax=ax1, iii=iii, jjj=jjj, rect = (loop)+1,colorbar = True, v_min= aaa,v_max=bbb, log_scale=log_scale)
        plot.subtitle(comp_ttl)
        # Diff : 
        #cube = prep_cube_diff(ref_file_name, comp_file_name, var)
        #if Meth=="Surface":
        #    cube = prep_cube_surf(cube)
        #elif Meth=="Vert_Inv":
        #    cube = prep_cube_vert_inv(cube)
        ### get diff min max    
        #if loop == 1 :    
        #    #aaa_d = cube.data.min() 
        #    #bbb_d = cube.data.max()
        #    #ccc = np.max([- aaa_d, bbb_d])
        #    aaa_d = np.percentile(cube.data, 5)
        #    bbb_d = np.percentile(cube.data, 95)
        #    ccc = (-aaa_d + bbb_d)/2
        ### plot : 
        #if Meth=="Transect":
        #    if loop == 1 :
        #        aaa_d = np.percentile(cube.data[0,:,:,262], 5) 
        #        bbb_d = np.percentile(cube.data[0,:,:,262], 95)
        #        #ccc = np.max([- aaa_d, bbb_d])
        #        ccc = (-aaa_d + bbb_d)/2
        #    plot = subplot_proj_Orcagrid(cube[0,:,:,262], fig = fig, iii=iii, jjj=jjj, rect=(2 * loop)+2, Proj = False, Epipel_plot = Epipel_plot, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        #else:
        #    plot = subplot_proj_Orcagrid(cube[0,:,:], fig = fig, iii=iii, jjj=jjj, rect=(2 * loop)+2, colorbar = True, cmap = cmap1, v_min= -ccc,v_max=ccc)
        #plot.subtitle(comp_ttl + " minus " + ref_ttl)
        # ++ loop
        loop = loop + 1
    ## End loop    
        
   
    plt.suptitle(cube.long_name, size=16)
    plt.savefig(f_name)
    plt.show()
    
    





def prep_cube_Atl_Pac(cube, runlist, run_nb=None,NEW=False, MONTHLY=None) :
    '''
    Function : prep_cube_Atl_Pac(cube)
    - get a cube
    - regrid mask atl
    - regrid mask pac
    - regrid it cube
    - make mask-Atl --> Long avg
    - make mask-Pac --> Long avg
    - stick them both nicely together
    '''
    import os

    dirr="SAVE_NC_FILE"  ## Directory with all saved regrided nc files
    ##
    ##first check dir exists and create it if needed:
    if os.path.isdir(dirr) != True :
        os.mkdir(dirr)

    if MONTHLY is not None :
        finalname = dirr + "/regridded_" + cube.var_name + "_" + runlist[run_nb] + "_ALT_PAC_m" + MONTHLY + ".nc"
        NEW=True
    #src_cubi = prep_cube(rundict_ptrc[runlist[0]],"DIC")
    ATL_MASK = "all_atl"
    PAC_MASK = "all_pac"
    ##
    ## Atl region :
    ## Read the nc file if it exists :
    fname = dirr + "/regridded_" + cube.var_name + "_" + runlist[run_nb] + "_" + ATL_MASK + "_mask"
    if MONTHLY is not None :
        fname = fname + "_2D.nc"
    else :
        fname = fname + "_3D.nc"
    print("check if file exists : ", fname)
    if os.path.isfile(fname) and not NEW :
        # file exists : read it
        print("read regridded file : ",fname)
        reg_atl = iris.load_cube(fname)
        #print(reg_atl)
    else :
        ## Mask atl :
        atl_cubi = prep_Masked_cube(cube, ATL_MASK)
        #  regrid :
        reg_atl = regrid_ORCA(atl_cubi, runlist, run_nb=run_nb, MASK=ATL_MASK, NEW=NEW)
    ##
    # Check :
    # print(reg_atl.data.shape)
    ##
    ## PAC region :
    ## Read the nc file if it exists :
    fname = dirr + "/regridded_" + cube.var_name + "_" + runlist[run_nb] + "_" + PAC_MASK + "_mask"
    if MONTHLY is not None :
        fname = fname + "_2D.nc"
    else :
        fname = fname + "_3D.nc"
        #
    print("check if file exists : ", fname)
    if os.path.isfile(fname) and not NEW :
        # file exists : read it
        print("read regridded file : ",fname)
        reg_pac = iris.load_cube(fname)
        #print(reg_pac)
    else :
        ## Mask pac :
        pac_cubi = prep_Masked_cube(cube, PAC_MASK)
        #  regrid :
        reg_pac = regrid_ORCA(pac_cubi, runlist, run_nb=run_nb, MASK=PAC_MASK, NEW=NEW)
    ##
    # Zonal average :
    long_atl = reg_atl.collapsed('longitude', iris.analysis.MEAN)
    long_pac = reg_pac.collapsed('longitude', iris.analysis.MEAN)
    ##
    #print(long_atl)
    ## Here we create an empty Cube big enough to include both Ocean mean :
    # create it
    #print("area shape     : ", area.shape)
    #xx = reg_cube.coord(axis='x')
    #yy = long_atl.coord(axis='y')
    yy1 = np.arange(-180, 180, 1)
    yy = iris.coords.DimCoord(yy1, standard_name='latitude', units='degrees')
    #print(reg_cube)
    if MONTHLY is None :
        ### 3D
        tt = long_atl.coords(axis='t')[0]
        zz = long_atl.coord(axis='z')
        ##
        ttdim = long_atl.shape[0]
        zzdim = long_atl.shape[1]
        yydim = (long_atl.shape[2] * 2 )
        #xxdim = reg_cube.shape[2]
        #print(ttdim, zzdim, yydim)
        #print(yy)
        reg_cubi = iris.cube.Cube(np.ma.zeros((ttdim, zzdim, yydim)),
                    dim_coords_and_dims=[(tt, 0),
                                         (zz, 1),
                                         (yy, 2)])
        ## Fill the cube
        ## Atl first and Pacific after
        reg_cubi.data[:,:,0:180] = np.flip(long_atl.data,2)
        reg_cubi.data[:,:,180:] = long_pac.data
    else :
        ### 2D
        ### 3D
        tt = long_atl.coords(axis='t')[0]
        #zz = long_atl.coord(axis='z')
        ##
        ttdim = long_atl.shape[0]
        #zzdim = long_atl.shape[1]
        yydim = (long_atl.shape[1] * 2 )
        #xxdim = reg_cube.shape[2]
        #print(ttdim, zzdim, yydim)
        #print(yy)
        reg_cubi = iris.cube.Cube(np.ma.zeros((ttdim, yydim)),
                    dim_coords_and_dims=[(tt, 0),
                                         #(zz, 1),
                                         (yy, 1)])
        ## Fill the cube
        ## Atl first and Pacific after
        reg_cubi.data[:,0:180] = np.flip(long_atl.data,1)
        reg_cubi.data[:,180:] = long_pac.data

    reg_cubi.name = cube.name
    reg_cubi.long_name = cube.long_name
    reg_cubi.units = cube.units
    print(reg_cubi.shape)

    if MONTHLY is not None :
        print("saving file : ",finalname)
        iris.save(reg_cubi, finalname)
    return reg_cubi




def prep_cube_Atl_Pac_hov_from_monthly(file_name, var, runlist, run_nb=None, REF_COMP=None, Obs=False, MASK_ZERO=True ,NEW_regrid=False) :
    ##
    import os
    ##
    '''
    prepare an atl-pac-averaged surface section of a whole year (month by month)
    nc file, ready to be plotted.
    the file is saved and ready to use to avoid having to redo the long regridding/averaging work.

    call like :
    cube = prep_cube_Atl_Pac_hov_from_monthly(file_name,
                                              var,
                                              run_nb=None,
                                              REF_COMP=None,
                                              Obs=False,
                                              MASK_ZERO=True ,
                                              NEW_regrid=False)
    '''
    ##
    dirr="SAVE_NC_FILE"  ## Directory with all saved regrided nc files
    ##
    ##first check dir exists and create it if needed:
    if os.path.isdir(dirr) != True :
        os.mkdir(dirr)
    ##
    finalname = dirr + "/Hovm_" + var + "_" + runlist[run_nb] + "_ALT_PAC_.nc"
    if os.path.isfile(finalname) :
        # file exists : read it
        print("read ready Hovmoell file : ",finalname)
        cube = iris.load_cube(finalname)
    else :
        ## assume monthly file name ends with *_T_m?.nc
        print(" reading monthly files for hovmoeller plot : ")
    #    print(file_name[:-4])
        for mm in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']:
            ## check for the if monthly section already exists -
            fname = dirr + "/regridded_" + var + "_" + runlist[run_nb] + "_ALT_PAC_m" + mm + ".nc"
            if os.path.isfile(fname) :
                # file exists : read it
                print("read monthly file : ",fname)
                cube1m = iris.load_cube(fname)
            else :
                # file doesn't exist - do the regrid/avg work.
                base=file_name[:-4]+mm+".nc"
                print(base)
                cube1m = prep_cube(base, var, REF_COMP=REF_COMP, Obs=Obs,MASK_ZERO=MASK_ZERO )
                ## extract surface:
                try :
                    cube1m = prep_cube_surf(cube1m)
                    if cube1m is None :
                        print("cube ref surface gives None")
                except :
                    print("already 2D -- ", cube1m.data.shape ," no need to extract ")
                ##
                ## Does mask, regrid and stick atl pand PAC
                ## might not work here because of surface only fields... let see how that goes but easy to fix
                #try :
                cube1m = prep_cube_Atl_Pac(cube1m, runlist, run_nb=run_nb, NEW=NEW_regrid, MONTHLY=mm)
                if cube1m is None :
                    print("cube ref Atl_Pac_Lat gives None")
                #except :
                #    print("Atl-Pac average not working... -- ", cube1m.data.shape )
            if mm == '1' :
                cube_list = iris.cube.CubeList([cube1m])
            else :
                cube_list.append(cube1m)

        cube = cube_list.concatenate()[0]
        iris.save(cube,finalname)
        ##
    print(cube)
    #
    return cube






#def prep_cube_transect(FDDT) :
    #temp = prep_cube(rundict_diad[runlist[0]], "FASTN")
    #print(temp)
    #FDDT = prep_cube(filename, var)
    #print(FDDT)
    #cube = temp.copy()
#    cube= FDDT.extract(iris.Constraint(longitude=-30)) ## extract mid Atlantic transect -- no constraint on
    #cube.name = FDDT.name
    #cube.long_name = FDDT.long_name
    #cube.units = FDDT.units
#
#    return cube

def prep_cube_transect(FDDT, trans) :
    if trans == "atlantic" :
        cube = FDDT.extract(iris.Constraint(latitude=330))
    elif trans == "pacific" :
        cube = FDDT.extract(iris.Constraint(latitude=168.5))
    #
    return cube

def prep_cube_extract_depth(fPOC, depth, oneD=False):
    ##
    ## B coeff is the  POC flux at 1000 to 100m ratio
    ## 1- needs the POC flux and extract both layer fields
    #
    #fPOC =  prep_cube(file_name,"DET_FLUX_C")
    #
    if depth == 1 :
        kk = 0
    elif depth == 100 :
        kk = 23
    elif depth == 200 :
        kk = 30
    elif depth == 500 :
        kk = 39
    elif depth == 1000 :
        kk = 46
    #   
    ## extract depth
    surf_slice = cube.extract(iris.Constraint(depth=kk) ) 
    if surf_slice is None :     
        surf_slice = cube.extract(iris.Constraint(coord_values={'Vertical T levels':lambda cell: cell == kk}))
    #
    return surf_slice

def prep_cube_vert_inv(FDDT, oneD=False) :
    # Sure there is a better way but : 
    cube = prep_cube_surf(FDDT, oneD=oneD)
    cube.data = FDDT.data.sum(axis=1) ## Sum on the vertical level
    print("cube shape = ", cube.data.shape)
    cube.name = FDDT.name
    cube.long_name = FDDT.long_name
    cube.units = FDDT.units

    return cube




def prep_mesh_mask(nemogrid="ORCA1_NEMO42", machine="NOC_linux") :
    '''
    prep_mesh_mask(nemogrid="ORCA1_NEMO42", machine="NOC_linux")
    returns grid adapted meshfile and regional masks
    -- machine could be :
    NOC_linux
    ARCHER2
    -- nemogrids are :
    ORCA1_NEMO42
    ORCA1_NEMO36

    '''
    if machine=="NOC_linux" :
        if nemogrid == "ORCA1_NEMO42" :
            path="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/"
            meshfile = path + "mesh_mask.nc"
            maskfile = path + "eORCA100_masks.nc"
        elif nemogrid == "ORCA1_NEMO36" :
            path="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO36/"
            meshfile = path + "mesh_mask.nc"
            maskfile = path + "eORCA100_masks.nc"
    elif machine=="ARCHER2" :
        if nemogrid == "ORCA1_NEMO42" :
            path="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/"
            meshfile = path + "mesh_mask.nc"
            maskfile = path + "eORCA100_masks.nc"
        elif nemogrid == "ORCA1_NEMO36" :
            path="/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO36/"
            meshfile = path + "mesh_mask.nc"
            maskfile = path + "eORCA100_masks.nc"

    return meshfile, maskfile




def prep_Masked_cube(cube,MASK,nemogrid="ORCA1_NEMO42", machine="NOC_linux") :
    ##

    #get correct mesh and mask :
    [meshfile, maskfile]=prep_mesh_mask(nemogrid=nemogrid, machine=machine)

    #copy cube in to avoid masking several times :
    FDDT = cube.copy()
    ##read the mask : (mask is 2D x:y only)
    cmask = prep_cube(maskfile,MASK)
    #will mask the cube now : 
    if FDDT.data.ndim == 4 :
        print("Masking 4D -- ",FDDT.data.shape," -- Mask is ", MASK)
        ttdim = FDDT.shape[0]
        zzdim = FDDT.shape[1]
        for tt in np.arange(ttdim):
            for zz in np.arange(zzdim):
                bibi = np.ma.masked_where(cmask.data !=0, FDDT[tt,zz,:,:].data)
                FDDT.data[tt,zz,:,:] = bibi
    elif FDDT.data.ndim == 3 :
        print("Masking 3D -- ",FDDT.data.shape," -- Mask is ", MASK)
        ## don't care if first dim is tt or zz : 
        ttdim = FDDT.shape[0]
        for tt in np.arange(ttdim):
            bibi = np.ma.masked_where(cmask.data != 0, FDDT[tt,:,:].data)
            FDDT.data[tt,:,:] = bibi
    elif FDDT.data.ndim == 2 :
        # already 2D
        print("Masking 2D -- ",FDDT.data.shape," -- Mask is ", MASK)
        bibi = np.ma.masked_where(cmask.data !=0, FDDT.data)
        FDDT.data = bibi
    else : 
        print("data shape not understood : ", FDDT.data.shape)
    return FDDT


def prep_cube_surf(FDDT, oneD=False,) :
    ##
    if FDDT.data.ndim == 4 :
        cube = FDDT.extract(iris.Constraint(depth=0) ) 
        if cube is None :     
            cube = FDDT.extract(iris.Constraint(coord_values={'Vertical T levels':lambda cell: cell == 0}))
        return cube
    elif FDDT.data.ndim == 3 :
        ## check the first dim is depth (more than one layer)
        ## if not, no need to extract
        if FDDT.data.shape[0] <= 12 :
            print("already 2D -- ",FDDT.data.shape," -- no need to extract --")
            return FDDT
        else :
            # real 3D
            cube = FDDT.extract(iris.Constraint(depth=0) ) 
            if cube is None :     
                cube = FDDT.extract(iris.Constraint(coord_values={'Vertical T levels':lambda cell: cell == 0}))   
            return cube
    elif FDDT.data.ndim == 2 :
        # already 2D
        print("already 2D -- ",FDDT.data.shape," -- no need to extract --")
        return FDDT
    else :
        print("data shape not understood : ", FDDT.data.shape)






def prep_cube_Obs(FDDT, diad_file) :
    temp = prep_cube(diad_file, "FASTN")
    #print(temp)
    #FDDT = prep_cube(filename, var)
    #print(FDDT)
    if FDDT.data.ndim == 4 :
        cube = temp.copy()
        cube.data = FDDT[:,0,:,:].data ## extract the surface level
        cube.name = FDDT.name
        cube.long_name = FDDT.long_name
        cube.units = FDDT.units
        return cube
    elif FDDT.data.ndim == 3 :
        ## check the first dim is depth (more than one layer)
        ## if not, no need to extract
        if FDDT.data.shape[0] == 1 :
            print("already 2D -- ",FDDT.data.shape," -- no need to extract --")
            return FDDT
        else :
            # real 3D
            cube = temp.copy()
            cube.data = FDDT[:,0,:,:].data ## extract the surface level
            cube.name = "obs"
            cube.long_name = "obs"
            cube.units = "obs units"
            return cube
    elif FDDT.data.ndim == 2 :
        # already 2D
        #print("already 2D -- ",FDDT.data.shape," -- no need to extract --")
        cube = temp.copy()
        cube.data[0,:,:] = FDDT.data ## extract the surface level
        cube.name = FDDT.name
        cube.long_name = FDDT.long_name
        cube.units = FDDT.units
        return cube
    else :
        print("data shape not understood : ", FDDT.data.shape)








def prep_cube(file_name, var, YYear=13, REF_COMP=None, Obs=False, oneD=False, MASK_ZERO=False) :
    '''
    prep_cube(file_name, var, YYear=13, REF_COMP=None, Obs=False, oneD=False, MASK_ZERO=False)
    
    Function that read a variable from a nc or mat file, and return an iris cube.
    can be called with only cube = prep_cube("filename","var"), but has also some extra options and build up variables
    the function is called with : 
    -- Filename : the file where to read the variable from
    -- var      : the variable to read. can ve a specific variable from the file, or a combination of MEDUSA's variables.
                  of which : 
                  "PP"        - total primary prod
                  "ML_PP"     - Mixed layer PP
                  "SUB_LM_PP" - sub Mixed Layer PP
                  "ZOO_growth"- with PAC included if PAC is there
                  "PHYTO"     - phytoplankton population
                  "ZOO"       - zooplankton population
                  "PLKT"      - PHYTO + ZOO 
                  "CHL"       - PHYTO + ZOO 
                  "DETFLUX_100m "
                  "DETFLUX_200m "
                  "DETFLUX_500m "
                  "DETFLUX_1000m "
                  ....
                  ....
                  
    -- YYear   : probably old unused option
    -- REF_COMP: A reference filename use to compare the variable with.
                 The resulting cube is the diff of the variable from (file_name - REF_COMP )
                 Do not fill if not needed 
    -- Obs     : Set to True if the file to read is an observation : 
                 in an automated script, like subplot_diff, where one might want to only use MEDUSA's var_name, 
                 the obs known variable are automatically used instead, and alternative file_name is looked for if needed.
    -- oneD    : to set as True if the file is One dimention (Z dim or T-Z dim only)
    -- MASK_ZERO: set to True if you want to mask the 0.0 values automatically.
    '''
    print("reading", var," in ",file_name, "REF_COMP :", REF_COMP )
    if file_name[-3:] ==".nc" :
        #print("read nc file", file_name)
        if Obs :
            var = nc_obs2D_varname(var)
            try :
                print("read Obs nc file", file_name, var)
                cubi = read_cube(file_name, var)
            except :
                print("read Obs nc file", rundict_diad[runlist[0]], var)
                cubi = read_cube(rundict_diad[runlist[0]], var)
            ## Make the cube comparable with the NEMO files
            cube = prep_cube_Obs(cubi, REF_COMP)
        else :
            if var == "PP" :
                PRN = prep_cube(file_name, "PRN", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "PRD", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data
                cube.data = cube.data * 6.625 * 12.011 * 1e-3
                cube.var_name = "PP"
                cube.long_name = "Total Primary production"
                cube.units = "g-C/m2/d"
            elif var == "ML_PP" :
                PRN = prep_cube(file_name, "ML_PRN", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "ML_PRD", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data
                cube.data = cube.data * 6.625 * 12.011 * 1e-3
                cube.var_name = "ML_PP"
                cube.long_name = "Mixed Layer Primary production"
                cube.units = "g-C/m2/d"
            elif var == "SUB_ML_PP" :
                PP    = prep_cube(file_name, "PP"   , oneD=oneD, MASK_ZERO=MASK_ZERO)
                ML_PP = prep_cube(file_name, "ML_PP", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube  = PP.copy()
                cube.data = PP.data - ML_PP.data
                #cube.data = cube.data * 6.625 * 12.011 * 1e-3
                cube.var_name = "SUB_ML_PP"
                cube.long_name = "Sub-Surface Primary production"
                cube.units = "g-C/m2/d"
            elif var == "ZOO_growth" :
                ZIG = prep_cube(file_name, "ZI_GROW", oneD=oneD, MASK_ZERO=MASK_ZERO)
                ZEG = prep_cube(file_name, "ZE_GROW", oneD=oneD, MASK_ZERO=MASK_ZERO)
                try :
                    ZPG = prep_cube(file_name, "ZP_GROW", oneD=oneD, MASK_ZERO=MASK_ZERO)
                except :
                    ZPG = ZIG.copy()
                    ZPG.data = ZPG.data * 0.0
                cube = ZIG.copy()
                cube.data = ZIG.data + ZEG.data + ZPG.data
                cube.long_name = "Total Zooplankton Growth"
            elif var == "PHYTO" :
                PRN = prep_cube(file_name, "PHN_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "PHD_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data
                cube.data = cube.data * 6.625 * 12.011 * 1e-3
                cube.long_name = "Integrated Phyto biomass"
                cube.units = "g-C/m2"
            elif var == "ZOO" :
                PRN = prep_cube(file_name, "ZMI_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "ZME_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data
                cube.data = cube.data * 5.625 * 12.011 * 1e-3
                cube.long_name = "Integrated Zoo biomass"
                cube.units = "g-C/m2"
            elif var == "PLKT" :
                PRN = prep_cube(file_name, "PHN_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "PHD_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                ZIG = prep_cube(file_name, "ZMI_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                ZEG = prep_cube(file_name, "ZME_E3T", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRN.data = PRN.data * 6.625 * 12.011 * 1e-3
                PRD.data = PRD.data * 6.625 * 12.011 * 1e-3
                ZIG.data = ZIG.data * 5.625 * 12.011 * 1e-3
                ZEG.data = ZEG.data * 5.625 * 12.011 * 1e-3
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data + ZIG.data + ZEG.data
                cube.long_name = "Integrated plankton biomass"
                cube.units = "g-C/m2"
            elif var == "CHL" :
                PRN = prep_cube(file_name, "CHN", oneD=oneD, MASK_ZERO=MASK_ZERO)
                PRD = prep_cube(file_name, "CHD", oneD=oneD, MASK_ZERO=MASK_ZERO)
                cube = PRN.copy()
                cube.data = PRN.data + PRD.data
                cube.var_name = "CHL" 
                cube.long_name = "Total Chlorophyll"
                cube.units = "mg-C/m3"
            elif var == "DETFLUX_100m":
                ## 2D cube template.
                # FDT_100; FDT_200; FDT_200; FDT_1000;  (fast N only)
                # SDC_100; SDC_200; SDC_200; SDC_1000;  (slow C only)
                # DETFLUX3 (tot N flux)
                # FDS_NIT3; FD1_NIT3; FDS_CAR3; FD1_CAR3 ; 3D N or C slow or fast
                ### ==> DETFLUX3
                DETF3 = prep_cube(file_name, "DETFLUX3", oneD=oneD, MASK_ZERO=MASK_ZERO) ## 3D
                ### Extract 2D from specific layer :
                cube = prep_cube_extract_depth(DETF3, 100)
                cube.long_name = "sinking N flux at 100m depth"
            elif var == "DETFLUX_200m":
                ### ==> DETFLUX3
                DETF3 = prep_cube(file_name, "DETFLUX3", oneD=oneD, MASK_ZERO=MASK_ZERO) ## 3D
                ### Extract 2D from specific layer :
                cube = prep_cube_extract_depth(DETF3, 200)
                cube.long_name = "sinking N flux at 200m depth"
            elif var == "DETFLUX_500m":
                ### ==> DETFLUX3
                DETF3 = prep_cube(file_name, "DETFLUX3", oneD=oneD, MASK_ZERO=MASK_ZERO) ## 3D
                ### Extract 2D from specific layer :
                cube = prep_cube_extract_depth(DETF3, 500)
                cube.long_name = "sinking N flux at 500m depth"
            elif var == "DETFLUX_1000m":
                ### ==> DETFLUX3
                DETF3 = prep_cube(file_name, "DETFLUX3", oneD=oneD, MASK_ZERO=MASK_ZERO) ## 3D
                ### Extract 2D from specific layer :
                cube = prep_cube_extract_depth(DETF3, 1000)
                cube.long_name = "sinking N flux at 1000m depth"
            elif var == "TRANSF_EFF" :
                ## flux at 100m :
                DET100 = prep_cube(file_name, "DETFLUX_100m", oneD=oneD, MASK_ZERO=MASK_ZERO)
                ## flux at 1000m :
                DET1000 = prep_cube(file_name, "DETFLUX_1000m", oneD=oneD)
                ## Transfer eff = flux 1000 / flux 100
                cube = DET100.copy()
                cube.data = DET1000.data / DET100.data * 100.0
                cube.long_name = "Transfer Efficiency"
                cube.units = "%"
            elif var == "TOT_C" :
            #     z3d(:,:,:) = (trn(:,:,:,jpmphn) * xthetapn)  +  &
            #                 (trn(:,:,:,jpmphd) * xthetapd)  +  &
            #                 (trn(:,:,:,jpmzmi) * xthetazmi) +  &
            #                 (trn(:,:,:,jpmzme) * xthetazme) +  &
            #                 trn(:,:,:,jpmdtc) + trn(:,:,:,jpmdic)
            #ENDIF
            #z2d(:,:)   = zn_sed_c(:,:) + zn_sed_ca(:,:)
            #zsum3d     = glob_sum( 'trcrrst', z3d(:,:,:) * zvol(:,:,:) )
            #zsum2d     = glob_sum( 'trcrrst', z2d(:,:) * zarea(:,:) )
            #! total tracer, and delta
            #zinvt      = zsum3d + zsum2d
            #### for each of the pool containinh C :
                CC1  = prep_cube(file_name, "PHN_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC1a = prep_cube_vert_inv(CC1)
                #
                CC2  = prep_cube(file_name, "PHD_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC2a = prep_cube_vert_inv(CC2)
                #
                CC3  = prep_cube(file_name, "ZMI_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC3a = prep_cube_vert_inv(CC3)
                #
                CC4  = prep_cube(file_name, "ZME_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC4a = prep_cube_vert_inv(CC4)
                #
                CC5  = prep_cube(file_name, "DTC_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC5a = prep_cube_vert_inv(CC5)
                #
                CC6  = prep_cube(file_name, "DIC_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC6a = prep_cube_vert_inv(CC6)
                ## dont forget benthic C on diad file :
                file_d      = diad_replace(file_name)
                if REF_COMP is None :
                    file_d_comp = None
                else :
                    file_d_comp = diad_replace(REF_COMP)
                CC7  = prep_cube(file_d, "BEN_C" , REF_COMP=file_d_comp, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC8  = prep_cube(file_d, "BEN_CA", REF_COMP=file_d_comp, oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                cube = CC7.copy()
                cube.data = 6.625 * ( CC1a.data + CC2a.data + CC3a.data + CC4a.data) + CC5a.data + CC6a.data + CC7.data + CC8.data
                cube.long_name = "Total carbon"
                cube.units = "mmol-C/m2"
            elif var == "TOT_C_org" :
                CC1  = prep_cube(file_name, "PHN_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC1a = prep_cube_vert_inv(CC1)
                #
                CC2  = prep_cube(file_name, "PHD_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC2a = prep_cube_vert_inv(CC2)
                #
                CC3  = prep_cube(file_name, "ZMI_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC3a = prep_cube_vert_inv(CC3)
                #
                CC4  = prep_cube(file_name, "ZME_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC4a = prep_cube_vert_inv(CC4)
                #
                CC5  = prep_cube(file_name, "DTC_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC5a = prep_cube_vert_inv(CC5)
                #
                cube = CC1a.copy()
                cube.data = 6.625 * ( CC1a.data + CC2a.data + CC3a.data + CC4a.data) + CC5a.data
                cube.long_name = "Total org carbon"
                cube.units = "mmol-C/m2"
            elif var == "TOT_C_ben" :
                ## dont forget benthic C on diad file :
                #file_d = file_name[:-9]+"diad"+file_name[-5:]
                CC7  = prep_cube(file_name, "BEN_C", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC8  = prep_cube(file_name, "BEN_CA", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                cube = CC7.copy()
                cube.data = CC7.data + CC8.data
                cube.long_name = "Total benthic carbon"
                cube.units = "mmol-C/m2"
            elif var == "TOT_A" :
            #         z3d(:,:,:) = trn(:,:,:,jpmalk)
            #  z2d(:,:)   = zn_sed_ca(:,:) * 2.0
            #  zsum3d     = glob_sum( 'trcrrst', z3d(:,:,:) * zvol(:,:,:) )
            #  zsum2d     = glob_sum( 'trcrrst', z2d(:,:) * zarea(:,:) )
            #  ! total tracer, and delta
            #  zinvt      = zsum3d + zsum2d
                CC1  = prep_cube(file_name, "ALK_E3T", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC1a = prep_cube_vert_inv(CC1)
                #
                file_d      = diad_replace(file_name) ## same replace
                if REF_COMP is None :
                    file_d_comp = None
                else :
                    file_d_comp = diad_replace(REF_COMP) ## same replace
                CC7  = prep_cube(file_d, "BEN_CA", REF_COMP=file_d_comp, oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                cube = CC7.copy()
                cube.data = CC1a.data + 2 * CC7.data
                cube.long_name = "Total Alkalinity"
                cube.units = "mmol/m2"
            elif var == "TOT_A_ben" :
                CC7  = prep_cube(file_name, "BEN_CA", REF_COMP=REF_COMP, oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                cube = CC7.copy()
                cube.data = 2 * CC7.data
                cube.long_name = "Total Alkalinity"
                cube.units = "mmol/m2"
            elif var == "RATIO_C_A" :
                ## Tot
                CC1a  = prep_cube(file_name, "TOT_C", oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC1b  = prep_cube(REF_COMP , "TOT_C", oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                CC2a  = prep_cube(file_name, "TOT_A", oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC2b  = prep_cube(REF_COMP , "TOT_A", oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                ## diff TOT_C :
                CC1 = CC1a.copy()
                CC1.data = CC1a.data - CC1b.data
                ## diff TOT_A :
                CC2 = CC2a.copy()
                CC2.data = CC2a.data - CC2b.data
  
                cube = CC2.copy()
                try:
                    cube.data = CC1.data / CC2.data
                except :
                    cube.data = 0.0

                cube.long_name = "Ratio of delta Tot C on delta Tot Alk"
                cube.units = "  "
            elif var == "DPCO2" :
            #         z3d(:,:,:) = trn(:,:,:,jpmalk)
            #  z2d(:,:)   = zn_sed_ca(:,:) * 2.0
            #  zsum3d     = glob_sum( 'trcrrst', z3d(:,:,:) * zvol(:,:,:) )
            #  zsum2d     = glob_sum( 'trcrrst', z2d(:,:) * zarea(:,:) )
            #  ! total tracer, and delta
            #  zinvt      = zsum3d + zsum2d
                CC1  = prep_cube(file_name, "ATM_PCO2", oneD=oneD, MASK_ZERO=MASK_ZERO)
                CC2  = prep_cube(file_name, "OCN_PCO2", oneD=oneD, MASK_ZERO=MASK_ZERO)
                #
                cube = CC2.copy()
                cube.data = CC2.data - CC1.data
                cube.long_name = "Ocn-Atm PCO2 difference"
                cube.units = "uatm"
            elif (var == "nav_lon") or (var == "longitude") :
                daata = read_cube("/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc", "glamt")
                temp  = prep_cube(rundict_diad[runlist[1]], "FASTN")
                cube = temp.copy()
                cube.data[0,...] = daata.data[0,...]
                print("cube shape = ", cube.data.shape)
                cube.name = daata.name
                cube.long_name = "longitude"
                cube.units = "degree"
            elif (var == "nav_lat") or (var == "latitude") :
                daata = read_cube("/noc/users/jpp1m13/WORKING/UKESM/MESH/eORCA1/NEMO42/mesh_mask.nc", "gphit")
                temp  = prep_cube(rundict_diad[runlist[1]], "FASTN")
                cube = temp.copy()
                cube.data[0,...] = daata.data[0,...]
                print("cube shape = ", cube.data.shape)
                cube.name = daata.name
                cube.long_name = "latitude"
                cube.units = "degree"
            else :
                #print(" REF_COMP :", REF_COMP)
                if REF_COMP is None :
                    #print(" REF_COMP is none loop")
                    cube = read_cube(file_name, var, MASK_ZERO=MASK_ZERO)
                    ## fill the halo :
                    if oneD == False :
                        cube.data[...,-1] = cube.data[...,1]
                        cube.data[...,0] = cube.data[...,1]
                        #bibi.data[...,0] = bibi.data[...,1]
                        cube.data.mask[...,-1] = cube.data.mask[...,1]
                        cube.data.mask[...,0] = cube.data.mask[...,1]
                        ### mask 0.0
                        #cube.data = ma.masked_where(cube.data == 0.0, cube.data )
                        #print('got masked ', cube)
                        #print('extracted: ', cube)
                        #plt.rcParams.update({'font.size': 18})
                        #plot_proj_Orcagrid(cube[0,:,:])
                else :
                    print(" REF_COMP filled loop")
                    cube = prep_cube_diff(REF_COMP, file_name, var, oneD=oneD, MASK_ZERO=MASK_ZERO)
    #else :
    elif file_name[-3:] == "mat" :
        print("read mat file", file_name)
        f = h5py.File(file_name,'r')
        f.keys()
        #print("--", obsvar)
        try:
            obsvar = mat_obs2D_varname(var)
            data = f.get(obsvar)
            obs2D = True
            data = np.array(data) # For converting to a NumPy array
            print(data.shape)
            ndim = data.ndim
            newdata = np.swapaxes(data,ndim-2,ndim-1)
            nshape  = newdata.shape
        except :
            obsvar = mat_obs3D_varname(var)
            data = f.get(obsvar)
            obs2D = False
            data = np.array(data) # For converting to a NumPy array
            print(data.shape)
            ndim = data.ndim
            newdata = np.swapaxes(data,ndim-2,ndim-1)
            nshape  = newdata.shape
        #
        #print("reading obs ", obsvar )
        #data = np.array(data) # For converting to a NumPy array
        #print(data.shape)
        #ndim = data.ndim
        #newdata = np.swapaxes(data,ndim-2,ndim-1)
        #nshape  = newdata.shape
        print(newdata.ndim)
        #
        print(rundict_diad[runlist[1]], var)
        try :
            temp = prep_cube(rundict_diad[runlist[1]], var, oneD=oneD)
        except :
            temp = prep_cube(rundict_ptrc[runlist[1]], var, oneD=oneD)
        print("temp : ", temp.data.shape)
        if obs2D :
            try :
                bibi = prep_cube_surf(temp)
            except :
                print(temp.data.shape)
                bibi = temp
            print("bibi : ",bibi.data.shape)
            cube = bibi.copy()
            #print(cube)
            if ndim == 2 :
                print("getting 2D obs data in cube")
                #
                #plt.contourf(newdata[:,:])
                #plt.show()
                try :
                    cube.data[0,:,:] = newdata[:,:]
                except :
                    cube.data[:,:] = newdata[:,:]
                #
            elif ndim == 3 :
                print("getting 3D obs data in cube")
                # need to select a month --
                # default : annual mean
                #
                #plt.contourf(newdata[YYear-1,:,:])
                #plt.show()
                try :
                    cube.data[0,:,:] = newdata[YYear-1,:,:]
                except :
                    cube.data[:,:] = newdata[YYear-1,:,:]
                #
            elif ndim == 4 :
                print("getting 4D obs data in cube")
                ## NPP - also select the npp algo --
                # default is the mean of the 3 algo
                # need to select a month --
                # default : annual mean
                #
                #plt.contourf(newdata[3,YYear-1,:,:])
                #plt.show()
                try :
                    cube.data[0,:,:] = newdata[3,YYear-1,:,:]
                except :
                    cube.data[:,:] = newdata[3,YYear-1,:,:]
                #
        else :
            # 3D obs to fill in cube
            # obs shape: 13,75,332,362
            #print(temp.data.shape)
            bibi = temp
            print("3D bibi : ",bibi.data.shape)
            cube = bibi.copy()
            #print(cube)
            if ndim == 3 :
                print("getting 3D obs data in cube")
                try :
                    cube.data[0,:,:,:] = newdata[:,:,:]
                except :
                    cube.data[:,:,:] = newdata[:,:,:]
                #
            elif ndim == 4 :
                print("getting 4D obs data in cube")
                # need to select a month --
                # default : annual mean
                #
                #plt.contourf(newdata[3,YYear-1,:,:])
                #plt.show()
                try :
                    cube.data[0,:,:,:] = newdata[YYear-1,:,:,:]
                except :
                    cube.data[:,:,:] = newdata[YYear-1,:,:,:]
        #
        cube.data.mask = bibi.data.mask
        ### manage if needed :
        #if var == "OCN_PCO2" :
        #    cube.data = cube.data - 380.0
        if var == "CO2FLUX" :
            cube.data = cube.data * -1 * 1e15 * 1e-3 / 12.011 / 360

        #
    #cube.data = np.ma.masked_invalid(cube.data )
    #
    ###############
    if REF_COMP is not None and not Obs :
        cube.long_name = "delta " + cube.long_name
    ## Correct units if needed -
    #if var == "SILic" : 
    #    cube.units = "mmol-Si/m3"
    #if var == "ALKalin" :
    #    cube.units = "meq/m3"

    return cube


def prep_cube_diff(file_name_ref, file_name, var, oneD=False) :
    #print(var)
    print("prep_cube_diff : ", file_name, "ref file ", file_name_ref, "var : ",var)
    cube_ref = prep_cube(file_name_ref, var, oneD=oneD)
    #
    cube = prep_cube(file_name, var, oneD=oneD)
    #
    cube_diff = cube.copy()
    cube_diff.data = cube.data - cube_ref.data

    return cube_diff






def read_cube(nc_P, var, mask = True, hard_mask = True, MASK_ZERO=False) :
    v_TT = iris.Constraint(cube_func=(lambda c: c.var_name == var))
    TT = iris.load(nc_P, constraints=v_TT)[0]
    if mask == True :
        SS = np.ma.array(TT.data, mask= np.absolute(TT.data) >1e17, hard_mask=hard_mask)
        if MASK_ZERO : 
            SS = np.ma.masked_where(TT.data==0, SS)
        #
        cube = TT.copy(SS)
    else :
        cube = TT.copy()
    return cube





def plot_timeseries_from_file(var, runlist, run_exp_list, region_list = [None, "NorthOrtho"], YY=None, SHOW=False) :
    '''
    YY is number of year to plot 
    '''
    from os.path import exists
    ##
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1, 1, 1)
    ##
    reg_names="reg"
    get_name = True
    for exp in run_exp_list :
        ###
        ##1 - read saved time_serie if exist
        for region in region_list :
            if region is None :
                tm_name=var+"_"+runlist[exp]+"_time_serie"
                label = runlist[exp]
                if get_name :
                    reg_names = reg_names + "_global"
            elif region == "NorthOrtho" :
                tm_name=var+"_"+runlist[exp]+"_"+region+"_time_serie"
                label = runlist[exp] + " Arctic"
                if get_name :
                    reg_names = reg_names + "_Arctic"
            else :
                tm_name=var+"_"+runlist[exp]+"_"+region+"_time_serie"
                label = runlist[exp] + " " + region
                if get_name :
                    reg_names = reg_names + "_" + region
            print(tm_name)
            get_name = False
            if exists(tm_name+".npy") :
                tm = np.load(tm_name+".npy")
                print(tm.shape)
                tm[tm>1e10]=np.nan
                if YY is None : 
                    bb = plt.plot(tm[0,:],tm[1,:], axes = ax, label = label)
                else : 
                    bb = plt.plot(tm[0,:YY],tm[1,:YY], axes = ax, label = label)
    ##
    ax.set_xlabel("Years")
    ax.set_ylabel(var + " (in mmol m-3)")
    ax.legend()
    #
    f_name=var + "_river_timeseries" + reg_names
    plt.savefig(f_name)
    if SHOW :
        plt.show()






class subplot_proj_regional(object):
    """
    Plot projected maps
    """

    def __init__(self, cube,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, Proj = True,
                 Epipel_plot = False, cmap = None, norm = None,
                 BIOMASK = None, ADDCTOUR = None,
                 colorbar = False, location=None, longit = None, latit=None,
                 xrev = False, ytick_loc = "left", **attrs):

        """
        Tri-Polar Grid Projected Plotting
        =================================

        This example demonstrates cell plots of data on the semi-structured ORCA2 model
        grid.

        First, the data is projected into the PlateCarree coordinate reference system.

        Second four pcolormesh plots are created from this projected dataset,
        using different projections for the output image.
        --Projection list:
        Mollweide
        PlateCarree
        NorthPolarStereo
        SouthPolarStereo
        Orthographic
        SouthOrtho
        NorthOrtho
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as pltc
        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA
        from matplotlib import ticker
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # Load data
        #filepath = iris.sample_data_path('orca2_votemper.nc')
        #cube = iris.load_cube(filepath)

        # Choose plot projections
        if Proj == True :
            projections = {}
            projections['Mollweide'] = ccrs.Mollweide()
            projections['Robinson'] = ccrs.Robinson(central_longitude=0, globe=None)
            projections['Mollweide_pac'] = ccrs.Mollweide(central_longitude=-170)
            projections['Mollweide_ind'] = ccrs.Mollweide(central_longitude= 90)
            projections['PlateCarree'] = ccrs.PlateCarree()
            projections['NorthPolarStereo'] = ccrs.NorthPolarStereo()
            projections['SouthPolarStereo'] = ccrs.SouthPolarStereo()
            projections['Orthographic'] = ccrs.Orthographic(central_longitude=-90,
                                                        central_latitude=45)
            projections['SouthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                        central_latitude=-90)
            projections['NorthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                        central_latitude=+90)
            projections['OrthoAtl'] = ccrs.Orthographic(central_longitude=-30,
                                                        central_latitude=15)
            projections['OrthoPac'] = ccrs.Orthographic(central_longitude=160,
                                                        central_latitude=7)
            projections['OrthoInd'] = ccrs.Orthographic(central_longitude=80,
                                                        central_latitude=-20)
            projections['Amazon']  = ccrs.PlateCarree(central_longitude=-50)
            projections['Bengal']  = ccrs.PlateCarree(central_longitude=90)
            projections['UKshelv'] = ccrs.PlateCarree(central_longitude=-2)
            pcarree = projections['PlateCarree']

            if 'proj' in attrs :
                proj = attrs['proj']
            else:
                name='Robinson'


            if proj is None :
                name='Robinson'
            elif proj == 'Amazon' or proj == 'Bengal' or proj == 'UKshelv' :
                name = proj
            else :
                name=proj

        ctour = attrs['contour'] if 'contour' in attrs else False
        log_sc = attrs['log_scale'] if 'log_scale' in attrs else False

        aaa = cube.data.min()
        bbb = cube.data.max()
        #ttmp = prep_cube(rundict_ptrc[runlist[0]],"SIL")
        if latit is None :
            print("latit from cube")
            latit = cube.coord('latitude').points
        if longit is None :
            print("longit from cube")
            longit = cube.coord('longitude').points
        #print(longit)
        #print(longit.mask )
        #print("################")
        #print(latit)
        #print(latit.mask)
        #print("################")
        #longit(cube.data.mask) = np.nan
        #latit(cube.data.mask) = np.nan

        print('min var = ',aaa,'; max var = ',bbb)
        ###

        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        if ax is None :
            if Proj == True :
                ax = fig.add_subplot(iii, jjj, rect, projection=projections[name])
                #fig.add_subplot(ax)
                #ax = plt.subplot(projection=projections[name])
                # Set limits
            else :
                ax = fig.add_subplot(iii, jjj, rect)
            if Proj == True :
                boulet=1
                #ax.set_global()

        if cmap :
            mycmap = cmap
        elif ('compar' in attrs) :
            mycmap = plt.cm.get_cmap("BrBG")# jet, seismic, spectral
        else :
            mycmap = plt.cm.get_cmap("jet")# jet, seismic, spectral

        # plot with Iris quickplot pcolormesh

        loc_in_attr = location is not None

        if colorbar and not loc_in_attr :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                elif log_sc :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax,transform=ccrs.PlateCarree())
                else :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            cube.data = np.log10(cube.data)
                            log_lev = np.log10(levels)
                            bb = plt.contour(longit,latit,cube.data, log_lev, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                        else:
                            bb = plt.contour(longit,latit,cube.data, levels, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                        ax.clabel(bb)
                    else :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                else :
                    bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())


        else :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                elif log_sc :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax,transform=ccrs.PlateCarree())
                else :
                    if ctour :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())
                    else :
                        bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax,transform=ccrs.PlateCarree())

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            cube.data = np.log10(cube.data)
                            log_lev = np.log10(levels)
                            bb = plt.contour(longit,latit,cube.data, log_lev, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                        else:
                            bb = plt.contour(longit,latit,cube.data, levels, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                        ax.clabel(bb)
                    else :
                        bb = plt.contourf(longit,latit,cube.data, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                else :
                    bb = plt.pcolormesh(longit,latit,cube.data, cmap=mycmap, axes = ax,transform=ccrs.PlateCarree())
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())


        #if location is not None :
                #location = attrs['location']
        cbar = plt.colorbar(bb, ax=ax, location='bottom', shrink=0.7, aspect=10)

                #divider = make_axes_locatable(ax)
                ## creating new axes on the right
                ## side of current axes(ax).
                ## The width of cax will be 5% of ax
                ## and the padding between cax and ax
                ## will be fixed at 0.05 inch.
                #colorbar_axes = divider.append_axes(location,
                #                                    size="10%",
                #                                    pad=0.1)
                ## Using new axes for colorbar
                #cbar = plt.colorbar(bb, cax=colorbar_axes)

                ## colorbar label :
        cbar.set_label(cube.units)



        # Draw coastlines
        try:
            ax.coastlines()
        except:
            print("no coastline to draw")
        if Epipel_plot == True :
            ax.set_ylim([500, 0.0])
        self.ax = ax         # Graphical axes
        if xrev :
            ax.yaxis.tick_right()
        if ytick_loc == "right" :
            ax2.invert_xaxis()


    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_ylabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)




class subplot_timeseries(object):
    """
    Plot projected maps
    """

    def __init__(self, cube,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, Proj = False,
                 Epipel_plot = False, cmap = None, norm = None,
                 BIOMASK = None, ADDCTOUR = None,
                 colorbar = False, location=None, longit = None, latit=None, outfreq="1M",
                 xrev = False, ytick_loc = "left", label = None, **attrs):

        """
        Tri-Polar Grid Projected Plotting
        =================================

        This example demonstrates cell plots of data on the semi-structured ORCA2 model
        grid.

        First, the data is projected into the PlateCarree coordinate reference system.

        Second four pcolormesh plots are created from this projected dataset,
        using different projections for the output image.
        --Projection list:
        Mollweide
        PlateCarree
        NorthPolarStereo
        SouthPolarStereo
        Orthographic
        SouthOrtho
        NorthOrtho
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as pltc
        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA
        from matplotlib import ticker
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # Load data
        #filepath = iris.sample_data_path('orca2_votemper.nc')
        #cube = iris.load_cube(filepath)



        ctour = attrs['contour'] if 'contour' in attrs else False
        log_sc = attrs['log_scale'] if 'log_scale' in attrs else False

        aaa = cube.data.min()
        bbb = cube.data.max()
        #ttmp = prep_cube(rundict_ptrc[runlist[0]],"SIL")

        print('min var = ',aaa,'; max var = ',bbb)
        ###

        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            fig = plt.figure(figsize=(20,10))

        if ax is None :
            ax = fig.add_subplot(iii, jjj, rect)

        if cmap :
            mycmap = cmap
        elif ('compar' in attrs) :
            mycmap = plt.cm.get_cmap("BrBG")# jet, seismic, spectral
        else :
            mycmap = plt.cm.get_cmap("jet")# jet, seismic, spectral

        # plot with Iris quickplot pcolormesh

        loc_in_attr = location is not None
        #print("ndim = ", cube.ndim)
        time_long = cube.shape[0]
        day = np.arange(time_long)

        ##
        ## Want time axe in year.
        if outfreq=="1M" :
            longit = day /12
        elif outfreq == "1D" :
            longit = day /360
        elif outfreq == "5D" :
            longit = day * 5 / 360
        elif outfreq == "10D" :
            longit = day * 10 / 360

        if cube.ndim > 1 :
            deep = cube.coord('Vertical T levels').points
            latit  = deep

        cubbi = np.transpose(cube.data)

        if cube.ndim == 1 :
            bb = plt.plot(longit,cubbi, axes = ax, label = label)

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                ax.set_ylim([v_min, v_max])
            if log_sc :
                plt.semilogy()

        elif colorbar and not loc_in_attr :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                elif log_sc :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                else :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            cube.data = np.log10(cube.data)
                            log_lev = np.log10(levels)
                            bb = plt.contour(longit,latit,cubbi, log_lev, cmap=mycmap, axes = ax)
                        else:
                            bb = plt.contour(longit,latit,cubbi, levels, cmap=mycmap, axes = ax)
                        ax.clabel(bb)
                    else :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, axes = ax)
                else :
                    bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, axes = ax)
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())


        else :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                elif log_sc :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                else :
                    if ctour :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            cube.data = np.log10(cube.data)
                            log_lev = np.log10(levels)
                            bb = plt.contour(longit,latit,cubbi, log_lev, cmap=mycmap, axes = ax)
                        else:
                            bb = plt.contour(longit,latit,cubbi, levels, cmap=mycmap, axes = ax)
                        ax.clabel(bb)
                    else :
                        bb = plt.contourf(longit,latit,cubbi, cmap=mycmap, axes = ax)
                else :
                    bb = plt.pcolormesh(longit,latit,cubbi, cmap=mycmap, axes = ax)
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())


        #if location is not None :
                #location = attrs['location']
        if cube.ndim > 1 :
            cbar = plt.colorbar(bb, ax=ax, location='bottom', shrink=0.7, aspect=10)

                #divider = make_axes_locatable(ax)
                ## creating new axes on the right
                ## side of current axes(ax).
                ## The width of cax will be 5% of ax
                ## and the padding between cax and ax
                ## will be fixed at 0.05 inch.
                #colorbar_axes = divider.append_axes(location,
                #                                    size="10%",
                #                                    pad=0.1)
                ## Using new axes for colorbar
                #cbar = plt.colorbar(bb, cax=colorbar_axes)

                ## colorbar label :
            cbar.set_label(cube.units)
            # Draw coastlines
            try:
                ax.coastlines()
            except:
                print("no coastline to draw")
            if Epipel_plot == True :
                ax.set_ylim([500, 0.0])
            else :
                ax.set_ylim([1000, 0.0])
        else :
            ax.legend()

        self.ax = ax         # Graphical axes
        if xrev :
            ax.yaxis.tick_right()
        if ytick_loc == "right" :
            ax2.invert_xaxis()


    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_ylabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)







class subplot_proj_Orcagrid(object):
    """
    Plot projected maps
    """

    def __init__(self, cube,
                 fig=None, ax=None, iii=1, jjj=1, rect=1, Proj = True,
                 Epipel_plot = False, cmap = None, norm = None,
                 BIOMASK = None, ADDCTOUR = None, COAST_LINE=True, 
                 colorbar = False, location=None,
                 Atl_Pac_section = False,
                 xrev = False, ytick_loc = "left", **attrs):

        """
        Tri-Polar Grid Projected Plotting
        =================================

        This example demonstrates cell plots of data on the semi-structured ORCA2 model
        grid.

        First, the data is projected into the PlateCarree coordinate reference system.

        Second four pcolormesh plots are created from this projected dataset,
        using different projections for the output image.
        --Projection list:
        Mollweide
        PlateCarree
        NorthPolarStereo
        SouthPolarStereo
        Orthographic
        SouthOrtho
        NorthOrtho
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as pltc
        import cartopy.crs as ccrs
        import iris
        import iris.analysis.cartography
        import iris.plot as iplt
        import iris.quickplot as qplt
        import mpl_toolkits.axisartist.floating_axes as FA
        from matplotlib import ticker
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # Load data
        #filepath = iris.sample_data_path('orca2_votemper.nc')
        #cube = iris.load_cube(filepath)

        # Choose plot projections
        if Proj == True :
            projections = {}
            projections['Mollweide'] = ccrs.Mollweide()
            projections['Robinson'] = ccrs.Robinson(central_longitude=0, globe=None)
            projections['Mollweide_pac'] = ccrs.Mollweide(central_longitude=-170)
            projections['Mollweide_ind'] = ccrs.Mollweide(central_longitude= 90)
            projections['PlateCarree'] = ccrs.PlateCarree()
            projections['NorthPolarStereo'] = ccrs.NorthPolarStereo()
            projections['SouthPolarStereo'] = ccrs.SouthPolarStereo()
            projections['Orthographic'] = ccrs.Orthographic(central_longitude=-90,
                                                        central_latitude=45)
            projections['SouthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                        central_latitude=-90)
            projections['NorthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                        central_latitude=+90)
            projections['OrthoAtl'] = ccrs.Orthographic(central_longitude=-30,
                                                        central_latitude=15)
            projections['OrthoPac'] = ccrs.Orthographic(central_longitude=160,
                                                        central_latitude=7)
            projections['OrthoInd'] = ccrs.Orthographic(central_longitude=80,
                                                        central_latitude=-20)
            projections['Amazon']  = ccrs.Orthographic(central_longitude=-50,
                                                        central_latitude=2)
            projections['Bengal']  = ccrs.Orthographic(central_longitude=90,
                                                        central_latitude=16)
            projections['UKshelv'] = ccrs.Orthographic(central_longitude=-2,
                                                        central_latitude=55)
            pcarree = projections['PlateCarree']

            if 'proj' in attrs :
                proj = attrs['proj']
            else:
                proj = None

            if proj is None :
                name='Mollweide'
                # Transform cube to target projection
                new_cube, extent = iris.analysis.cartography.project(cube, pcarree,
                                                             nx=400, ny=200)
                if BIOMASK :
                    biocubes = iris.load(BIOMASK)
                    biocube_ref = biocubes[0].copy()
                    cube_ref, extent = iris.analysis.cartography.project(biocube_ref, pcarree,
                                                             nx=400, ny=200)
                if ADDCTOUR :
                    ctour_ref, extent = iris.analysis.cartography.project(ADDCTOUR, pcarree,
                                                             nx=400, ny=200)
            elif proj == 'Amazon' or proj == 'Bengal' or proj == 'UKshelv' :
                name = proj
                # Transform cube to target projection
                new_cube, extent = iris.analysis.cartography.project(cube, pcarree,
                                                             nx=50, ny=30)
                if BIOMASK :
                    biocubes = iris.load(BIOMASK)
                    biocube_ref = biocubes[0].copy()
                    cube_ref, extent = iris.analysis.cartography.project(biocube_ref, pcarree,
                                                             nx=50, ny=30)
                if ADDCTOUR :
                    ctour_ref, extent = iris.analysis.cartography.project(ADDCTOUR, pcarree,
                                                             nx=50, ny=30)
            else :
                name=proj
                # Transform cube to target projection
                new_cube, extent = iris.analysis.cartography.project(cube, pcarree,
                                                             nx=400, ny=200)
                if BIOMASK :
                    biocubes = iris.load(BIOMASK)
                    biocube_ref = biocubes[0].copy()
                    cube_ref, extent = iris.analysis.cartography.project(biocube_ref, pcarree,
                                                             nx=400, ny=200)
                if ADDCTOUR :
                    ctour_ref, extent = iris.analysis.cartography.project(ADDCTOUR, pcarree,
                                                             nx=400, ny=200)

        ctour = attrs['contour'] if 'contour' in attrs else False
        log_sc = attrs['log_scale'] if 'log_scale' in attrs else False

        aaa = cube.data.min()
        bbb = cube.data.max()

        print('min var = ',aaa,'; max var = ',bbb)
        ###

        #fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title

        if fig is None:
            if Atl_Pac_section :
                if cube.ndim == 3 :
                    fig = plt.figure(figsize=(20,10))
                    plt.rcParams.update({'font.size': 16})
                elif cube.ndim == 2 :
                    fig = plt.figure(figsize=(10,20))
                    plt.rcParams.update({'font.size': 14})
            else : 
                fig = plt.figure(figsize=(20,16))
                plt.rcParams.update({'font.size': 18})

        if ax is None :
            if Proj == True :
                ax = fig.add_subplot(iii, jjj, rect, projection=projections[name])
                #fig.add_subplot(ax)
                #ax = plt.subplot(projection=projections[name])
                # Set limits
                if proj !='Amazon' or proj != 'Bengal' or proj != 'UKshelv' :
                    ax.set_global()
            else :
                ax = fig.add_subplot(iii, jjj, rect)
                new_cube = cube.copy()
        else :
            if Proj == True :
                boulet=1
                #ax.set_global()
            else :
                new_cube = cube.copy()

        if cmap :
            mycmap = cmap
        elif ('compar' in attrs) :
            mycmap = plt.cm.get_cmap("BrBG")# jet, seismic, spectral
        else :
            mycmap = plt.cm.get_cmap("jet")# jet, seismic, spectral

        # plot with Iris quickplot pcolormesh

        loc_in_attr = location is not None

        if colorbar and not loc_in_attr :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = qplt.contourf(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = qplt.pcolormesh(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                elif log_sc :
                    if ctour :
                        bb = qplt.contourf(new_cube, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                    else :
                        bb = qplt.pcolormesh(new_cube, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                else :
                    if ctour :
                        bb = qplt.contourf(new_cube, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = qplt.pcolormesh(new_cube, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            new_cube.data = np.log10(new_cube.data)
                            log_lev = np.log10(levels)
                            bb = qplt.contour(new_cube, log_lev, cmap=mycmap, axes = ax)
                        else:
                            bb = qplt.contour(new_cube, levels, cmap=mycmap, axes = ax)
                        ax.clabel(bb)
                    else :
                        bb = qplt.contourf(new_cube, cmap=mycmap, axes = ax)
                else :
                    bb = qplt.pcolormesh(new_cube, cmap=mycmap, axes = ax)
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())

            if BIOMASK :
                iplt.contour(cube_ref, colors = 'black', axes = ax)
            if ADDCTOUR :
                iplt.contour(ctour_ref, colors = 'black', axes = ax)


        else :

            if ('v_min' in attrs) and ('v_max' in attrs):
                v_min = attrs['v_min']
                v_max = attrs['v_max']
                # ax = plt.axes()
                ## cmap = BrBG

                if norm :
                    if ctour :
                        bb = iplt.contourf(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = iplt.pcolormesh(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max, axes = ax)
                elif log_sc :
                    if ctour :
                        bb = iplt.contourf(new_cube, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                    else :
                        bb = iplt.pcolormesh(new_cube, cmap=mycmap, norm=pltc.LogNorm(vmin=v_min, vmax=v_max), axes = ax)
                else :
                    if ctour :
                        bb = iplt.contourf(new_cube, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)
                    else :
                        bb = iplt.pcolormesh(new_cube, cmap=mycmap, vmin=v_min, vmax=v_max, axes = ax)

            else:
                # ax = plt.axes()
                if ctour :
                    if 'levels' in attrs:
                        levels = attrs['levels']
                        if log_sc :
                            new_cube.data = np.log10(new_cube.data)
                            log_lev = np.log10(levels)
                            bb = iplt.contour(new_cube, log_lev, cmap=mycmap, axes = ax)
                        else:
                            bb = iplt.contour(new_cube, levels, cmap=mycmap, axes = ax)
                        ax.clabel(bb)
                    else :
                        bb = iplt.contourf(new_cube, cmap=mycmap, axes = ax)
                else :
                    bb = iplt.pcolormesh(new_cube, cmap=mycmap, axes = ax)
                #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())

            if BIOMASK :
                iplt.contour(cube_ref, colors = 'black', axes = ax)

            if ADDCTOUR :
                iplt.contour(ctour_ref, colors = 'black', axes = ax)

            if location is not None :
                #location = attrs['location']
                #
                #divider = make_axes_locatable(ax)
                #cax = divider.append_axes(location, size="10%", pad=0.05)
                #cbar = fig.colorbar(bb, cax=cax)
                #
                #cbar = fig.colorbar(bb, ax=ax, location=location)
                if Proj : 
                    shrink=0.5
                else : 
                    shrink=1.0
                cbar = fig.colorbar(bb, ax=ax, location=location, shrink=shrink, aspect=10)
                #divider = make_axes_locatable(ax)
                ## creating new axes on the right
                ## side of current axes(ax).
                ## The width of cax will be 5% of ax
                ## and the padding between cax and ax
                ## will be fixed at 0.05 inch.
                #colorbar_axes = divider.append_axes(location,
                #                                    size="10%",
                #                                    pad=0.1)
                ## Using new axes for colorbar
                #cbar = plt.colorbar(bb, cax=colorbar_axes)

                ## colorbar label :
                if cube.units != "unknown" or cube.units != "none":  
                    cbar.set_label(cube.units)
                else :
                    var=cube.var_name
                    print("var : ",var)
                    print(cube.name)
                    print(cube.long_name)
                    if var == "TOT_SHALK" : 
                        units = "meq m-2"
                    elif var == "SIL" : 
                        units = "mmol-Si/m3"
                    elif var == "ALK" :
                        units = "meq/m3"
                    cbar.set_label(units)

                if 'remove_cb' in attrs :
                    cbar.remove()


        # Draw coastlines
        if COAST_LINE :
            try:
                ax.coastlines()
            except:
                print("no coastline to draw")
        ##
        ## Specific section Layout : 
        if Atl_Pac_section : 
            #yy1 = cube.coord(axis='y').points
            #yysize = yy1.shape[0]
            #rest = yysize - 180
            yy = np.arange(-180,180, 30 )
            yy_labels = ["90°","60°","30°","0°","-30°","-60°", "-90°","-60°","-30°", "0°","30°","60°"]

            ## test if z coord is here.  if it is, the plot is the Atl_Pac section.
            ##                         otherwise, it is the hovmoeller
            try :
                ZZZ = cube.coord(axis='z')
                if ZZZ is not None :
                    section=True
                else :
                    ## not z coord : it not a section
                    section=False
            except :
                ## not z coord : it not a section
                section=False


            if  section :
                #yy[181:] = np.arange(90,90 - (rest-1), -1)
                ax.xaxis.set_ticks(yy, yy_labels)
                secax = ax.secondary_xaxis('top')
                yy2 = np.arange(-180,180, 90 )
                yy2_labels = ["Arctic","Atlantic", "Antactica", "Pacific"]
                secax.xaxis.set_ticks(yy2, yy2_labels)
            else :
                #yy[181:] = np.arange(90,90 - (rest-1), -1)
                ax.yaxis.set_ticks(yy, yy_labels)
                secax = ax.secondary_yaxis('right')
                yy2 = np.arange(-180,180, 90 )
                yy2_labels = ["Arctic","Atlantic", "Antactica", "Pacific"]
                secax.yaxis.set_ticks(yy2, yy2_labels, rotation=-90, va='center')
        ##    
        if Epipel_plot == True :
            ax.set_ylim([500, 0.0])
        self.ax = ax         # Graphical axes
        if xrev :
            ax.yaxis.tick_right()
        if ytick_loc == "right" :
            ax2.invert_xaxis()


    def subtitle(self, title):
        self.ax.set_title(title)

    def y_label(self, title):
        self.ax.set_ylabel(title)
        #self.ax.text(-1.5, 0.5, title, horizontalalignment='left',
        #             verticalalignment='center', rotation=90 )

    def x_label(self, title):
        self.ax.set_xlabel(title)









def plot_proj_Orcagrid(cube, cmap = None, norm = None, BIOMASK = None, ADDCTOUR = None, **attrs) :
    """
    Tri-Polar Grid Projected Plotting
    =================================

    This example demonstrates cell plots of data on the semi-structured ORCA2 model
    grid.

    First, the data is projected into the PlateCarree coordinate reference system.

    Second four pcolormesh plots are created from this projected dataset,
    using different projections for the output image.
    --Projection list:
    Mollweide
    PlateCarree
    NorthPolarStereo
    SouthPolarStereo
    Orthographic
    SouthOrtho
    NorthOrtho
    """

    import matplotlib.pyplot as plt

    import cartopy.crs as ccrs
    import iris
    import iris.analysis.cartography
    import iris.plot as iplt
    import iris.quickplot as qplt


    def main():
        # Load data
        #filepath = iris.sample_data_path('orca2_votemper.nc')
        #cube = iris.load_cube(filepath)

        # Choose plot projections
        projections = {}
        projections['Mollweide'] = ccrs.Mollweide()
        projections['Mollweide_pac'] = ccrs.Mollweide(central_longitude=-170)
        projections['Mollweide_ind'] = ccrs.Mollweide(central_longitude= 90)
        projections['PlateCarree'] = ccrs.PlateCarree()
        projections['NorthPolarStereo'] = ccrs.NorthPolarStereo()
        projections['SouthPolarStereo'] = ccrs.SouthPolarStereo()
        projections['Orthographic'] = ccrs.Orthographic(central_longitude=-90,
                                                    central_latitude=45)
        projections['SouthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                    central_latitude=-90)
        projections['NorthOrtho'] = ccrs.Orthographic(central_longitude=0,
                                                    central_latitude=+90)
        projections['OrthoAtl'] = ccrs.Orthographic(central_longitude=-30,
                                                    central_latitude=15)
        projections['OrthoPac'] = ccrs.Orthographic(central_longitude=160,
                                                    central_latitude=7)
        projections['OrthoInd'] = ccrs.Orthographic(central_longitude=80,
                                                    central_latitude=-20)
        projections['Amazon']  = ccrs.Orthographic(central_longitude=2,
                                                    central_latitude=-50)
        projections['Bengal']  = ccrs.Orthographic(central_longitude=16,
                                                    central_latitude=90)
        projections['UKshelv'] = ccrs.Orthographic(central_longitude=55,
                                                    central_latitude=-2)
        pcarree = projections['PlateCarree']

        # Transform cube to target projection
        new_cube, extent = iris.analysis.cartography.project(cube, pcarree,
                                                         nx=400, ny=200)
        if BIOMASK :
            biocubes = iris.load(BIOMASK)
            biocube_ref = biocubes[0].copy()
            cube_ref, extent = iris.analysis.cartography.project(biocube_ref, pcarree,
                                                         nx=400, ny=200)
        if ADDCTOUR :
            ctour_ref, extent = iris.analysis.cartography.project(ADDCTOUR, pcarree,
                                                         nx=400, ny=200)

        if 'proj' in attrs:
            name=attrs['proj']
        else:
            name='Mollweide'
        ctour = attrs['contour'] if 'contour' in attrs else False
        log_sc = attrs['log_scale'] if 'log_scale' in attrs else False

        aaa = cube.data.min()
        bbb = cube.data.max()

        print('min var = ',aaa,'; max var = ',bbb)
        ###

        fig = plt.figure(figsize=(20,10))
        #fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title
        ax = plt.subplot(projection=projections[name])
        # Set limits
        ax.set_global()

        if cmap :
            mycmap = cmap
        elif ('compar' in attrs) :
            mycmap = plt.cm.get_cmap("BrBG")# jet, seismic, spectral
        else :
            mycmap = plt.cm.get_cmap("jet")# jet, seismic, spectral

        # plot with Iris quickplot pcolormesh

        if ('v_min' in attrs) and ('v_max' in attrs):
            v_min = attrs['v_min']
            v_max = attrs['v_max']
            # ax = plt.axes()
            ## cmap = BrBG
            if norm :
                if ctour :
                    qplt.contourf(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max)
                else :
                    qplt.pcolormesh(new_cube, cmap=mycmap, norm = norm, vmin=v_min, vmax=v_max)
            else :
                if ctour :
                    qplt.contourf(new_cube, cmap=mycmap,vmin=v_min, vmax=v_max)
                else :
                    qplt.pcolormesh(new_cube, cmap=mycmap,vmin=v_min, vmax=v_max)
        else:
            # ax = plt.axes()
            if ctour :
                if 'levels' in attrs:
                    levels = attrs['levels']
                    if log_sc :
                        new_cube.data = np.log10(new_cube.data)
                        log_lev = np.log10(levels)
                        bb = qplt.contour(new_cube, log_lev, cmap=mycmap)
                    else:
                        bb = qplt.contour(new_cube, levels, cmap=mycmap)
                    ax.clabel(bb)
                else :
                    qplt.contourf(new_cube, cmap=mycmap)
            else :
                qplt.pcolormesh(new_cube, cmap=mycmap)
            #bibi = ax.pcolormesh(Y, X, var, transform=ccrs.PlateCarree())

        if BIOMASK :
            iplt.contour(cube_ref, colors = 'black')

        if ADDCTOUR :
            iplt.contour(ctour_ref, colors = 'black')

        # Draw coastlines
        ax.coastlines()

        fnam = attrs['f_name'] if 'f_name' in attrs else 'test_plot'
        #ax.set_xlabel(Xax, fontsize=16)
        #ax.set_ylabel(Yax, fontsize=16)
        #plt.show()
        plt.savefig(fnam)

        iplt.show()

    if __name__ == '__main__':
        main()







def bbox_extract_2Dcoords(cube, ii_min=-10, ii_max=10, jj_min=-10, jj_max=10):
    """
    Extract a sub-set of a cube inside a lon, lat bounding box
    bbox=[lon_min lon_max lat_min lat_max].
    NOTE: This is a work around too subset an iris cube that has
    2D lon, lat coords.

    """
    minmax = lambda x: (np.min(x), np.max(x))

    lons = cube.coord('longitude').points
    lats = cube.coord('latitude').points
    #lons = wrap_lon180(lons)

    inregion = np.logical_and(np.logical_and(lons > ii_min,
                                             lons < ii_max),
                              np.logical_and(lats > jj_min,
                                             lats < jj_max))
    region_inds = np.where(inregion)
    imin, imax = minmax(region_inds[0])
    jmin, jmax = minmax(region_inds[1])
    return cube[..., imin:imax+1, jmin:jmax+1]

    #sub_cube = bbox_extract_2Dcoords(cube, bbox)

    #print('Original {} and new {} horizontal shapes'.format(cube.shape[1:],
    #                                                    sub_cube.shape[1:]))







