yr = "YY"

def ptrcf(mm, YYYY="1859") :
    return "eORCA1_MED_UKESM_1m_"+YYYY+"0101_"+YYYY+"1230_ptrc_T_"+YYYY+mm+"-"+YYYY+mm+".nc"

def diadf(mm, YYYY="1859") :
    return "eORCA1_MED_UKESM_1m_"+YYYY+"0101_"+YYYY+"1230_diad_T_"+YYYY+mm+"-"+YYYY+mm+".nc"

def gridf(mm, YYYY="1859") :
    return "eORCA1_MED_UKESM_1m_"+YYYY+"0101_"+YYYY+"1230_grid_T_"+YYYY+mm+"-"+YYYY+mm+".nc"


def rundict_file(yr, runlist,rundict):
    '''
    function to get all run-to-be-compared's files name and path.
    --
    !! ==> Function to be amended for all new experiments !! 
    --
    rundict_file is called with : 
    [rundict_ptrc,rundict_diad,rundict_grid] = rundict_file(yr, runlist, rundict)  
    where : 
    -- yr is the year to be used in the file names
    -- runlist is a list of names to be given to the different runs to be compared.
    the names are used as key in the dictionaries, and as run name in the plots.
    ex :  
    runlist = ["Obs",
           "No-River",
           "River-Flx",
           "River-conc",
           "OMIP-NO_RIV",
           "OMIP_RIV_FLX",
           "No-River+10",
           "No-River+10_hybrid",
           "OMIP_NO_RIV+10",
           "OMIP_NO_RIV+10_hybrid"]
    -- rundict is a dictionary with the path where to find the runs nc-files (medusa's ptrc and diad, and nemo grid T files) we'll read later.
    ex : 
    rundict = {
      runlist[0] : "OBS",
      runlist[1] : "NO_RIV",
      runlist[2] : "RIV_FLX",
      runlist[3] : "RIV_CONC",
      runlist[4] : "OMIP_NO_RIV",
      runlist[5] : "OMIP_RIV_FLX",
      runlist[6] : "NO_RIV_+10",
      runlist[7] : "NO_RIV_+10_hybrid",
      runlist[8] : "OMIP_NO_RIV_+10",
      runlist[9] : "OMIP_NO_RIV_+10_hybrid"
    }
    -- then the function returns 3 dictionaries 
    rundict_ptrc,rundict_diad,rundict_grid 
    with all experiments run path and file names for the 
       - ptrc files (main MEDUSA outputs)
       - diad files (MEDUSA diagnostics )
       - grid files (nemo grid T files with T, S and mld,.... variables.
    '''
    import glob
    rundict_ptrc = {
       runlist[0] : rundict[runlist[0]] +"/"+ "eORCA100-4.2_tracers.nc",
       #runlist[1] : rundict[runlist[1]] +"/"+ "ORCA1_MED42_UKESM_1y_19750101_19791230_ptrc_T_1975-1975.nc", #clim_2000s_ptrcT.nc",
       runlist[1] : glob.glob(rundict[runlist[1]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[2] : glob.glob(rundict[runlist[2]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[3] : glob.glob(rundict[runlist[3]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[4] : glob.glob(rundict[runlist[4]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[5] : glob.glob(rundict[runlist[5]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[6] : glob.glob(rundict[runlist[6]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[7] : glob.glob(rundict[runlist[7]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[8] : glob.glob(rundict[runlist[8]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_ptrcT.nc", 
       runlist[9] : glob.glob(rundict[runlist[9]] +"/"+ "ORCA1_MED42_UKESM_1y_*_ptrc_T_"+str(yr)+"-"+str(yr)+".nc") #clim_2000s_ptrcT.nc",  
    }
    ##
    rundict_diad = {
       runlist[0] : rundict[runlist[0]] +"/"+ "eORCA100-4.2_NPP_CHL_tracers.nc",
       runlist[1] : glob.glob(rundict[runlist[1]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc",
       runlist[2] : glob.glob(rundict[runlist[2]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc",
       runlist[3] : glob.glob(rundict[runlist[3]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc"
       runlist[4] : glob.glob(rundict[runlist[4]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc"
       runlist[5] : glob.glob(rundict[runlist[5]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc",
       runlist[6] : glob.glob(rundict[runlist[6]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc",
       runlist[7] : glob.glob(rundict[runlist[7]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc"
       runlist[8] : glob.glob(rundict[runlist[8]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_diadT.nc"
       runlist[9] : glob.glob(rundict[runlist[9]] +"/"+ "ORCA1_MED42_UKESM_1y_*_diad_T_"+str(yr)+"-"+str(yr)+".nc") #clim_2000s_diadT.nc" 
    }
    ##
    rundict_grid = {
       runlist[0] : glob.glob(rundict[runlist[0]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc",
       runlist[1] : glob.glob(rundict[runlist[1]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc"
       runlist[2] : glob.glob(rundict[runlist[2]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc",
       runlist[3] : glob.glob(rundict[runlist[3]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc"
       runlist[4] : glob.glob(rundict[runlist[4]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc",
       runlist[5] : glob.glob(rundict[runlist[5]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc"
       runlist[6] : glob.glob(rundict[runlist[6]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc",
       runlist[7] : glob.glob(rundict[runlist[7]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc"
       runlist[8] : glob.glob(rundict[runlist[8]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc"), #clim_2000s_gridT.nc",
       runlist[9] : glob.glob(rundict[runlist[9]] +"/"+ "ORCA1_MED42_UKESM_1y_*_grid_T_"+str(yr)+"-"+str(yr)+".nc") #clim_2000s_gridT.nc
    }
    return rundict_ptrc,rundict_diad,rundict_grid


