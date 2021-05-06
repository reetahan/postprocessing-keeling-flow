import xarray as xr
import pandas as pd
import numpy as np
import gc
import time
import sys
from tqdm import tqdm

# how to use this script?
# python get_csv_new.py 3586369 train 0 99
        
def get_aero_data(job_id, sc_id):
    p1 = "/data/keeling/a/zzheng25/d/partmc-mam4-high-lat/"
    p2 = "/scenarios/scenario_"
    p3 = "/out/urban_plume_var.nc"
    aero_nc = xr.open_dataset(p1+job_id+p2+sc_id+p3)
    aero_df = aero_nc[["tot_so4_conc", "tot_no3_conc", "tot_cl_conc", "tot_nh4_conc", "tot_msa_conc",
        "tot_aro1_conc", "tot_aro2_conc", "tot_alk1_conc", "tot_ole1_conc", "tot_api1_conc",
        "tot_api2_conc", "tot_lim1_conc", "tot_lim2_conc", "tot_co3_conc", "tot_na_conc",
        "tot_ca_conc", "tot_oin_conc", "tot_oc_conc", "tot_bc_conc", "tot_h2o_conc"]].to_dataframe().reset_index()
    
    aero_df["Mass_so4"] = aero_df["tot_so4_conc"]
    aero_df["Mass_bc"] = aero_df["tot_bc_conc"]
    aero_df["Mass_ncl"] = aero_df["tot_na_conc"]+aero_df["tot_cl_conc"]
    #aero_df["Mass_ncl"] = aero_df["tot_na_conc"]+aero_df["tot_cl_conc"]+aero_df["tot_so4_conc"]
    aero_df["Mass_dst"] = aero_df["tot_oin_conc"]
    aero_df["Mass_pom"] = aero_df["tot_oc_conc"]
    aero_df["Mass_soa"] = aero_df["tot_aro1_conc"]+aero_df["tot_aro2_conc"]+aero_df["tot_alk1_conc"]+aero_df["tot_ole1_conc"]+aero_df["tot_api1_conc"]+aero_df["tot_api2_conc"]+aero_df["tot_lim1_conc"]+aero_df["tot_lim2_conc"]
    
    return aero_df[["Mass_so4", "Mass_bc", "Mass_ncl", "Mass_dst", "Mass_pom",
                    "Mass_soa", "time"]]

# extract the environmental variables and gas species from partmc raw output
def single_scenario(job_id, sc_id, time_length):
    p1 = "/data/keeling/a/zzheng25/d/partmc-mam4-high-lat/"
    p2 = "/scenarios/scenario_"
    p3 = "/out/urban_plume_0001_"
    
    # define the lists
    time_ls = []
    rh_ls = []
    gas_ls = []
    temperature_ls = []
    sza_ls = []
    lat_ls = []
    
    file_ls = []
    sc_ls = []
    # define the flag for gas species
    gas_name_flag = True

    for i in range(1,1+time_length):
        file = "%08i" %i  #define the file number
        print(p1+job_id+p2+sc_id+p3+file+".nc")
        ds = xr.open_dataset(p1+job_id+p2+sc_id+p3+file+".nc")
        
        sc_ls.append(sc_id)
        file_ls.append(file)
        
        time_ls.append(np.array(ds["time"]))
        rh_ls.append(np.array(ds["relative_humidity"]))
        temperature_ls.append(np.array(ds["temperature"]))
        sza_ls.append(np.array(ds["solar_zenith_angle"]))
        lat_ls.append(np.array(ds["latitude"]))
        gas_ls.append(ds.gas_mixing_ratio.to_pandas())

        if gas_name_flag:
            gas_names_ls = ds.gas_species.names.split(",")
            gas_name_flag = False

        del ds
        gc.collect()

    df = pd.concat(gas_ls,axis=1).transpose().reset_index().drop(['index'],axis=1)
    df.columns = gas_names_ls
    
    df["sc_id"] = sc_ls
    df["time_id"] = file_ls
    df["time"] = np.array(time_ls)
    df["RELHUM"] = rh_ls
    df["T"] = temperature_ls
    df["SZA"] = sza_ls
    df["lat"] = lat_ls

    
    df["SOAG_SRF"] = (df["ARO1"]+df["ARO2"]+df["ALK1"]+df["OLE1"]+df["API1"]+df["API2"]+df["LIM1"]+df["LIM2"])
    df["DMS_SRF"] =  df["DMS"] 
    df["H2SO4_SRF"] =  df["H2SO4"]  
    df["H2O2_SRF"] = df["H2O2"]   
    df["SO2_SRF"] = df["SO2"]  
    df["O3_SRF"] = df["O3"] 

    df["NO_SRF"] = df["NO"] 
    df["NO3_SRF"] = df["NO3"] 
    df["NOX_SRF"] = df["NO_SRF"] + df["NO2"] + df["NO3_SRF"]
    df["CO_SRF"] = df["CO"] 
    df["C2H6_SRF"] = df["C2H6"] 
    df["ETH_SRF"] = df["ETH"] 
    df["OLET_PAR_SRF"] = df["OLET"] + df["PAR"] 
    df["3PAR_SRF"] = 3 * df["PAR"] 
    df["4PAR_SRF"] = 4 * df["PAR"] 
    df["5PAR_SRF"] = 5 * df["PAR"] 
    df["6PAR_SRF"] = 6 * df["PAR"] 
    df["TOL_SRF"] = df["TOL"] 
    df["XYL_SRF"] = df["XYL"] 
    df["CH3OH_SRF"] = df["CH3OH"] 
    df["ALD2_SRF"] = df["ALD2"] 
    df["AONE_SRF"] = df["AONE"] 


    df_final = df[["sc_id", "time_id", "time",
                   "SOAG_SRF", "DMS_SRF", "H2SO4_SRF", "O3_SRF", "H2O2_SRF", "SO2_SRF", 
                   "NO_SRF","NO3_SRF","NOX_SRF","CO_SRF","C2H6_SRF","ETH_SRF",
                   "OLET_PAR_SRF","3PAR_SRF","4PAR_SRF","5PAR_SRF","6PAR_SRF", "TOL_SRF",
                   "XYL_SRF","CH3OH_SRF","ALD2_SRF","AONE_SRF","T", "RELHUM", "SZA", "lat"]]
    
    return df_final

def multi_scenarios(job_id, sc_number_start, sc_number_end, time_length):
    vari_ls = []
    p1 = "/data/keeling/a/zzheng25/d/partmc-mam4-high-lat/"
    p2 = "/scenarios/scenario_"
    p3 = "/out/urban_plume_chi.nc"
    
    for i in range(sc_number_start, sc_number_end+1):
        sc_id = "%04i" %i  #define the file number
        print("scenario id:", sc_id)
        
        # get mixing state metrics
        msm_nc = xr.open_dataset(p1+job_id+p2+sc_id+p3)
        msm = msm_nc[['chi_abd']].to_dataframe().reset_index()
        
        # get environmental and gas variables
        gas_env = single_scenario(job_id, sc_id, time_length)
        print("\n")
        
        # get aero variables
        aero = get_aero_data(job_id, sc_id)
        
        # merge environmental, gas, and mixing state metrics within scenario
        gas_env_msm = pd.merge(gas_env, msm, how = "left", on = "time")
        gas_env_msm_aero = pd.merge(gas_env_msm, aero, how = "left", on = "time")
        
        vari_ls.append(gas_env_msm_aero)
        
        del msm_nc, msm, gas_env, aero, gas_env_msm, gas_env_msm_aero
        gc.collect()
        
    # concat the df from different scenarios    
    df_final = pd.concat(vari_ls).reset_index(drop=True) 
    
    return df_final
    
job_id = sys.argv[1]
csv_name = sys.argv[2]
sc_number_start_type = int(sys.argv[3])
sc_number_end_type = int(sys.argv[4])
df = multi_scenarios(job_id, sc_number_start=sc_number_start_type, sc_number_end=sc_number_end_type, time_length=25)
df_train = df[:int(len(df)/2)]
df_test = df[int(len(df)/2):]
df_train.to_csv("/data/keeling/a/rm10/d/partmc-mam4-high-lat/train_v3.csv",index=False)
df_test.to_csv("/data/keeling/a/rm10/d/partmc-mam4-high-lat/test_v3.csv",index=False)
