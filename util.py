# %%
import numpy as np
import math
import xarray as xr
from global_land_mask import globe
# %% Subfunction
# This function is used for modify the gas emission at certain row
def modify_gas_emit_func(flist_row, coefficient):
    flist_split_list = flist_row.split()
    flist_split_np=np.array(flist_split_list[1:],dtype=np.float64)
    #print(flist_split_np)
    
    # multiply by coefficient
    flist_split_np=np.array(flist_split_np * coefficient,dtype=np.str)
    #print(flist_split_np)
    
    # convert back to list
    flist_split_new_list_np = list(flist_split_np)
    flist_split_new_list = [str.upper("{:.4}".format(float(x))) for x in flist_split_new_list_np]

    # add the "\n" and the begining element of the original list (e.g., rate)
    flist_split_new_list.append("\n")
    flist_split_new_list.insert(0,flist_split_list[0])
    
    # join back the list
    flist_row_new = ' '.join(flist_split_new_list)
    return flist_row_new 

# This function is used for modify the temperature at a constant profile
def modify_temp_func(doy,year,lat, ss_option):
    doy = int(float(doy)) 
    year = int(float(year))
    year = 2009 + year
    lat = round_latlong(lat)
    lon = get_longitude(lat, ss_option)
    print(doy)
    print(year)
    print(lat)
    print(lon)

    is_leap = False
    if((year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)):
        is_leap = True
    '''
    print('entering modify_temp_func')
    flist_split_list = flist_row.split()
    print('flist_split_list:')
    print(flist_split_list)
    print(len(flist_split_list))
    flist_split_np=np.array(flist_split_list[1:],dtype=np.float64)
    #print(flist_split_np)
    
    # multiply by coefficient
    print(coefficient)
    flist_split_np=np.array(np.repeat(coefficient,len(flist_split_np)),dtype=np.str)
    print(flist_split_np)
    
    # convert back to list
    flist_split_new_list = list(flist_split_np)

    # add the "\n" and the begining element of the original list (e.g., rate)
    flist_split_new_list.append("\n")
    flist_split_new_list.insert(0,flist_split_list[0])
    
    # join back the list
    flist_row_new = ' '.join(flist_split_new_list)
    print('final flist_row')
    print(flist_row_new)
    '''

    p1 = '../data/'
    p3 = '_data.nc'

    filename = p1 + str(year) + p3
    ds = xr.open_mfdataset(filename)
    ds = ds['t2m']
    print('Grabbing temps')
    
    data = ds
    adj_lat = (90 - int(lat))*4
    adj_lon = -1
    if(lon < 0):
        adj_lon = (180 + (180 - (-1 * int(lon)))) * 4
    else:
        adj_lon = int(lon) * 4
    
    '''
    months_year = [31,28,31,30,31,30,31,31,30,31,30,31]
    date = 0
    for i in range(month-1):
        date = date + months_year[i]
    date = date + day
    '''
    date = doy 
    time_zone = round(lon/15)
    print('Calculated date, time')
    mean_list = None
    '''
    if(time_zone == 0):
        mean = data.isel(days=date, lats=adj_lat, longs=adj_lon)
        #print(mean["__xarray_dataarray_variable__"].values)
        mean_list = mean["__xarray_dataarray_variable__"].values
        #mean_arr = mean.to_array(dim="added")
        #mean_arrs = mean_arr["added"]
    elif(time_zone > 0):
        mean_1 = data.isel(days=date-1, hours=np.arange(24-time_zone,24), lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date, hours=np.arange(0,24-time_zone) , lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)
    else:
        mean_1 = data.isel(days=date, hours=np.arange(-1*time_zone,24),lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date+1, hours=np.arange(0,-1*time_zone),lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)
    '''
    #print(ds)
    #print('=======')
    #print(date)
    #print(adj_lat)
    #print(adj_lon)
    #print(time_zone)
    mean_list = get_temp_profile(ds,date,year,adj_lat,adj_lon,time_zone,is_leap)
    final_str = 'temp '
    for entry in mean_list:
        final_str = final_str + str(entry) + ' '
    final_str = final_str + '\n'
    return final_str
    
def mt_np_arr():
    a = np.array([0])
    a = np.delete(a, 0)
    return a

def round_latlong(val):
    val = float(val)
    mult = 0.25
    rem = val % mult
    if(rem >= mult/2):
        val = val + mult - rem
    else:
        val = val - rem
    return val

def get_longitude(lat, ss_option):
    lat = float(lat)
    lon = None
    land_list = mt_np_arr()
    sea_list = mt_np_arr()

    for i in np.arange(-180,180,0.25):
        if(globe.is_land(lat,i)):
            land_list = np.append(land_list, i)
        else:
            sea_list = np.append(sea_list, i)

    print(land_list)
    print(sea_list)
    if(ss_option == None):
        lon = np.random.choice(land_list)
    else:
        if(np.random.random_sample() < 2.0/3.0):
            lon = np.random.choice(sea_list)
        else:
            lon = np.random.choice(land_list)
    return lon

def get_temp_profile(data,doy,year,lat,lon,tz,is_leap):

    last_day = 365
    if is_leap:
        last_day = 366
    day_info = None

    if (doy == 0 and tz > 6):
        p1 = '../data/'
        p3 = '_data.nc'
        filename = p1 + str(year-1) + p3
        ds_last = xr.open_mfdataset(filename)
        ds_last = ds['t2m']
        amt = tz - 6
        indices_1 = np.arange((last_day-1)*24 + (24-amt),(last_day-1)*24 + 24)
        indices_2 = np.arange(0, 24-amt)
        day_info_1 = ds_last.isel(time=indices_1,latitude=lat,longitude=lon)
        day_info_2 = data.isel(time=indices_2,latitude=lat,longitude=lon)
        day_info = day_info_1 + day_info_2


    elif (doy == last_day and tz < 6):
        p1 = '../data/'
        p3 = '_data.nc'
        filename = p1 + str(year+1) + p3
        ds_next = xr.open_mfdataset(filename)
        ds_next = ds['t2m']
        amt = 6 -tz
        indices_1 = np.arange((last_day-1)*24 + amt, (last_day-1)*24 + 24)
        indices_2 = np.arange(0, 24-amt)
        day_info_1 = data.isel(time=indices_1,latitude=lat,longitude=lon)
        day_info_2 = ds_next.isel(time=indices_2,latitude=lat,longitude=lon)
        day_info = day_info_1 + day_info_2

    else:
        #if(time_zone == 0):
        indices = np.arange(24*(doy-1)+6-tz,24*doy+6-tz)
        #mean = data.isel(days=date, lats=adj_lat, longs=adj_lon)
        print(indices)
        print(lat)
        print(lon)
        day_info = data.isel(time=indices,latitude=lat,longitude=lon)
        #print(mean["__xarray_dataarray_variable__"].values)
        
        #mean_arr = mean.to_array(dim="added")
        #mean_arrs = mean_arr["added"]
    day_info = day_info.data
    '''
    elif(time_zone > 0):
        mean_1 = data.isel(days=date-1, hours=np.arange(24-time_zone,24), lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        day_info = data[time=doy,latitude=lat,longitude=lon]
        day_info = day_info["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date, hours=np.arange(0,24-time_zone) , lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)
    else:
        mean_1 = data.isel(days=date, hours=np.arange(-1*time_zone,24),lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date+1, hours=np.arange(0,-1*time_zone),lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)
    day_info = data[time=doy,latitude=lat,longitude=lon]
    '''
    return day_info

    
# %% function aero_emit
def modify_aero_emit_comp_carbo(directory, matrix):
    f=open(directory+"/aero_emit_comp_carbo.dat", "r+")
    flist=f.readlines()
    # modify the matrix here
    flist[1] = "BC              " + "{:.4}".format(matrix[24]) + "\n"
    flist[2] = "OC              " + "{:.4}".format(1-matrix[24]) + "\n"    
    f=open(directory+"/aero_emit_comp_carbo.dat", "w+")
    f.writelines(flist)
    f.close()
    
def modify_aero_emit_comp_ss1(directory, matrix):
    f=open(directory+"/aero_emit_comp_ss1.dat", "r+")
    flist=f.readlines()
    # modify the matrix here
    flist[1] = "OC              " + "{:.4}".format(matrix[28]) + "\n"
    flist[2] = "Na              " + "{:.4}".format((1-matrix[28])*0.3856) + "\n"
    flist[3] = "Cl              " + "{:.4}".format((1-matrix[28])*0.5389) + "\n"
    flist[3] = "SO4             " + "{:.4}".format((1-matrix[28])*0.0755) + "\n"
    f=open(directory+"/aero_emit_comp_ss1.dat", "w+")
    f.writelines(flist)
    f.close()

def modify_aero_emit_comp_ss2(directory, matrix):
    # modify the filename here
    f=open(directory+"/aero_emit_comp_ss2.dat", "r+")
    flist=f.readlines()
    # modify the matrix here
    flist[1] = "OC              " + "{:.4}".format(matrix[32]) + "\n"
    flist[2] = "Na              " + "{:.4}".format((1-matrix[32])*0.3856) + "\n"
    flist[3] = "Cl              " + "{:.4}".format((1-matrix[32])*0.5389) + "\n"
    flist[3] = "SO4             " + "{:.4}".format((1-matrix[32])*0.0755) + "\n"
    f=open(directory+"/aero_emit_comp_ss2.dat", "w+")
    f.writelines(flist)
    f.close()


def modify_aero_emit_dist(directory, matrix, ss_option, dust_option):
    f=open(directory+"/aero_emit_dist.dat", "r+")
    flist=f.readlines()
    
    # carbo
    flist[4] = "num_conc " + "{:.4}".format(matrix[23]) + "                     # particle number density (#/m^2)\n"
    flist[5] = "geom_mean_diam " + "{:.4}".format(matrix[21]) + "                # geometric mean diameter (m)\n"
    flist[6] = "log10_geom_std_dev " + "{:.4}".format(math.log10(matrix[22])) + "           # log_10 of geometric std dev of diameter\n"
    
    if ss_option != None:
        # ss1
        flist[12] = "num_conc " + "{:.4}".format(matrix[27]) + "                     # particle number density (#/m^2)\n"
        flist[13] = "geom_mean_diam " + "{:.4}".format(matrix[25]) + "                # geometric mean diameter (m)\n"
        flist[14] = "log10_geom_std_dev " + "{:.4}".format(math.log10(matrix[26])) + "           # log_10 of geometric std dev of diameter\n"
        
        # ss2
        flist[20] = "num_conc " + "{:.4}".format(matrix[31]) + "                     # particle number density (#/m^2)\n"
        flist[21] = "geom_mean_diam " + "{:.4}".format(matrix[29]) + "                # geometric mean diameter (m)\n"
        flist[22] = "log10_geom_std_dev " + "{:.4}".format(math.log10(matrix[30])) + "           # log_10 of geometric std dev of diameter\n"
    else:
        for ii in range(8, 24):
            flist[ii] = ""
    
    if dust_option != None:
        # dust1
        flist[28] = "num_conc " + "{:.4}".format(matrix[35]) + "                     # particle number density (#/m^2)\n"
        flist[29] = "geom_mean_diam " + "{:.4}".format(matrix[33]) + "                # geometric mean diameter (m)\n"
        flist[30] = "log10_geom_std_dev " + "{:.4}".format(math.log10(matrix[34])) + "           # log_10 of geometric std dev of diameter\n"
        
        # dust2
        flist[36] = "num_conc " + "{:.4}".format(matrix[38]) + "                     # particle number density (#/m^2)\n"
        flist[37] = "geom_mean_diam " + "{:.4}".format(matrix[36]) + "                # geometric mean diameter (m)\n"
        flist[38] = "log10_geom_std_dev " + "{:.4}".format(math.log10(matrix[37])) + "           # log_10 of geometric std dev of diameter\n"
    else:
        for ii in range(24,40):
            flist[ii] = ""

    f=open(directory+"/aero_emit_dist.dat", "w+")
    f.writelines(flist)
    f.close()

# %% function gas
def modify_gas_emit(directory, matrix, DMS_option):
    f=open(directory+"/gas_emit.dat", "r+")
    flist=f.readlines()
    for i in range(7,24):
        flist[i]=modify_gas_emit_func(flist[i], matrix[i-3])
    if DMS_option != None:
        flist[24]=modify_gas_emit_func(flist[24], matrix[40])
    else:
        flist[24]=""
    f=open(directory+"/gas_emit.dat", "w+")
    f.writelines(flist)
    f.close()

# %% function environment
def modify_temp(directory, matrix, ss_option):
    print('Entering modify_temp with params:')
    print('directory:')
    print(directory)
    print('matrix:')
    print(matrix)
    f=open(directory+"/temp.dat", "r+")
    flist=f.readlines()
    print(flist)
    flist[3]=modify_temp_func("{:.6}".format(matrix[2]),"{:.6}".format(matrix[3]), 
        "{:.6}".format(matrix[1]), ss_option)
    flist[2] = 'time 0\t3600\t7200\t10800\t14400\t18000\t21600\t25200\t28800\t32400\t36000\t39600\t43200\t46800\t50400\t54000\t57600\t61200\t64800\t68400\t72000\t75600\t79200\t82800\t86400\n'
    f=open(directory+"/temp.dat", "w+")
    print('flist!')
    print(flist)
    f.writelines(flist)
    f.close()
    print('exiting modify_temp')

# %% Make spec file
def make_spec(directory, scenario_num, matrix):
    f=open(directory+"/urban_plume_init.spec", "r+")
    flist=f.readlines()   
   
    # modify the matrix here
    # flist[1] = "output_prefix out/urban_plume_" + str(scenario_num).zfill(4) + "   # prefix of output files \n"
    
    flist[27] = "rel_humidity " + "{:.4}".format(matrix[0]) + "               # initial relative humidity (1) \n"
    print("Relative Humidity", "{:.4}".format(matrix[0]))
   
    flist[28] = "latitude  " + "{:.4}".format(matrix[1]) + "                      # latitude (degrees, -90 to 90) \n"
    print("latitude  " + "{:.4}".format(matrix[1]))
    #flist[29] = "longitude " + "{:.4}".format(matrix[0]) + "                     # longitude (degrees, -180 to 180) \n"   
   
    flist[32] = "start_day " + str(int(matrix[2])) + "                   # start day of year (UTC) \n"
    print("start_day " + str(int(matrix[2])))
    
    f=open(directory+"/urban_plume_init.spec", "w+")
    f.writelines(flist)
    f.close()

# %% Make spec file
def make_spec_restart(directory, scenario_num, matrix):
    f=open(directory+"/urban_plume_restart.spec", "r+")
    flist=f.readlines()   
   
    # modify the matrix here
    # flist[1] = "output_prefix out/urban_plume_" + str(scenario_num).zfill(4) + "   # prefix of output files \n"
    restart_time_stamp = "%08i" % int(matrix[39])
    
    flist[5] = "restart_file out_init/urban_plume_0001_" + restart_time_stamp + ".nc \n"
    print("restart_file out_init/urban_plume_0001_" + restart_time_stamp + ".nc \n")
   
    flist[21] = "rel_humidity " + "{:.4}".format(matrix[0]) + "               # initial relative humidity (1) \n"
    print("Relative Humidity", "{:.4}".format(matrix[0]))
   
    flist[22] = "latitude  " + "{:.4}".format(matrix[1]) + "                      # latitude (degrees, -90 to 90) \n"
    print("latitude  " + "{:.4}".format(matrix[1]))
    #flist[29] = "longitude " + "{:.4}".format(matrix[0]) + "                     # longitude (degrees, -180 to 180) \n"   
   
    flist[26] = "start_day " + str(int(matrix[2])) + "                   # start day of year (UTC) \n"
    print("start_day " + str(int(matrix[2])))    


    f=open(directory+"/urban_plume_restart.spec", "w+")
    f.writelines(flist)
    f.close()

def gen_temp_func(month,day, lat, lon):
    is_on_land = globe.is_land(lat, lon)
    land_str = ""
    if(is_on_land):
        land_str = "(Point on Land)"
    else:
        land_str = "(Point on Water)"

    filename = 'global_temp_avgs.nc'
    ds = xr.open_mfdataset(filename)
    print('Grabbing temps')
    
    data = ds
    adj_lat = (90 - lat)*4
    adj_lon = -1
    if(lon < 0):
        adj_lon = (180 + (180 - (-1 * lon))) * 4
    else:
        adj_lon = lon * 4
    
    months_year = [31,28,31,30,31,30,31,31,30,31,30,31]
    date = 0
    for i in range(month-1):
        date = date + months_year[i]
    date = date + day
    
    time_zone = round(lon/15)
    print('Calculated date, time')
    mean_list = None
    if(time_zone == 0):
        mean = data.isel(days=date, lats=adj_lat, longs=adj_lon)
        #print(mean["__xarray_dataarray_variable__"].values)
        mean_list = mean["__xarray_dataarray_variable__"].values
        #mean_arr = mean.to_array(dim="added")
        #mean_arrs = mean_arr["added"]
    elif(time_zone > 0):
        mean_1 = data.isel(days=date-1, hours=np.arange(24-time_zone,24), lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date, hours=np.arange(0,24-time_zone) , lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)
    else:
        mean_1 = data.isel(days=date, hours=np.arange(-1*time_zone,24),lats=adj_lat, longs=adj_lon)
        mean_list_1 = mean_1["__xarray_dataarray_variable__"].values
        mean_2 = data.isel(days=date+1, hours=np.arange(0,-1*time_zone),lats=adj_lat, longs=adj_lon)
        mean_list_2 = mean_2["__xarray_dataarray_variable__"].values
        print('Patching parts of days together')
        mean_list = mt_np_arr()
        mean_list = np.append(mean_list,mean_list_1)
        mean_list = np.append(mean_list,mean_list_2)

    return mean_list