import datetime as dt
from datetime import timedelta
import time
import os
import shutil
import requests
import warnings
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from bs4 import BeautifulSoup

import step1_download_NCEP as step1
import step2_download_CMEMS_copernicus_marine_client as step2
import step3_gen_tide_TPOX9_last_TMD as step3
import step4_gen_Inverted_Barometer as step4
import step5_make_wave_forcing as step5
import step6_make_wind_forcing as step6
import step7_run_SWAN as step7
import step8_postprocess_output as step8 
# import step9_make_flood_risk as step9
# import step10_run_XBeach as step10
# import step11_run_SFINCS as step11
# import step12_postprocess_SFINCS as step12


# import step10_ingest2GUI as step10


def list_available_runs(url):
    session = requests.Session()
    retry = Retry(connect=5, backoff_factor=1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount(url, adapter)

    try:
        req = session.get(url).text
        session.close() 
        soup = BeautifulSoup(req, 'html.parser')     
        x = (soup.find_all('a'))
        runs = []
        for i in x:
            file_name = i.extract().get_text()
            runs.append(int(file_name))
    except:
        
        runs = []
        print('Keep working on making the dowinloading process more robust')
        
      
    return(runs)



def delete_ndaysbefore(now,ndays):
    folder_name='../runs/'
    flist = os.listdir(folder_name)
    for d in flist:
        rundate=dt.datetime.strptime(d+'0000',"%Y%m%d%H%M%S")
        if rundate<now-timedelta(ndays):
            shutil.rmtree('../runs/'+ d)
    return()

def delete_nmonthsbefore(now,nmonths):
    folder_name='../archives/'
    flist = os.listdir(folder_name)
    for d in flist:
        nb = d.split("_")[0]
        rundate=dt.datetime.strptime(nb+'0000',"%Y%m%d%H%M%S")
        if rundate<now-timedelta(nmonths * 30):
            shutil.rmtree('../archives/'+ d)
    return()
            
    

###############################################################################


# Suppress only deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


#- Find the latest available run in nomads.ncep.noaa.gov
now = dt.datetime.utcnow()
url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
runs = list_available_runs(url)
if len(runs)==0:
    now = dt.datetime.utcnow()-timedelta(1)
    url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
    runs = list_available_runs(url)

#- Define the run to be used
runs=sorted(runs)
now = now.replace(hour=runs[-1],minute=0,second=0,microsecond=0)

#- Find the llatest available run in nomads.ncep.noaa.gov
# now = dt.datetime.utcnow()
# url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
# runs = list_available_runs(url)

# if not runs:
#     now = dt.datetime.utcnow()
#     now = now - dt.timedelta(days=1) 
#     url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
#     runs = list_available_runs(url)

    
# if len(runs)==1:
#     now = dt.datetime.utcnow()-timedelta(1)
#     url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
#     runs = list_available_runs(url)
#     #- Define the run to be used
#     now = now.replace(hour=runs[-1],minute=0,second=0,microsecond=0)
# else:
#     runs=sorted(runs)
#     now = now.replace(hour=runs[len(runs)-2],minute=0,second=0,microsecond=0)



# delete previous runs older than 14 days
try:
    delete_ndaysbefore(now,14)
except Exception as e:
    print(e)
    

#now = dt.datetime(2025, 2, 20, 12, 0)

forecast_path = '../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") 
try:
    os.mkdir(forecast_path)
except OSError as error:
    print(error)    


start_time = time.time()

step1.download_NCEP(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step2.download_CMEMS_copernicus_marine_client(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step3.gen_tide(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step4.gen_IB(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step5.make_waves(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step6.make_winds(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step7.run_SWAN(now)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
step8.postprocess_SWAN(now,1)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")
# step9.make_flood_points(now,1,1)
# elapsed_time = time.time() - start_time
# print(f"Elapsed time: {elapsed_time:.6f} seconds")
# step10.run_XBeach_profiles(now)
# elapsed_time = time.time() - start_time
# print(f"Elapsed time: {elapsed_time:.6f} seconds")
# step11.run_SFINCS(now)
# elapsed_time = time.time() - start_time
# print(f"Elapsed time: {elapsed_time:.6f} seconds")
# step12.plot_SFINCS(now)
# elapsed_time = time.time() - start_time
# print(f"Elapsed time: {elapsed_time:.6f} seconds")




# # delete previous archived swan outputs older than 3 months
# try:
#     delete_nmonthsbefore(now,3)
# except Exception as e:
#     print(e)
    
# step9.archive_output(now)
# step10.ingest2GUI(now)

