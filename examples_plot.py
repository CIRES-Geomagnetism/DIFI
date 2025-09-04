import datetime
import numpy as np
import matplotlib.pyplot as plt
from DIFI.getSQfield import getSQfield

#This example shows the DIFI output at Honolulu, Hawaii
#for a week in June 2015
single_geodetic_lat = 21.3099
single_lon = -157.8581 #negative for west longitude
single_altitude = 0.

#DIFI's date and time inputs are lists, floats, or numpy arrays
#(year,month,day,hours,minutes)

start_datetime = datetime.datetime(2015,6,1,0,0,0)
end_datetime = datetime.datetime(2015,6,7,0,0,0)

#Build up lists for each DIFI date and time input
#by adding one hour to the date until it reaches the end date
years,months,days = [],[],[]
hours = []
dtimes = [] #Also accumulate the datetimes for the plot x axis
dtime = start_datetime
n_times = 0
while dtime < end_datetime:
    n_times+=1
    dtime += datetime.timedelta(hours=1)
    dtimes.append(dtime)
    years.append(dtime.year)
    months.append(dtime.month)
    days.append(dtime.day)
    hours.append(dtime.hour)
    
# DIFI can take lists and floats as inputs but can also use
# numpy arrays, as shown below
lat_arr = np.ones((n_times,))*single_geodetic_lat
lon_arr = np.ones((n_times,))*single_lon
alt_arr = np.ones((n_times,))*single_altitude

year_arr = np.array(years)
month_arr = np.array(months)
day_arr = np.array(days)
hour_arr = np.array(hours)

f = plt.figure(figsize=(8,6), dpi=200)
axx = f.add_subplot(3,1,1)
axy = f.add_subplot(3,1,2)
axz = f.add_subplot(3,1,3)

for model_name in ['difi8','xdifi2']:

    Bxyz = getSQfield(lat_arr,
                        lon_arr,
                        year_arr,
                        month_arr,
                        day_arr,
                        hour = hour_arr,
                        h = alt_arr,
                        model_name=model_name)
    
    Bx = Bxyz['X']
    By = Bxyz['Y']
    Bz = Bxyz['Z']

    axx.plot(dtimes,Bx,label=f'{model_name}',marker='.')
    axy.plot(dtimes,By,label=f'{model_name}',marker='.')
    axz.plot(dtimes,Bz,label=f'{model_name}',marker='.')

axx.set_ylabel('Bx [nT]')
axy.set_ylabel('By [nT]')
axz.set_ylabel('Bz [nT]')

for ax in [axx,axy,axz]:
    ax.legend()
    ax.grid(True)

#format plot to use a single x axis label set
f.autofmt_xdate()
f.suptitle(f'DIFI Output Timeseries\n Latitude {single_geodetic_lat:.02f}\n Longitude {single_lon:.02f}')
f.savefig('example_plot.png')
plt.show()