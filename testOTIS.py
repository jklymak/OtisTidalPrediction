import otis_tide_pred as otp
import numpy as np
import matplotlib.pyplot as plt
import datetime

modfile = './DATA/Model_haw'


# test ordinal dates
dates = [datetime.datetime(2001, 4, 3),
         datetime.datetime(2001, 4, 4),
         datetime.datetime(2001, 4, 5)]

lon = np.array([198, 199, ])
lat = np.array([21, 19])

dates = np.array([datetime.datetime.toordinal(d) for d in dates])
h, u, v = otp.tide_pred(modfile, lon, lat, dates,conlist=None)


# test datetime dates
dates = [datetime.datetime(2001, 4, 3),
         datetime.datetime(2001, 4, 4),
         datetime.datetime(2001, 4, 5)]
lon = np.array([198, 199, ])
lat = np.array([21, 19])

h, u, v = otp.tide_pred(modfile, lon, lat, dates,conlist=None)

# test numpy datetime64 dates
dates = np.arange(np.datetime64('2001-04-03'),
                  np.datetime64('2001-05-03'), dtype='datetime64[h]' )

lon = np.array([198, 199, ])
lat = np.array([21, 19])

h, u, v = otp.tide_pred(modfile, lon, lat, dates,conlist=None)

print(np.shape(h))
fig, ax = plt.subplots()
ax.plot(dates, u, dates, v)
plt.show()
