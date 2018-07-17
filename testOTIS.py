import otis_tide_pred as otp
import numpy as np
import matplotlib.pyplot as plt

modfile = './DATA/Model_haw'

dates = np.arange(np.datetime64('2001-04-03'),
                  np.datetime64('2001-05-03'), dtype='datetime64[h]' )

lon = np.array([198, 199])
lat = np.array([21, 19])

h, u, v = otp.tide_pred(modfile, lon, lat, dates, z=None,conlist=None)

fig, ax = plt.subplots()
ax.plot(dates, u, dates, v)
plt.show()
