# OSU Tidal predictions in python

Gary Egbert's group at Oregon State University produce  inversions of
open-ocean tidal signals to give predictions of ocean tides.  There are
a variety of solutions available at http://volkov.oce.orst.edu/tides/otps.html

This software is based on a chain of software to decode these files.  The
files themselves are straightforward, but to get accurate tidal predictions
the amplitudes and phases they provide must be corrected by astronomical
tidal constants that depend on time.  Hence, there is a bit of extra
processing that is necessary (and somewhat arcane), provided by this package.

This package is extracted from https://github.com/ofringer/suntanspy and
work that Matt Rayson (@mrayson) did to code up the older matlab code.  The
only changes are that we dropped `suntanspy`-specific libraries for
interpolation in favour of `scipy.interpolate`, and use standard
`datetime` functions to handle the dates.  

## Installing

This doesn't have a `PyPI` package yet, but its pretty easy to install into
a `pip` or `conda` environment (if you aren't using `pip`, I'm not sure
how to easily install, but putting `otis_tide_pred.py` in your local
directory will work)

```
git clone https://github.com/jklymak/OtisTidalPrediction.git
cd OtisTidalPrediction
pip install -e .
```

## Using

Once installed, you should be able to `import otis_tide_pred as otp`.  

An example usage is in `testOTIS.py`.

```python
import otis_tide_pred as otp
import numpy as np
import matplotlib.pyplot as plt

modfile = './DATA/Model_haw'

dates = np.arange(np.datetime64('2001-04-03'),
                  np.datetime64('2001-05-03'), dtype='datetime64[h]' )

lon = np.array([198, 199, ])
lat = np.array([21, 19])

h, u, v = otp.tide_pred(modfile, lon, lat, dates, z=None,conlist=None)

fig, ax = plt.subplots()
ax.plot(dates, u, dates, v)
plt.show()
```

Note that we have had to modify the paths in the file `./Data/Model_haw`.
This is an odd annoyance of the OSU data packaging.  
