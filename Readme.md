# DIFI-7 IONOSPHERE MAGNETIC FIELD MODEL

The Dedicated Ionospheric Field Inversion (DIFI) model is a Swarm-based, global model of the Sq and equatorial electrojet magnetic fields at mid and low-latitudes. It describes variations with local time, season and solar flux and separates primary and induced magnetic fields.
The latest version of DIFI is DIFI-7. Please see the [DIFI-7](https://geomag.colorado.edu/difi-7) for the detail.

## Installation
Our team will release to PyPI soon. If you would like to install DIFI at this time, please clone our git repository or use the following instructions to install DIFI in your virtual environment.

1. git clone the library
2. At the top folder in step 1, use `pip install .` to install the DIFI to your virtual environment

## Run the DIFI Module with Python

DIFI 7 is developed and tested in Python 3.11.2 environment

## Quick Start

### Get magnetic vectors from single point

```python
from DIFI import getSQfield
lat, lon, year, month, day, hour = 20.5, 100.5, 2024, 6, 6, 0
B = getSQfield(lat, lon, year, month, day, hour=hour)
```

It will return the results in dict type
```python
{'Z': array([10.59279794]), 'Y': array([28.49755808]), 'X': array([-1.42726185])}
```
### Get magnetic vectors from array inputs 

```python
from DIFI import getSQfield
year = 2023
month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
day = 15

lat = [80, 81, 82, 83, 84, 85, 86, 87, 88, 89]
lon = [163.0, 163.1, 163.2, 163.3, 163.4, 163.5, 163.6, 163.7, 163.8, 163.9]

B = getSQfield(lat, lon, year, month, day)
```
It will return the results
```python
{'Z': array([-1.11753492, -2.17806267, -2.36170017, -1.55958987, -0.4729657 ,
        0.02314406, -0.62739855, -1.84115162, -2.36710727, -1.87638184,
       -0.8895632 , -0.48355758]), 
 'Y': array([-0.53981766,  0.02365691,  1.0244612 ,  2.68627321,  4.66764832,
        6.15517171,  5.83115627,  3.98613964,  1.96799105,  0.7158783 ,
       -0.06700398, -0.52570463]), 
 'X': array([-1.63857938, -0.99863218, -0.70018449, -1.01992701, -1.70769538,
       -2.239885  , -2.03735071, -1.2146875 , -0.52421696, -0.60781745,
       -1.30196673, -1.82828954])}
```

