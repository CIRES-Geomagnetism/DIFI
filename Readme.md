# DIFI-8 IONOSPHERE MAGNETIC FIELD MODEL

The Dedicated Ionospheric Field Inversion (DIFI) model is a Swarm-based, global model of the Sq and equatorial electrojet magnetic fields at mid and low-latitudes. It describes variations with local time, season and solar flux and separates primary and induced magnetic fields.
The latest version of DIFI is DIFI-8. Please see the [DIFI-8](https://geomag.colorado.edu/index.php/node/785) for the detail. Find info on xDIFI2 [here](https://geomag.colorado.edu/index.php/node/786).

## Prerequirements

- Officially Python 3.9, 3.10, 3.11. 3.12 3.13
- geomaglib >= 1.2.1
  - install by `pip install geomaglib>=1.2.1`
- numpy 


## Installation
Our team will release the DIFI Python API to PyPI soon. If you would like to install DIFI at this time, please clone our git repository or use the following instructions to install DIFI in your virtual environment.

`pip install git+https://github.com/CIRES-Geomagnetism/DIFI.git@main`


## Quick Start

### Get magnetic vectors from single point

```python
from DIFI import getSQfield
lat, lon, year, month, day, hour = 20.5, 100.5, 2024, 6, 6, 0
B = getSQfield(lat, lon, year, month, day, hour=hour)
```

It will return the results in dict type
```python
{'Z': array([13.19036336]), 'Y': array([28.73130052]), 'X': array([0.65322314])}
```
### Get magnetic vectors from array inputs 

```python
import random
from DIFI import getSQfield

year = 2023
month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
day = 15

N = len(month)
lat = [random.uniform(-90, 90) for _ in range(N)]
lon = [random.uniform(-180, 180) for _ in range(N)]

B = getSQfield(lat, lon, year, month, day, model_name="xdifi2")
print(B)
```
It will return the results
```python
{'Z': array([-2.89438186, 18.22932278, -0.66058193,  7.91297099,  8.00560867,
       -2.84424521,  0.28690831,  5.21360212,  1.30044321,  5.78345537,
       -0.15228731,  0.17291943]), 
 'Y': array([  2.47943185,  36.00415516,   1.29842625,  22.99920884,
         4.12382535,   2.21563329,  -5.13570097,   6.71462628,
        -9.61616475,   0.52689522, -10.66703289,  -0.21119736]), 
 'X': array([ 2.77922129, 45.06151382,  1.14500601, -9.76631768, -1.57326838,
        4.39040841, -7.6578116 , 11.15431048, -5.66268463, -6.03204388,
       -4.6054734 ,  4.19473274])}
```

