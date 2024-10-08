# DIFI-7 IONOSPHERE MAGNETIC FIELD MODEL

The Dedicated Ionospheric Field Inversion (DIFI) model is a Swarm-based, global model of the Sq and equatorial electrojet magnetic fields at mid and low-latitudes. It describes variations with local time, season and solar flux and separates primary and induced magnetic fields.
The latest version of DIFI is DIFI-7. Please see the [DIFI-7](https://geomag.colorado.edu/difi-7) for the detail.

## Run the DIFI Module with Python

- DIFI 7 is developed and tested in Python 3.11.2 environment

## Install

```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ DIFI==1.1.0
```
TestPyPI only have numpy==1.9, so we need the flag `--extra-index-url` to search the latest versions of numpy frm PyPI
## Quick Start

### Get magnetic vectors from single point

```python
from DIFI import getSQfield
lat, lon, year, month, day, hour = 20.5, 100.5, 2024, 6, 6, 0
B = getSQfield(lat, lon, year, month, day, hour=hour)
```

It will return the results in dict type
```python
{'Z': array([-17.96973597]), 'Y': array([-17.27277951]), 'X': array([16.2994898])}
```
### Get time series of magnetic vectors 

```python
from DIFI import getSQfield
lat, lon = 20.5, 100.5
year = 2023
month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
day = 15
B = getSQfield(lat, lon, year, month, day)
```



