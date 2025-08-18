import random
from DIFI import getSQfield
import matplotlib.pyplot as plt
import numpy as np

def list_inputs():
    year = 2023
    month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    day = 15

    N = len(month)
    lat = [random.uniform(-90, 90) for _ in range(N)]
    lon = [random.uniform(-180, 180) for _ in range(N)]

    B = getSQfield(lat, lon, year, month, day, model_name="xdifi2")

    print(B)

def point_inputs():
    lat, lon, year, month, day, hour = 20.5, 100.5, 2024, 6, 6, 0
    B = getSQfield(lat, lon, year, month, day, hour=hour)

    print(B)

def main():
    print("Example 1: List of inputs")
    list_inputs()

    print("\nExample 2: Point inputs")
    point_inputs()

if __name__ == "__main__":
    # def getSQfield(lat: Union[float, list], lon: Union[float, list], year: Union[int, list], 
    #                month: Union[int, list], day: Union[int, list], hour: Union[int, list]=0, 
    #                minutes: Union[int, list]=0, h: Union[float, list]=0,r: Union[float, list]=0, 
    #                f107_1: Optional[Union[float, list]]=None, model_name: Optional[Union[str]]="xdifi2", 
    #                geoc:Optional[bool] = False, return_geoc:Optional[bool] = False) -> dict:

    N_datapoints = 1000
    lat = 0*np.ones(N_datapoints)
    lon = 20*np.ones(N_datapoints)
    alt = np.ones(N_datapoints)
    year = np.int64(np.linspace(2023,2023,N_datapoints))
    month = np.int64(np.linspace(4,4,N_datapoints))
    day = np.int64(np.linspace(6,13,N_datapoints))
    y = 12*np.sin(np.linspace(0,14*np.pi, N_datapoints)) + 12
    hour = np.int64(24*np.sin(day*np.pi))
    plt.plot(np.linspace(0,2*np.pi, N_datapoints),y, label = 'artificial')
    plt.plot(day,hour, label = "supposed to be clock")
    plt.legend()
    plt.show()
    # hour = np.int64(np.linspace(0,23,N_datapoints))

    minute = np.int64(np.linspace(0,0,N_datapoints))
    second = np.int64(np.linspace(0,0,N_datapoints))
    out_b = getSQfield(lat, lon, year, month, day, hour = hour)
    plt.plot(day, out_b['Z'])
    plt.show()
    plt.plot(day, out_b['X'])
    plt.show()
    plt.plot(day, out_b['Y'])
    plt.show()
    # main()
