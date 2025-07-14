import random
from DIFI import getSQfield

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
    main()
