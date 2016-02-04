import geomag
import numpy as np
import sys
import os

start_lat = -60
start_lon = -180
resolution = 10

lat_entries = (abs(start_lat) * 2) / resolution
lon_entries = (abs(start_lon) * 2) / resolution

lat = start_lat
lon = start_lon

i = 0
j = 0

dec = np.zeros(shape=(lat_entries+1, lon_entries+1), dtype='int')
inc = np.zeros(shape=(lat_entries+1, lon_entries+1), dtype='int')
ti  = np.zeros(shape=(lat_entries+1, lon_entries+1), dtype='int')

while i <= lat_entries:
    while j <= lon_entries:
        dec[i,j], inc[i,j], ti[i,j] = geomag.all(lat, lon)
        lon += resolution
        j += 1
    lon = start_lon
    lat += resolution
    j = 0
    i += 1

def print_matrix(matrix, rows, columns):
    i = 0
    j = 0
    while i < rows:
        row = "    { "
        while j < columns:
            row += "{:>3}".format(str(matrix[i,j]))
            if j < (columns - 1):
                row += ","
            j += 1
        row += "},"
        print(row)
        j = 0
        i += 1

# Print header
print("#ifndef _GEO_MAG_GENERATED_H_")
print("#define _GEO_MAG_GENERATED_H_")
print("")

# Print defines
print("#define SAMPLING_RES "     + str(resolution))
print("#define SAMPLING_MIN_LAT " + str(start_lat))
print("#define SAMPLING_MAX_LAT " + str(abs(start_lat)))
print("#define SAMPLING_MIN_LON " + str(start_lon))
print("#define SAMPLING_MAX_LON " + str(abs(start_lon)))
print("")

# Print tables
print("static const int8_t declination_table[13][37] = \\")
print("{")
print_matrix(dec, lat_entries, lon_entries)
print("};")

print("static const int8_t inclination_table[13][37] = \\")
print("{")
print_matrix(inc, lat_entries, lon_entries)
print("};")

print("static const int8_t magnitude_table[13][37] = \\")
print("{")
print_matrix(ti, lat_entries, lon_entries)
print("};")

# Closing definte
print("")
print("#endif /* _GEO_MAG_GENERATED_H_ */")

# Matrix printing
#np.set_printoptions(threshold='nan', linewidth='nan')
#print(dec)
#print(inc)
#print(ti)

