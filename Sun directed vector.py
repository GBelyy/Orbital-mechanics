import re as reg
import math as mth
import numpy as np

def Jul_date(day, month, year, hour, min, sec):
    return 367 * year - int(7 * (year + int((month+9)/12)) / 4) + int(275 * month/9)\
           + day + 1721013.5 + ((sec/60 + min)/60 + hour)/24

AU = 149597870.7  # астрономическая единица

# Paris location
#  longitude = 2.35
#  latitude = 48.86
#  time zone = UTC+3
#  16/08/2021,12:00:00

longitude = float(input("Local longitude of POI,degrees "))
latitude = float(input("Local latitude of POI,degrees "))
zone = float(input("Type time zone in UTC: "))
epoch = [float(item) for item in reg.split(r'[,/:]', input("Epoch in UTC format "))]

#  epoch -> julian date
JD = Jul_date(epoch[0],epoch[1],epoch[2],epoch[3],epoch[4],epoch[5])

#  julian centuries
T = (JD - 2451545.0)/36525

#  obliquity of the ecliptic
o_e = 23.439291 - 0.013004*T

#  mean longitude of the Sun
l_M = 280.46 + 36000.771 * T

#  mean anomaly for the Sun
M = 357.5291092 + 35999.05034 * T

#  ecliptic longitude
l_e = l_M + 1.914666471* mth.sin(mth.radians(M))\
      + 0.019994643*mth.sin(mth.radians(2*M))

#  ecliptic altitude
phi_e = 0.000333

#  module of sun directed vector
r = 1.000140612 - 0.016708617*mth.cos(mth.radians(M))\
    - 0.000139589*mth.cos(mth.radians(2*M))

r_vec = [r * mth.cos(mth.radians(l_e)),
         r * mth.cos(mth.radians(o_e)) * mth.sin(mth.radians(l_e)),
         r * mth.sin(mth.radians(o_e)) * mth.sin(mth.radians(l_e))]
print("Vector of sun direction is (in km)")
print([round(i * AU,5) for i in r_vec])

# Unit vector of sun direction
print("Unit vector of sun direction is (in AU)")
print([round(i,5) for i in r_vec/np.linalg.norm(r_vec)])

# Sun elevation angle
# Local Standard Time Meridian (LSTM)
LSTM = 15 * zone

# Equation of Time (EoT)
B = 360/365 * (JD - Jul_date(1,1,epoch[2],0,0,0) - 81)
EoT = 9.87 * mth.sin(mth.radians(2*B)) - 7.53 * mth.cos(mth.radians(B)) - 1.5 * mth.sin(mth.radians(B))

# Time Correction Factor (TC),minutes
TC = 4*(longitude - LSTM) + EoT

# Local Solar Time (LST),hours
LST = epoch[3] + (epoch[4] + epoch[5]/60)/60 + TC/60

# Hour Angle (HRA),degrees
HRA = 15 * (LST - 12)

# Declination angle
decl = -23.45 * mth.cos(360/365 * (JD - Jul_date(1,1,epoch[2],0,0,0) + 10))

#Elevation angle
theta = mth.degrees(mth.asin(mth.sin(mth.radians(decl)) * mth.sin(mth.radians(latitude)) \
        + mth.cos(mth.radians(decl)) * mth.cos(mth.radians(latitude)) * mth.cos(mth.radians(HRA))))
print('Sun elevation angle is ', round(theta, 5), 'degrees')



