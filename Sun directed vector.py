import re as reg
import math as mth
import numpy as np

def Jul_date(day, month, year, hour, min, sec):
    return 367 * year - int(7 * (year + int((month+9)/12)) / 4) + int(275 * month/9)\
           + day + 1721013.5 + ((sec/60 + min)/60 + hour)/24

AU = 149597870.7  # astronomical unit [km]
R_e = 6378.1363  # Earth radius [km]

# Paris location
#  longitude = 2.35
#  geodetic latitude = 48.86
#  time zone = UTC+3
#  16/08/2021,12:00:00
#  average altitude in Paris = 0.03 km

longitude = mth.radians(float(input("Local longitude of POI,degrees ")))
latitude = mth.radians(float(input("Local(geodetic) latitude of POI,degrees ")))
altitude = float(input("Average altitude of POI,km "))
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

#  module of geocentric sun directed vector
r = 1.000140612 - 0.016708617*mth.cos(mth.radians(M))\
    - 0.000139589*mth.cos(mth.radians(2*M))

# Geocentric position vector of the Sun, [km]
r_vec = np.array([r * AU * mth.cos(mth.radians(l_e)),
                  r * AU * mth.cos(mth.radians(o_e)) * mth.sin(mth.radians(l_e)),
                  r * AU * mth.sin(mth.radians(o_e)) * mth.sin(mth.radians(l_e))])

# eccentricity of th Earth
e = 0.01670862 - 0.000042037*T - .0000001236*(T**2) + 0.00000000004*(T**3)
# radius of curvature in the meridian, [km]
C_0 = R_e/mth.sqrt(1 - (e**2) * (mth.sin(latitude))**2)
# Coordinates of POI in geocentric frame, [km]
poi_vec = np.array([(C_0 + altitude) * mth.cos(latitude) * mth.cos(longitude),
                    (C_0 + altitude) * mth.cos(latitude) * mth.sin(longitude),
                    (C_0 * (1 - e**2) + altitude) * mth.sin(latitude)])

# Sun directed vector from POI
res_vec = r_vec - poi_vec

print("Vector of sun direction is (in km)")
print([round(i,5) for i in res_vec])

# Unit vector of sun direction
print("Unit vector of sun direction ")
print([round(i,5) for i in res_vec/np.linalg.norm(res_vec)])


# Sun elevation angle
C = 360  # 2pi radians in degrees
td = 24  # length of a day in hours
tilt = 23.5  # tilt of Earth on its axes
d_year = 365  # length of a year in days
d_sol = 173  # Julian day of Summer Solstice

days = JD - Jul_date(1,1,epoch[2],epoch[3],epoch[4],epoch[5]) + 1  # Julian day of a current day from the beginning of the year
time = epoch[3] + (epoch[4] + epoch[5]/60)/60  # time in UTC
decl = mth.radians(tilt * mth.cos(mth.radians(C*(days - d_sol)/d_year))) # solar declination angle
theta = mth.asin(mth.sin(latitude)*mth.sin(decl) - mth.cos(latitude)\
                 *mth.cos(decl)*mth.cos(mth.radians(C*time/td) - longitude)) * 180/mth.pi

print('Sun elevation angle is ', round(theta, 5), 'degrees')

""""
#???????????? ???? ????????????

longitude = radians(2.35) # float(input("Local longitude of POI,degrees "))
latitude = radians(48.86) # float(input("Local latitude of POI,degrees "))
epoch = [float(item) for item in reg.split(r'[,/:]', '16/08/2021,00:00:00')]
C = 360  # 2pi radians in degrees
td = 24  # length of a day in hours
tilt = 23.5  # tilt of Earth on its axes
d_year = 365  # length of a year in days
d_sol = 173  # Julian day of Summer Solstice

JD = Jul_date(epoch[0],epoch[1],epoch[2],epoch[3],epoch[4],epoch[5]) - Jul_date(1,1,epoch[2],epoch[3],epoch[4],epoch[5]) + 1
t = []
angle = []
for i in range(0,25):
    time = i + (epoch[4] + epoch[5]/60)/60
    decl = radians(tilt * cos(radians(C*(JD - 173)/365)))
    theta = asin(sin(latitude)*sin(decl) - cos(latitude)*cos(decl)*cos(radians(C*time/24) - longitude)) * 180/pi
    angle.append(theta)
    t.append(i)
fig, ax = plt.subplots()
ax.plot(t, angle)
ax.set_xlabel('Time at 16 Avg, hours')
ax.set_ylabel('Sun elevation angle, degrees')
ax.set_xlim([min(t), max(t)])
ax.set_ylim([round(min(angle)-5,-1), floor(max(angle)+5)])
ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))

ax.grid(which='major',
        color = 'k')

ax.grid(which='minor',
        color = 'gray',
        linestyle = ':')
        
ax.minorticks_on()
fig.set_figwidth(8)
fig.set_figheight(8)

plt.show()
"""
