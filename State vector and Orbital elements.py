import numpy as np
import math as mth

MU = 398800.4418

X = input("Type the components of the state vector with a space:")  # state vector(r1, r2, r3, v1, v2, v3) in IJK
r = np.array(X.split(' ')[:3], dtype='float')  # coordinate vector
v = np.array(X.split(' ')[3:6], dtype='float')  # velocity vector

h = np.cross(r, v)  # angular momentum vector(defined for all orbits)
n = np.cross([0, 0, 1], h)  # node vector

# Eccentricity
e_vec = ((np.dot(v, v) - MU/np.linalg.norm(r)) * r - np.dot(r, v) * v )/ MU
e = np.linalg.norm(e_vec)

# Semimajor axis and semiparameter
if e == 1:  # case of parabolic orbits
    a = 'infinity'
    p = np.dot(h,h)/MU
else:
    a = (2/np.linalg.norm(r) - np.dot(v, v)/MU) ** (-1)
    p = a * (1 - e ** 2)

# Inclination
i = np.rad2deg(mth.acos(h[2]/np.linalg.norm(h)))

# Right ascension of the ascending node
if np.linalg.norm(n) == 0:
    if e_vec[1] < 0:
        RAAN = 'Undefined\ntrue longitude of periapsis(non-circular equatorial orbit),degrees = ' \
               + np.rad2deg(mth.acos(e_vec[0] / e))
    else:
        RAAN = 'Undefined\ntrue longitude of periapsis(non-circular equatorial orbit),degrees = ' \
               + 360.0 - np.rad2deg(mth.acos(e_vec[0] / e))
elif n[1] > 0:
    RAAN = np.rad2deg(mth.acos(n[0] / np.linalg.norm(n)))
else:
    RAAN = 360.0 - np.rad2deg(mth.acos(n[0] / np.linalg.norm(n)))

# Argument of perigee

if np.linalg.norm(n) == 0 or np.linalg.norm(e_vec) == 0:
    omega = 'Undefined'
elif e_vec[2] > 0:
    omega = np.rad2deg(mth.acos(np.dot(n,e_vec) / (e * np.linalg.norm(n))))
else:
    omega = 360.0 - np.rad2deg(mth.acos(np.dot(n,e_vec) / (e * np.linalg.norm(n))))

# True anomaly
if np.linalg.norm(e_vec) == 0:
    if r[2] > 0:
        t_a = 'Undefined\nargument of latitude(circular inclined orbit),degrees = ' \
              + np.rad2deg(mth.acos(np.dot(n,r) / (np.linalg.norm(n) * np.linalg.norm(r))))
    else:
        t_a = 'Undefined\nargument of latitude(circular inclined orbit),degrees = ' \
              + 360.0 - np.rad2deg(mth.acos(np.dot(n,r) / (np.linalg.norm(n) * np.linalg.norm(r))))
elif np.dot(r, v) > 0:
    t_a = np.rad2deg(mth.acos(np.dot(e_vec, r) / (e * np.linalg.norm(r))))
else:
    t_a = 360.0 - np.rad2deg(mth.acos(np.dot(e_vec, r)/(e * np.linalg.norm(r))))

print(' Semimajor axis, km = ', a, '\n',
      'Semiparameter , km = ', p, '\n',
      'Eccentricity = ', e, '\n',
      'Inclination,degrees = ', i ,'\n',
      'Right ascension of the ascending node,degrees = ', RAAN, '\n',
      'Argument of perigee,degrees  = ', omega, '\n',
      'True anomaly, degrees  = ',t_a)

if np.linalg.norm(n) == 0 and np.linalg.norm(e) == 0:
    if r[1] > 0:
        l = np.rad2deg(mth.acos(r[0] / np.linalg.norm(r)))
    else:
        l = 360.0 - np.rad2deg(mth.acos(r[0] / np.linalg.norm(r)))
    print('True longitude (circular equatorial orbit), degrees = ', l)

"""""
p = float(input("Semiparameter,km "))
e = float(input("Eccentricity "))
i = np.deg2rad(float(input("Inclination,degrees ")))
RAAN = np.deg2rad(float(input("Right ascension of the ascending node,degrees ")))
omega = np.deg2rad(float(input("Argument of perigee,degrees ")))
try:
    t_a = np.deg2rad(float(input("True anomaly, degrees ")))
except :
    try:
        l = np.deg2rad(float(input("True longitude, degrees ")))
    except :
        try:
            u = np.deg2rad(float(input("Argument of latitude, degrees ")))
        except:
            w_t = np.deg2rad(float(input("True longitude of periapsis, degrees ")))

if (omega, RAAN) == (0, 0):
    t_a = l
if omega == 0.0:
    t_a = u
if RAAN == 0.0:
    omega = w_t

# State vector in inertial system(PQW)

r_PQW = np.array([(p * mth.cos(t_a))/(1+e*mth.cos(t_a)),(p * mth.sin(t_a))/(1+e*mth.cos(t_a)),0],dtype='float')
v_PQW = np.array([-mth.sin(t_a)*(MU/p)**0.5,(e+mth.cos(t_a))*(MU/p)**0.5,0],dtype='float')
M = np.array([[mth.cos(RAAN)*mth.cos(omega) - mth.sin(RAAN)*mth.sin(omega)*mth.cos(i),
          -mth.cos(RAAN)*mth.sin(omega) - mth.sin(RAAN)*mth.cos(omega)*mth.cos(i),
          mth.sin(RAAN)*mth.sin(i)],
         [mth.sin(RAAN)*mth.cos(omega) + mth.cos(RAAN)*mth.sin(omega)*mth.cos(i),
          -mth.sin(RAAN)*mth.sin(omega) + mth.cos(RAAN)*mth.cos(omega)*mth.cos(i),
          -mth.cos(RAAN)*mth.sin(i)],
         [mth.sin(omega)*mth.sin(i),
          mth.cos(omega)*mth.sin(i),
          mth.cos(i)]],dtype='float')

print("Coordinate vector in IJK:", np.dot(M,r_PQW))
print("Velocity vector in IJK:", np.dot(M,v_PQW))
"""""