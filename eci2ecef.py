import math as mth
import re as reg
import numpy as np

#  ECI is J2000
#  ECEF is WGS84
#  ECI -> ECEF

w = 7.29211585275553e-5  # Average inertial rotation
arcs = 3600*180/mth.pi  # Arcseconds per radian
MJD_J2000 = 51544.5  # Modified Julian Date of J2000


def mod_jul_date(day, month, year, hour, minute, sec):
    return 367 * year - int(7 * (year + int((month+9)/12)) / 4) + int(275 * month/9) +\
           day + 1721013.5 + ((sec/60 + minute)/60 + hour)/24 - 2400000.5


def m_sid_time(day, month, year, hour, minute, sec):

    MJD = mod_jul_date(day, month, year, hour, minute, sec)
    secs = 86400  # Seconds per day
    MJD_0 = round(MJD,0)
    UTC = secs * (MJD - MJD_0)  # [seconds]
    T_0 = (MJD_0 - MJD_J2000)/36525
    T = (MJD - MJD_J2000)/36525
    gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UTC +\
           (0.093104 - 6.2e-6 * T) * T * T  # [seconds]
    return 2 * mth.pi * (gmst/secs - mth.floor(gmst/secs))  # [radians]


def precession_matrix(mjd_1, mjd_2):

    T = (mjd_1 - MJD_J2000)/36525
    dT = (mjd_2 - mjd_1)/36525

    # Precession angles
    zeta = -((2306.2181 + (1.39656 - 0.000139 * T) * T) +\
    ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / arcs

    z = -zeta - ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / arcs

    theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -\
    ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / arcs

    # Rotation matrices
    zeta_matrix = np.array([[mth.cos(zeta), mth.sin(zeta), 0.0],
                  [-mth.sin(zeta), mth.cos(zeta), 0.0],
                  [0.0, 0.0, 1.0]])
    z_matrix = np.array([[mth.cos(z), mth.sin(z), 0.0],
                [-mth.sin(z), mth.cos(z), 0.0],
                [0.0, 0.0, 1.0]])
    theta_matrix = np.array([[mth.cos(theta), 0.0, -mth.sin(theta)],
                    [0.0, 1.0, 0.0],
                    [mth.sin(theta), 0.0, mth.cos(theta)]])

    return (zeta_matrix.dot(z_matrix)).dot(theta_matrix)


def earth_rotation_matrix(s_time):
    return np.array([[mth.cos(s_time), mth.sin(s_time), 0.0],
                [-mth.sin(s_time), mth.cos(s_time), 0.0],
                [0.0, 0.0, 1.0]])  # Hour Angle matrix


def derivative_rotation(mjd, rot_matx):
    omega = w + 4.3e-15*((mjd - MJD_J2000)/36525)  # [rad/s]
    S = np.array([[0.0, 1.0, 0.0],
                  [-1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0]])
    d_rot_matx = omega * S.dot(rot_matx)
    return d_rot_matx


#  ECI vectors
r_eci = [float(item) for item in input("vector of coordinates, km (with space): ").split(' ')]
v_eci = [float(item) for item in input("velocity vector, km/s (with space): ").split(' ')]
epoch = [float(item) for item in reg.split(r'[,/:]', input("Epoch in UTC format "))]

mjd_utc = mod_jul_date(epoch[0], epoch[1], epoch[2], epoch[3], epoch[4], epoch[5])  # [days]
gmst = m_sid_time(epoch[0], epoch[1], epoch[2], epoch[3], epoch[4], epoch[5])  # [radians]

#  Matrices of transformation
P = precession_matrix(MJD_J2000, mjd_utc)
Theta = earth_rotation_matrix(mjd_utc)
dTheta = derivative_rotation(mjd_utc, Theta)
U = Theta.dot(P)
dU = dTheta.dot(P)

#  ECEF vectors
r_ecef = U.dot(r_eci)
v_ecef = U.dot(v_eci) + dU.dot(r_eci)

print("Coordinate ECEF vector: ")
print(r_ecef)
print("Velocity ECEF vector: ")
print(v_ecef)