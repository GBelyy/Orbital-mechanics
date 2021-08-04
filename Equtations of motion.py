import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp


MU = 398800.4418
J2 = 1.082e-3
R_e = 6378.1363  # [km]

#(r[0]**2 + r[1]**2 + r[2]**2)**0.5
def rhs(s,r):
    return \
    [r[3],
     r[4],
     r[5],
     -MU * r[0] / (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 1.5 + \
     3 * MU * J2 * R_e ** 2 / (2 * (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 2.5) \
     * (5/(r[0]**2 + r[1]**2 + r[2]**2)-1)*r[0],
     -MU * r[1] / (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 1.5 + \
     3 * MU * J2 * R_e ** 2 / (2 * (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 2.5) \
     * (5/(r[0]**2 + r[1]**2 + r[2]**2)-1)*r[1],
     -MU * r[2] / (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 1.5 + \
     3 * MU * J2 * R_e ** 2 / (2 * (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 2.5) \
     * ((5 / (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) - 1) * r[2]-2)]
res = solve_ivp(rhs, (0, 1000), [0, 10, 10000, 0, 0, 0])

print(res.message)


fig = plt.figure()
ax = Axes3D(fig)
ax.plot(res.y[0],res.y[1],res.y[2], label ='parametric curve')
ax.set_xlabel('X, km')
ax.set_ylabel('Y, km')
ax.set_zlabel('Z, km')
plt.show()
