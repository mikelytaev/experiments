from scipy.integrate import RK45
import numpy as np
import math as fm

import matplotlib.pyplot as plt


def solve(theta_degrees, v0, epsilon):
    a_matrix = np.diag([1, 1], k=2)
    grav_acc = 10
    b_vector = np.array([0, 0, 0, -grav_acc])


    def f(x):
        v1 = x[2]
        v2 = x[3]
        return np.array([0, 0, -fm.sqrt(v1*v1 + v2*v2)*v1, -fm.sqrt(v1*v1 + v2*v2)*v2])


    def rhs(t, x):
        return a_matrix.dot(x) + b_vector + epsilon * f(x)


    v0x = fm.cos(theta_degrees * fm.pi / 180) * 500
    v0y = fm.sin(theta_degrees * fm.pi / 180) * 500
    rk =RK45(fun=rhs, t0=0, y0=np.array([0, 0, v0x, v0y]), t_bound=5000, vectorized=False, max_step=0.1)

    coords = []
    for i in range(10000):
        rk.step()
        #print(rk.t)
        #print(rk.y)
        coords += [(rk.y[0], rk.y[1])]
        if rk.y[1] < -0.0:
            break

    x1, y1, x2, y2 = coords[-2][0], coords[-2][1], coords[-1][0], coords[-1][1]
    x = -(x2-x1)/(y2-y1)*y1 + x1
    return coords, x


epsilon = 0.0001
v0 = 500

thetas = []
rs = []
all_coords = []
for theta_degrees in np.arange(10, 80, 0.05):
    coords, x = solve(theta_degrees, v0, epsilon)
    all_coords += coords
    print("theta=" + str(theta_degrees) + ", x="+ str(x))
    thetas += [theta_degrees]
    rs += [x]

plt.scatter(thetas, rs, s=1)
plt.grid(True)
plt.xlabel("Угол, градусы")
plt.ylabel("Место падения, м")
plt.show()

plt.scatter([c[0] for c in all_coords], [c[1] for c in all_coords], s=0.1)
plt.grid(True)
plt.xlim([0, 10000])
plt.xlabel("Расстояние, м")
plt.ylabel("Высота, м")
plt.show()

# plt.plot([c[0] for c in coords], [c[1] for c in coords])
# plt.grid(True)
# plt.show()