import numpy as np
import re
from mayavi import mlab
from scipy.spatial import ConvexHull
from itertools import groupby

def normalize(v: np.ndarray):
    norm = np.linalg.norm(v, axis=0)
    
    if norm < 1e-6:
        return np.zeros_like(v)
    
    return v / norm

class Surface:
    u_min = 0
    v_min = 0
    eps = 1e-5

    def __init__(self):
        pass

    def f(self, u, v):
        return np.zeros_like(u)

    def normal(self, u: np.ndarray, v: np.ndarray):
        du = (np.array(self.f(u + self.eps, v)) - np.array(self.f(u - self.eps, v))) / (2. * self.eps)
        dv = (np.array(self.f(u, v + self.eps)) - np.array(self.f(u, v - self.eps))) / (2. * self.eps)

        normal = normalize(np.cross(du, dv, axis=0))

        return normal
    
class BSpline(Surface):
    v_min = 2
    basis = 1 / 6 * np.array([[-1., 3., -3., 1.], [3., -6., 3., 0.], [-3., 0., 3., 0.], [1., 4., 1., 0.]])

    def __init__(self, path_to_cpts):
        with open(path_to_cpts) as f:
            lines = f.readlines()

        self.u_max, self.v_max = list(map(float, lines[0].split()))
        self.v_max -= 1

        self.points = np.zeros((int(self.v_max) + 1, int(self.u_max), 3))

        lines = list(map(lambda x: list(map(float, x.split())), lines[1:]))
        line_num = 0
        for v in range(int(self.v_max) + 1):
            for u in range(int(self.u_max)):
                self.points[v][u] = lines[line_num]
                line_num += 1

    def get_uv(self, u, v):
        while u < 0:
            u += self.u_max

        iu = int(u)
        iv = int(v)

        if iv == int(self.v_max):
            iv -= 1

        if (iv == 1):
            iv += 1

        u -= iu
        v -= iv

        iu %= int(self.u_max)

        return (u, v, iu, iv)
    
    def f(self, u, v):
        u, v, iu, iv = self.get_uv(u, v)

        U = np.array([u ** 3, u ** 2, u, 1.]) @ self.basis
        V = np.array([v ** 3, v ** 2, v, 1.]) @ self.basis

        iu_max = int(self.u_max)

        iu += iu_max

        P = np.array([
            V @ np.array([self.points[iv - 2][(iu - 2) % iu_max], self.points[iv - 1][(iu - 2) % iu_max], self.points[iv - 0][(iu - 2) % iu_max], self.points[iv + 1][(iu - 2) % iu_max]]),
            V @ np.array([self.points[iv - 2][(iu - 1) % iu_max], self.points[iv - 1][(iu - 1) % iu_max], self.points[iv - 0][(iu - 1) % iu_max], self.points[iv + 1][(iu - 1) % iu_max]]),
            V @ np.array([self.points[iv - 2][(iu + 0) % iu_max], self.points[iv - 1][(iu + 0) % iu_max], self.points[iv - 0][(iu + 0) % iu_max], self.points[iv + 1][(iu + 0) % iu_max]]),
            V @ np.array([self.points[iv - 2][(iu + 1) % iu_max], self.points[iv - 1][(iu + 1) % iu_max], self.points[iv - 0][(iu + 1) % iu_max], self.points[iv + 1][(iu + 1) % iu_max]]),
        ])

        return U @ P
    
    def generate_mesh(self, N=32):
        parametric_space = np.array(np.meshgrid(
            np.linspace(self.u_min, self.u_max, N),
            np.linspace(self.v_min, self.v_max, N),
        ))

        mesh = np.zeros((3, N, N))

        for i in range(N):
            for j in range(N):
                mesh[:, i, j] = self.f(*parametric_space[:, i, j])

        return mesh
    
    def generate_slice(self, v: float, N=32):
        parametric_space = np.linspace(self.u_min, self.u_max, N)

        curve = np.zeros((3, N))
        normals = np.zeros((3, N))

        for i in range(N):
            uv = (parametric_space[i], v)
            curve[:, i] = self.f(*uv)
            normals[:, i] = self.normal(*uv)

        return {
            "curve": curve,
            "normals": normals
        }
    
    def generate_slice_u(self, u: float, N=32):
        parametric_space = np.linspace(self.v_min, self.v_max, N)

        curve = np.zeros((3, N))
        normals = np.zeros((3, N))

        for i in range(N):
            uv = (u, parametric_space[i])
            curve[:, i] = self.f(*uv)
            normals[:, i] = self.normal(*uv)

        return {
            "curve": curve,
            "normals": normals
        }

fig = mlab.figure(size=(800,800))

surface = BSpline("../reference/vase-divot-deep.txt")
surface_mesh = surface.generate_mesh(512)

def plot():
    mapping = [
        ((15.7259, 8.30581),(15.8041, 8.34437)),
        ((15.625, 8.24797), (15.7815, 8.34445)),
        ((15.5232, 8.1641), (15.4201, 8.0958)),
        ((15.4425, 8.06735), (15.3638, 7.92585)),
        ((15.3893, 7.97641), (15.3497, 7.83903)),
        ((15.3529, 7.89411), (15.3393, 7.772)),
        ((15.3238, 7.8168), (15.3276, 7.71506)),
        ((15.298, 7.74207), (15.3124, 7.6641)),
        ((15.2737, 7.66863), (15.292, 7.6165)),
        ((15.25, 7.59566), (15.2645, 7.57011)),
    ]
    points = []
    mapped_points = []

    for original, mapped in mapping:
        points.append(surface.f(*original))
        mapped_points.append(surface.f(*mapped))

    mlab.points3d(*np.array(points).T, color=(1., 0., 0.), scale_factor=0.005)
    mlab.points3d(*np.array(mapped_points).T, color=(0., 0., 1.), scale_factor=0.005)
    mlab.mesh(*surface_mesh, color=(0.6, 0.65, 0.8))

    mlab.show()
plot()