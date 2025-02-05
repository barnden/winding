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
        du = (np.array(self.f(u + self.eps, v)) - np.array(self.f(u - self.eps, v))) / (
            2.0 * self.eps
        )
        dv = (np.array(self.f(u, v + self.eps)) - np.array(self.f(u, v - self.eps))) / (
            2.0 * self.eps
        )

        normal = normalize(np.cross(du, dv, axis=0))

        return normal


class BSpline(Surface):
    v_min = 2
    basis = (
        1
        / 6
        * np.array(
            [
                [-1.0, 3.0, -3.0, 1.0],
                [3.0, -6.0, 3.0, 0.0],
                [-3.0, 0.0, 3.0, 0.0],
                [1.0, 4.0, 1.0, 0.0],
            ]
        )
    )

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

        if iv == 1:
            iv += 1

        u -= iu
        v -= iv

        iu %= int(self.u_max)

        return (u, v, iu, iv)

    def f(self, u, v):
        u, v, iu, iv = self.get_uv(u, v)

        U = np.array([u**3, u**2, u, 1.0]) @ self.basis
        V = np.array([v**3, v**2, v, 1.0]) @ self.basis

        iu_max = int(self.u_max)

        iu += iu_max

        P = np.array(
            [
                V
                @ np.array(
                    [
                        self.points[iv - 2][(iu - 2) % iu_max],
                        self.points[iv - 1][(iu - 2) % iu_max],
                        self.points[iv - 0][(iu - 2) % iu_max],
                        self.points[iv + 1][(iu - 2) % iu_max],
                    ]
                ),
                V
                @ np.array(
                    [
                        self.points[iv - 2][(iu - 1) % iu_max],
                        self.points[iv - 1][(iu - 1) % iu_max],
                        self.points[iv - 0][(iu - 1) % iu_max],
                        self.points[iv + 1][(iu - 1) % iu_max],
                    ]
                ),
                V
                @ np.array(
                    [
                        self.points[iv - 2][(iu + 0) % iu_max],
                        self.points[iv - 1][(iu + 0) % iu_max],
                        self.points[iv - 0][(iu + 0) % iu_max],
                        self.points[iv + 1][(iu + 0) % iu_max],
                    ]
                ),
                V
                @ np.array(
                    [
                        self.points[iv - 2][(iu + 1) % iu_max],
                        self.points[iv - 1][(iu + 1) % iu_max],
                        self.points[iv - 0][(iu + 1) % iu_max],
                        self.points[iv + 1][(iu + 1) % iu_max],
                    ]
                ),
            ]
        )

        return U @ P

    def generate_mesh(self, N=32):
        parametric_space = np.array(
            np.meshgrid(
                np.linspace(self.u_min, self.u_max, N),
                np.linspace(self.v_min, self.v_max, N),
            )
        )

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

        return {"curve": curve, "normals": normals}

    def generate_slice_u(self, u: float, N=32):
        parametric_space = np.linspace(self.v_min, self.v_max, N)

        curve = np.zeros((3, N))
        normals = np.zeros((3, N))

        for i in range(N):
            uv = (u, parametric_space[i])
            curve[:, i] = self.f(*uv)
            normals[:, i] = self.normal(*uv)

        return {"curve": curve, "normals": normals}


fig = mlab.figure(size=(800, 800))

surface = BSpline("../reference/vase-divot-deep.txt")
surface_mesh = surface.generate_mesh(512)


def plot():
    mapping = [
        ((0.203187, -0.306702, -0.118525), (15.9064, 8.54333), (15.9769, 8.555)),
        ((0.203179, -0.309672, -0.127305), (15.8285, 8.51423), (15.9753, 8.54262)),
        ((0.203171, -0.312642, -0.136085), (15.7431, 8.47762), (15.9737, 8.52998)),
        ((0.203162, -0.315612, -0.144864), (15.6476, 8.42836), (15.972, 8.51701)),
        ((0.203154, -0.318582, -0.153644), (15.5413, 8.35722), (15.9704, 8.50363)),
        ((0.203146, -0.321552, -0.162423), (15.4339, 8.25698), (15.9689, 8.48973)),
        ((0.203137, -0.324522, -0.171203), (15.3506, 8.14228), (15.2191, 7.77859)),
        ((0.203129, -0.327492, -0.179983), (15.3, 8.03471), (15.244, 7.71557)),
        ((0.203121, -0.330462, -0.188762), (15.2714, 7.93854), (15.2616, 7.66875)),
        ((0.203113, -0.333431, -0.197542), (15.2533, 7.84927), (15.2708, 7.63201)),
        ((0.203104, -0.336401, -0.206321), (15.2406, 7.76385), (15.2719, 7.60135)),
        ((0.203096, -0.339371, -0.215101), (15.2318, 7.68096), (15.2654, 7.57413)),
        ((0.203088, -0.342341, -0.22388), (15.2267, 7.6003), (15.2509, 7.54843)),
    ]
    mapped2 = [
        (15.9582, 8.48405),
        (15.9358, 8.40163),
        (15.9102, 8.3202),
        (15.8805, 8.24006),
        (15.8455, 8.16149),
        (15.8043, 8.08474),
        (15.7558, 8.01001),
        (15.6992, 7.93737),
        (15.6346, 7.86656),
        (15.5626, 7.79717),
        (15.4843, 7.72867),
        (15.4009, 7.66046),
        (15.3142, 7.5919),
    ]

    points = []
    mapped_points = []
    mapped2_points = []
    world_points = []
    nor = []

    n0 = np.array([0.365308, -0.365989, -0.855921])
    n1 = np.array([0.613349, 0.188238, 0.767052])
    d = np.array([-0.00011612, -0.0415787, -0.122914])
    d_hat = np.array([-0.000894908, -0.320437, -0.947269])
    eta = n0 + n1
    eta /= np.linalg.norm(eta)
    tau = np.cross(eta, d_hat)
    eps = np.cross(tau, d_hat)
    eps /= np.linalg.norm(eps)

    for world, original, mapped in mapping:
        points.append(surface.f(*original))
        mapped_points.append(surface.f(*mapped))
        world_points.append(world)
        nor.append(eps)

    for p in mapped2:
        mapped2_points.append(surface.f(*p))

    quiver = np.array(mapped_points) - np.array(world_points)
    quiverw = np.array(world_points) - np.array(points)
    normal = np.array(nor)

    mlab.quiver3d(
        [world_points[0][0]],
        [world_points[0][1]],
        [world_points[0][2]],
        *n0,
        color=(1.0, 1.0, 1.0),
        scale_factor=0.025
    )
    mlab.quiver3d(
        [world_points[-1][0]],
        [world_points[-1][1]],
        [world_points[-1][2]],
        *n1,
        color=(1.0, 1.0, 1.0),
        scale_factor=0.025
    )
    mlab.quiver3d(
        [world_points[0][0]],
        [world_points[0][1]],
        [world_points[0][2]],
        *eta,
        color=(1.0, 1.0, 0.0),
        scale_factor=0.025
    )
    mlab.quiver3d(
        [world_points[0][0]],
        [world_points[0][1]],
        [world_points[0][2]],
        *d_hat,
        color=(1.0, 0.0, 1.0),
        scale_factor=0.025
    )
    mlab.quiver3d(
        *np.array(world_points).T, *normal.T, color=(0.0, 1.0, 1.0), scale_factor=0.25
    )

    mlab.points3d(
        *np.array(mapped2_points).T, color=(1.0, 1.0, 0.0), scale_factor=0.005
    )
    mlab.points3d(*np.array(points).T, color=(1.0, 0.0, 0.0), scale_factor=0.005)
    mlab.points3d(*np.array(mapped_points).T, color=(0.0, 0.0, 1.0), scale_factor=0.005)
    mlab.points3d(*np.array(world_points).T, color=(0.0, 1.0, 0.0), scale_factor=0.005)
    # mlab.quiver3d(*np.array(world_points).T, *quiver.T, color=(1., 1., 0.), scale_factor=1., scale_mode="vector")
    # mlab.quiver3d(*np.array(points).T, *quiverw.T, color=(0., 1., 1.), scale_factor=1., scale_mode="vector")
    mlab.mesh(*surface_mesh, color=(0.6, 0.65, 0.8))

    mlab.show()


plot()
