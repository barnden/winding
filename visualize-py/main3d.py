import numpy as np
import re
from mayavi import mlab
from matplotlib import pyplot as plt

class Surface:
    x_limits = [-1, 1]
    y_limits = [-1, 1]
    z_limits = [-1, 1]
    eps = 1e-5

    winding_path_file = "points.txt"
    def __init__(self, file=None):
        if file is not None:
            self.winding_path_file = file

    def normal(self, u, v):
        du = (np.array(self.f(u + self.eps, v)) - np.array(self.f(u - self.eps, v))) / (2. * self.eps)
        dv = (np.array(self.f(u, v + self.eps)) - np.array(self.f(u, v - self.eps))) / (2. * self.eps)

        normal = np.cross(du, dv, axis=0)
        normal /= np.linalg.norm(normal, axis=0)

        return normal

    def generate_mesh(self, N=32):
        parametric_space = np.meshgrid(
            np.linspace(self.u_min, self.u_max, N),
            np.linspace(self.v_min, self.v_max, N),
        )

        return self.f(*parametric_space)

    def winding_path(self):
        with open(self.winding_path_file) as f:
            data = list(map(lambda line: list(map(lambda x: float(x), re.sub(',', '', line).split())),f.readlines()))

        return data

class H1(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 1

    z_limits = [0, 1]

    winding_path_file = "hyperboloid.txt"

    @staticmethod
    def f(u, v):
        return (
            (v ** 2 - v + 0.4) * np.cos(u),
            (v ** 2 - v + 0.4) * np.sin(u),
            v
        )

class Vase(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = -1
    v_max = 1

    z_limits = [-1, 1]
    winding_path_file = "vase.txt"

    @staticmethod
    def f(u, v):
        return (
            (0.5 - 0.4 * v * (1. + v) * (1. - v)) * np.cos(u),
            (0.5 + 0.4 * v * (1. + v) * (1. - v)) * np.sin(u),
            v
        )

class Spring(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 4. * np.pi

    r1 = 1
    r2 = 0.3
    kh = 0.15

    z_limits = [0, 2]
    winding_path_file = "spring.txt"

    @staticmethod
    def f(u, v):
        return (
            (Spring.r1 + Spring.r2 * np.cos(u)) * np.cos(v),
            (Spring.r1 + Spring.r2 * np.cos(u)) * np.sin(v),
            -Spring.r2 * np.sin(u) + Spring.kh * v
        )
class Torus(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 2. * np.pi - 0.001

    r1 = 1
    r2 = 0.3

    z_limits = [-1, 1]
    winding_path_file = "torus.txt"

    @staticmethod
    def f(u, v):
        return (
            (Torus.r1 + Torus.r2 * np.cos(u)) * np.cos(v),
            (Torus.r1 + Torus.r2 * np.cos(u)) * np.sin(v),
            -Torus.r2 * np.sin(u)
        )

class TrefoilKnot(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 2. * np.pi - 0.001

    x_limits = [-3, 3]
    y_limits = [-3, 3]
    z_limits = [-3, 3]
    winding_path_file = "trefoil.txt"

    @staticmethod
    def f(u, v):
        p0 = np.array([
            np.sin(v) + 2. * np.sin(2. * v),
            np.cos(v) - 2. * np.cos(2. * v),
            -1. * np.sin(3. * v)
        ])

        t = np.array([
            np.cos(v) + 4. * np.cos(2. * v),
            -np.sin(v) + 4. * np.sin(2. * v),
            -3. * np.cos(3. * v)
        ])
        t /= np.linalg.norm(t, axis=0)

        e1 = np.array([0., 0., 1.])
        T = np.tensordot(e1, t, axes=[[0], [0]])

        e1 = e1[..., None, None].repeat(u.shape[0], 1).repeat(u.shape[1], 2)
        e1 = e1 - (T * t)
        e1 /= np.linalg.norm(e1, axis=0)

        e2 = np.cross(t, e1, axis=0)

        return p0 + 0.4 * np.cos(u) * e1 + 0.4 * np.sin(u) * e2

class BSpline(Surface):
    u_min = 0
    u_max = 2. * np.pi
    v_min = 2
    v_max = 2. * np.pi - 0.001

    x_limits = [-3, 3]
    y_limits = [-3, 3]
    z_limits = [-3, 3]
    winding_path_file = "trefoil.txt"
    control_point_path_file = "spline.txt"

    basis = 1 / 6 * np.array([[-1., 3., -3., 1.], [3., -6., 3., 0.], [-3., 0., 3., 0.], [1., 4., 1., 0.]])

    def __init__(self, winding_path, control_path):
        self.winding_path_file = winding_path
        self.control_point_path_file = control_path

        with open(self.control_point_path_file) as f:
            lines = f.readlines()

        self.u_max, self.v_max = list(map(float, lines[0].split()))
        self.v_max -= 1

        self.points = np.zeros((int(self.v_max) + 1, 3 * int(self.u_max), 3))

        lines = list(map(lambda x: list(map(float, x.split())), lines[1:]))
        k = 0
        for i in range(int(self.v_max) + 1):
            for j in range(int(self.u_max)):
                self.points[i][j] = lines[k]

                k += 1

    def get_uv(self, u, v):
        while u < 0.:
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

        U = np.array([u ** 3, u ** 2, u, 1.]) @ self.basis
        V = np.array([v ** 3, v ** 2, v, 1.]) @ self.basis

        iu += int(self.u_max)
        P = np.array([
            V @ np.array([self.points[iv - 2][(iu - 2) % int(self.u_max)], self.points[iv - 1][(iu - 2) % int(self.u_max)], self.points[iv - 0][(iu - 2) % int(self.u_max)], self.points[iv + 1][(iu - 2) % int(self.u_max)]]),
            V @ np.array([self.points[iv - 2][(iu - 1) % int(self.u_max)], self.points[iv - 1][(iu - 1) % int(self.u_max)], self.points[iv - 0][(iu - 1) % int(self.u_max)], self.points[iv + 1][(iu - 1) % int(self.u_max)]]),
            V @ np.array([self.points[iv - 2][(iu + 0) % int(self.u_max)], self.points[iv - 1][(iu + 0) % int(self.u_max)], self.points[iv - 0][(iu + 0) % int(self.u_max)], self.points[iv + 1][(iu + 0) % int(self.u_max)]]),
            V @ np.array([self.points[iv - 2][(iu + 1) % int(self.u_max)], self.points[iv - 1][(iu + 1) % int(self.u_max)], self.points[iv - 0][(iu + 1) % int(self.u_max)], self.points[iv + 1][(iu + 1) % int(self.u_max)]]),
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


if __name__ == "__main__":
    cmap = plt.get_cmap('viridis')
    cmaplist = np.array([cmap(i) for i in range(cmap.N)]) * 255

    # root = "surface-editing"
    # surfaces = [BSpline(f"{root}/path-{x}.txt", f"{root}/spline-step-{x}.txt") for x in range(12)]

    surfaces = [
        # BSpline("../build/path-0.txt", "surface-editing/spline-step-0.txt"),
        BSpline(None, "../build/spline-step-1.txt")
    ]

    for surface in surfaces:
        X = surface.generate_mesh(128)

        if surface.winding_path_file is not None:
            path = np.array(surface.winding_path())

        # paths = np.array(simulate())

        for i, title in [(4, "Surface Dist."), (5, "Path Dist."), (6, "Angle")]:
            mlab.figure(title, size=(1080, 720), fgcolor=(1., 1., 1.), bgcolor=(.1, .1, .1))

            if isinstance(surface, BSpline):
                J = surface.points.shape[1] // 3
                mlab.points3d(surface.points[:, :J, 0], surface.points[:, :J, 1], surface.points[:, :J, 2], scale_factor=.025)

                for j in range(surface.points.shape[0]):
                    mlab.plot3d(surface.points[j, :J, 0], surface.points[j, :J, 1], surface.points[j, :J, 2], tube_radius=0.003)

                for j in range(J):
                    mlab.plot3d(surface.points[:, j, 0], surface.points[:, j, 1], surface.points[:, j, 2], tube_radius=0.003)

            mlab.mesh(*X, color=(0.6, 0.65, 0.8))

            if surface.winding_path_file is not None:
                # 4: dist to surface, 5: path dist, 6: angles
                winding_plot = mlab.plot3d(*np.array(path)[:, 1:4].T, abs(np.array(path)[:, i].T), tube_radius=0.003)
                winding_plot.module_manager.scalar_lut_manager.lut.table = cmaplist

            view = (90., 60, 5, (0, 0, 0))
            mlab.view(*view, reset_roll=True)
            mlab.show()

            break
