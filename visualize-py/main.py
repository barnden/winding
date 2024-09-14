import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Surface:
    x_limits = [-1, 1]
    y_limits = [-1, 1]
    z_limits = [-1, 1]

    winding_path_file = "points.txt"

    def generate_mesh(self, N=32):
        parametric_space = np.meshgrid(
            np.linspace(self.u_min, self.u_max, N),
            np.linspace(self.v_min, self.v_max, N),
        )

        return self.f(*parametric_space)

    def winding_path(self):
        with open(self.winding_path_file) as f:
            points = [*map(lambda line: self.f(*map(lambda x: np.array([[float(x)]]), line.replace(",", "").split()[-2:])), f.readlines())]

        return points

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

for surface in [H1(), Vase(), Torus(), Spring(), TrefoilKnot()]: #[H1(), Vase(), Torus(), Spring()]:
    X = surface.generate_mesh(32)
    points = surface.winding_path()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlim3d(*surface.x_limits)
    ax.set_ylim3d(*surface.y_limits)
    ax.set_zlim3d(*surface.z_limits)
    ax.plot_surface(*X, linewidth=0.1, zorder=1)
    ax.plot([e[0] for e in points], [e[1] for e in points], [e[2] for e in points], markersize=10, color=[1, 0, 0], zorder=4)
    plt.show()
