import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class H1:
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 1

    z_min = 0
    z_max = 1

    @staticmethod
    def f(u, v):
        return (
            (v ** 2 - v + 0.4) * np.cos(u),
            (v ** 2 - v + 0.4) * np.sin(u),
            v
        )
    
    def generate_mesh(self):
        N = 32
        parametric_space = np.meshgrid(
            np.linspace(H1.u_min, H1.u_max, N),
            np.linspace(H1.v_min, H1.v_max, N)
        )

        return H1.f(*parametric_space)
    
    def winding_path(self):
        with open("hyperboloid.txt") as f:
            points = [*map(lambda line: H1.f(*map(float, line.replace(",", "").split()[-2:])), f.readlines())]

        return points


class Vase:
    u_min = 0
    u_max = 2. * np.pi
    v_min = -1
    v_max = 1

    z_min = -1
    z_max = 1

    @staticmethod
    def f(u, v):
        return (
            (0.5 - 0.4 * v * (1. + v) * (1. - v)) * np.cos(u),
            (0.5 + 0.4 * v * (1. + v) * (1. - v)) * np.sin(u),
            v
        )
    
    def generate_mesh(self):
        N = 32
        parametric_space = np.meshgrid(
            np.linspace(Vase.u_min, Vase.u_max, N),
            np.linspace(Vase.v_min, Vase.v_max, N)
        )

        return Vase.f(*parametric_space)
    
    def winding_path(self):
        with open("vase.txt") as f:
            points = [*map(lambda line: Vase.f(*map(float, line.replace(",", "").split()[-2:])), f.readlines())]

        return points
    
class Spring:
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 4. * np.pi

    r1 = 1
    r2 = 0.3
    kh = 0.15

    z_min = 0
    z_max = 2

    @staticmethod
    def f(u, v):
        return (
            (Spring.r1 + Spring.r2 * np.cos(u)) * np.cos(v),
            (Spring.r1 + Spring.r2 * np.cos(u)) * np.sin(v),
            -Spring.r2 * np.sin(u) + Spring.kh * v
        )
    
    def generate_mesh(self):
        N = 32
        parametric_space = np.meshgrid(
            np.linspace(Spring.u_min, Spring.u_max, N),
            np.linspace(Spring.v_min, Spring.v_max, N)
        )

        return Spring.f(*parametric_space)
    
    def winding_path(self):
        with open("spring.txt") as f:
            points = [*map(lambda line: Spring.f(*map(float, line.replace(",", "").split()[-2:])), f.readlines())]

        return points

class Torus:
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 2. * np.pi - 0.001
    z_min = -1.
    z_max = 1.

    r1 = 1
    r2 = 0.3

    @staticmethod
    def f(u, v):
        return (
            (Torus.r1 + Torus.r2 * np.cos(u)) * np.cos(v),
            (Torus.r1 + Torus.r2 * np.cos(u)) * np.sin(v),
            -Torus.r2 * np.sin(u)
        )

    def generate_mesh(self):
        N = 32
        parametric_space = np.meshgrid(
            np.linspace(Torus.u_min, Torus.u_max, N),
            np.linspace(Torus.v_min, Torus.v_max, N)
        )

        return Torus.f(*parametric_space)
    

    def winding_path(self):
        with open("torus.txt") as f:
            points = [*map(lambda line: Torus.f(*map(float, line.replace(",", "").split()[-2:])), f.readlines())]

        return points

class TrefoilKnot:
    u_min = 0
    u_max = 2. * np.pi
    v_min = 0
    v_max = 2. * np.pi * 0.001
    z_min = 0.
    z_max = 2.

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
        e1 = e1[..., None].repeat(32, 1)[..., None].repeat(32, 2)
        e1 = (e1 - T) * t
        e1 /= np.linalg.norm(e1, axis=0)

        e2 = np.cross(t, e1, axis=0)

        return p0 + 0.4 * np.cos(u) * e1 + 0.4 * np.sin(u) * e2
    
       
    def generate_mesh(self):
        N = 32
        parametric_space = np.meshgrid(
            np.linspace(TrefoilKnot.u_min, TrefoilKnot.u_max, N),
            np.linspace(TrefoilKnot.v_min, TrefoilKnot.v_max, N)
        )

        return TrefoilKnot.f(*parametric_space)


for surface in [H1(), Vase(), Torus(), Spring()]:
    X = surface.generate_mesh()
    points = surface.winding_path()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(surface.z_min, surface.z_max)
    ax.plot_surface(*X, linewidth=0.1, zorder=1)
    ax.plot([e[0] for e in points], [e[1] for e in points], [e[2] for e in points], markersize=10, color=[1, 0, 0], zorder=4)
    plt.show()