import numpy as np
import re
from mayavi import mlab

class Surface:
    x_limits = [-1, 1]
    y_limits = [-1, 1]
    z_limits = [-1, 1]
    eps = 1e-5

    winding_path_file = "points.txt"

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
            data = map(lambda x: re.sub('[,()]', '', x).split()[1:], f.readlines())
            data = [[int(x[0]), np.array([[float(x[1])]]), np.array([[float(x[2])]])] for x in data]

        points = []
        i = 0
        while i < len(data):
            r, *uv = data[i]
            points.append(self.f(*uv) + 0.001 * self.normal(*uv))
            i += r

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
    
class Vase2(Vase):
    winding_path_file = "vase-2.txt"

if __name__ == "__main__":
    for surface in [Vase2()]:
        X = surface.generate_mesh(128)
        path = surface.winding_path()

        mlab.figure(size=(1080, 720), fgcolor=(1., 1., 1.), bgcolor=(.1, .1, .1))
        mlab.mesh(*X, colormap='ocean')
        mlab.plot3d([e[0] for e in path], [e[1] for e in path], [e[2] for e in path], tube_radius=0.01)

        view = (90., 60, 5, (0, 0, 0))
        mlab.view(*view, reset_roll=True)
        mlab.show()