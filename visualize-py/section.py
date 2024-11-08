import numpy as np
import re
from mayavi import mlab
from scipy.spatial import ConvexHull

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

fig = mlab.figure(size=(800,800))

surface = BSpline("../build/ideformed-spline-step-0.txt")
# mlab.mesh(*surface.generate_mesh(128), color=(0.6, 0.65, 0.8))

"""
// Project into 2D plane containing the face
    Vec3 n = normal(state);
    Vec3 const& p0 = state.data0[m_v0];
    Vec3 p_prime = p - ((p - p0).dot(n) * n) - p0;

    Vec3 u = (state.data0[m_v1] - p0).normalized();
    Vec3 v = u.cross(n);

    return { p_prime.dot(u), p_prime.dot(v) };
"""

def compute_plane(section, N):
    p0 = section["curve"][:, 0]
    p1 = section["curve"][:, N // 3]
    p2 = section["curve"][:, (2 * N) // 3]
    e1 = normalize(p1 - p0)

    n = normalize(np.cross(e1, p2 - p0))

    return (n, p0, e1, np.cross(e1, n))

def ortho_projection(n, p):
    return normalize(n - (np.dot(n, p) * p))

def project_onto_plane(n, p0, e1, e2, p):
    p_prime = p - (np.dot(p - p0, n) * n) - p0

    return (np.dot(p_prime, e1), np.dot(p_prime, e2))

N = 128
for v in [8]: #range(int(surface.v_min), int(surface.v_max), 1):
    cpts = surface.points[v - 1]
    section = surface.generate_slice(v, N)
    plane = compute_plane(section, N)

    mlab.plot3d(*section["curve"], color=(1., 0., 0.), tube_radius=0.005)
    points = section["curve"].T
    projected_points = [project_onto_plane(*plane, p) for p in points]

    cv = ConvexHull(projected_points)
    convex_hull = np.array([points[i] for i in cv.vertices] + [points[cv.vertices[0]]]).T

    mlab.plot3d(*convex_hull, color=(1., 1., 0.), tube_radius=0.0055)

    for i in [13]:
        cp0 = cpts[i]
        cp1 = cpts[(i + 1) % len(cpts)]
        cp2 = cpts[(i + 2) % len(cpts)]
        cp3 = cpts[(i + 3) % len(cpts)]

        polygon = [cp0, cp1, cp2, cp3, cp0]

        mlab.plot3d(*np.array(polygon).T, color=(1., 1., 1.), tube_radius=0.005)

    i = 0
    for (point, normal) in zip(section["curve"].T, section["normals"].T):
        if (i := i + 1) % 4 == 0:
            x, y, z = point

            px, py, pz = ortho_projection(normal, plane[0]) / 5
            nx, ny, nz = normal / 5

            mlab.plot3d([x, x + nx], [y, y + ny], [z, z + nz], color=(0., 0., 1.), tube_radius=0.005)
            mlab.plot3d([x, x + px], [y, y + py], [z, z + pz], color=(0., 1., 0.), tube_radius=0.005)

mlab.show()