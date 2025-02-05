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

def group_consecutive(v):
    groups = []

    # see: https://stackoverflow.com/questions/73978318/splitting-a-list-on-non-sequential-numbers
    for _, g in groupby(enumerate(v), lambda x: x[0] - x[1]):
        groups.append([v for _, v in g])

    return groups

def plot():
    surface_values = np.zeros((surface_mesh.shape[1:]))

    if False:
        # for v in [1]:
        #     section = surface.generate_slice_u(v, 256)
        #     mlab.plot3d(*section["curve"], tube_radius=0.002)

        #     plane = compute_plane(section, 256)

        #     points = section["curve"].T
        #     projected_points = points
        #     # projected_points = [project_onto_plane(*plane, point) for point in points]
        #     hull = ConvexHull(projected_points).vertices

        #     hull_path = np.array([projected_points[i] for i in hull]).T
        #     mlab.plot3d(*hull_path, tube_radius=0.002, color=(1., 1., 0.))

        for iu in range(surface_mesh.shape[2]):
            section = {"curve": surface_mesh[:, :, iu]}

            plane = compute_plane(section, section["curve"].shape[-1])

            points = section["curve"].T
            projected_points = [project_onto_plane(*plane, p) for p in points]

            hull = ConvexHull(projected_points).vertices

            # hull_path = np.array([points[i] for i in hull]).T
            # mlab.plot3d(*hull_path, tube_radius=0.002, color=(1., 1., 0.))

            off_hull = set([*range(len(projected_points))]) - set(hull)
            concave_paths = group_consecutive(off_hull)

            values = np.zeros((len(points)))
            for path in concave_paths:
                if len(path) < 2:
                    continue

                start = points[path[0]]
                end = points[path[-1]]
                basis = normalize(end - start)
                S = [np.dot(point, basis) for point in map(lambda i: points[i], path)]

                values[path[0]:path[-1]] = .5

                increasing = True
                decreasing = True
                last = None

                for i in range(len(S) - 1):
                    c = S[i]
                    n = S[i + 1]

                    if increasing and c > n:
                        increasing = False

                        if not decreasing:
                            last = ["increasing", c, i]

                            continue

                    if decreasing and c < n:
                        decreasing = False

                        if not increasing:
                            last = ["decreasing", c, i]

                            continue

                    if (last is not None) and (last[0] == "increasing" and last[1] <= c):
                        values[path[last[2]]:path[i + 1]] = 1

                        last = None

            surface_values[:, iu] = values

    for iv in range(surface_mesh.shape[1]):
        section = {"curve": surface_mesh[:, iv, :]}

        plane = compute_plane(section, section["curve"].shape[-1])

        points = section["curve"].T
        projected_points = [project_onto_plane(*plane, p) for p in points]

        hull = ConvexHull(projected_points).vertices

        off_hull = set([*range(len(projected_points))]) - set(hull)
        concave_paths = group_consecutive(off_hull)

        values = np.zeros((len(points)))
        for path in concave_paths:
            if len(path) < 2:
                continue

            start = points[path[0]]
            end = points[path[-1]]
            basis = normalize(end - start)
            S = [np.dot(point, basis) for point in map(lambda i: points[i], path)]

            values[path[0]:path[-1]] = .5

            increasing = True
            decreasing = True
            last = None

            for i in range(len(S) - 1):
                c = S[i]
                n = S[i + 1]

                if increasing and c > n:
                    increasing = False

                    if not decreasing:
                        last = ["increasing", c, i]

                        continue

                if decreasing and c < n:
                    decreasing = False

                    if not increasing:
                        last = ["decreasing", c, i]

                        continue

                if (last is not None) and (last[0] == "increasing" and last[1] <= c):
                    values[path[last[2]]:path[i + 1]] = 1

                    last = None

        surface_values[iv] = np.maximum(surface_values[iv], values)

    mlab.mesh(*surface_mesh, scalars=surface_values, vmin=0, vmax=1)

    mlab.show()
plot()