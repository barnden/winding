import numpy as np
from mayavi import mlab

vertices = np.array([
    [0., 0., 0.],
    [0., 1., 0.],
    [1., 0., 0.],
    [1., 1., 0.],
    [0., 0., 1.],
    [0., 1., 1.],
    [1., 0., 1.],
    [1., 1., 1.],
])

faces = [
    [0, 1, 2],
    [3, 2, 1],
    [4, 6, 5],
    [7, 5, 6],
    [0, 2, 4],
    [6, 4, 2],
    [0, 4, 1],
    [5, 1, 4],
    [1, 5, 3],
    [7, 3, 5],
    [2, 3, 6],
    [7, 6, 3],
]

mlab.triangular_mesh(*vertices.T, faces, color=(1, 1, 1))
mlab.quiver3d(*centres.T, *normals.T)
mlab.show()