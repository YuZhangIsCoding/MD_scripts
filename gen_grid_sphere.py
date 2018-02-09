#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: gen_grid_sphere.py
# Description:  This is a python script that generates the grid for carbon ions
#               for the calculation of electrostatic potential.
#               For the mesh on the sphere, there are several models, including
#               UV sphere and icoshpere. 
#               UV mapping could be created by walking theta and phi.
#               Icosphere can be generate by winding up the isosahedron.
# Date:     01-25-2018  Created
#           01-27-2018  Added icosphere and UV sphere

import numpy as np
import argparse

def get_args():
    '''Get arguments from terminal, otherwise the default parameters will be used'''
    parser = argparse.ArgumentParser(description = 'Specify center and radius, etc')
    parser.add_argument('-c', '--center', dest = 'center', nargs = '*',
                        type = float, help = 'Centers of carbon onion')
    parser.add_argument('-r', '--radius', type = float, default = 1.018,
                        help = 'Radius of the carbon onion')
    parser.add_argument('-o', '--output', dest = 'output', default = 'gridElectrode1.dat',
                        help = 'Output filename')
    parser.add_argument('-s', '--subdivide', dest = 'sub', default = 0, 
                        type = int, help = 'Subdivision')
    parser.add_argument('--choice', dest = 'choice', choices = ('ico', 'uv'),
                       default = 'ico', help = 'Choice of mapping, ico or uv.')
    try:
        __IPYTHON__
        args = parser.parse_args([])
    except NameError:
        args = parser.parse_args()

    if args.center:
        centers = np.array(args.center).reshape(-1, 3)
    else:
        centers = [[4.386006683, 4.392999283, 4.37500211],
                [4.386006683, 4.392999283, 13.35800205]]
        centers = np.array(centers).reshape(-1, 3)
    print('The centers are:\n', centers)
    print('The radus is %.3f nm' %args.radius)
    print('Writing to file %s ...' %args.output)
    return centers, args.radius, args.output, args.choice, args.sub

class Vertex(object):
    '''Vertex class is each grid on the surface. Each vertex will have a unique index.
    '''
    index = 0
    midpoints = {}
    def __init__(self, x, y, z):
        '''Normalize the coordinates and assign index to a vertex'''
        length = np.sqrt(x**2+y**2+z**2)
        self.xyz = np.array([x, y, z])/length
        self.index = Vertex.index
        Vertex.index  += 1
    def get_index(self):
        return self.index
    def get_xyz(self):
        return self.xyz
    def dist_from(self, other):
        return np.sqrt(sum((self.get_xyz()-other.get_xyz())**2))
    def __lt__(self, other):
        return self.index < other.get_index()
    def __eq__(self, other):
        return self.index == other.get_index()
    def get_midpoint(self, other):
        key = (min(self, other).get_index(), max(self, other).get_index())
        if key not in Vertex.midpoints.keys():
            Vertex.midpoints[key] = Vertex(*(self.get_xyz()+other.get_xyz())/2)
        return Vertex.midpoints[key]
    def get_dist(v1, v2):
        return np.sqrt(sum((v1.get_xyz()-v2.get_xyz())**2))
    def get_coords(verts):
        return np.array([item.get_xyz() for item in verts])
    def clear_midpoints():
        '''Clear the dictionary'''
        Vertex.midpoints = {}
    def __str__(self):
        temp = (self.index, )+tuple(self.xyz)
        return '%d (%f, %f, %f)' %temp
    def __repr__(self):
        return str(self.index)

class Face(object):
    '''Face is compromised by 3 vertices.'''
    def __init__(self, v1, v2, v3):
        self.verts = [v1, v2, v3]
    def subdivide(self):
        mid = []
        for i in range(3):
            for j in range(i+1, 3):
                mid.append(self.verts[i].get_midpoint(self.verts[j]))
        subs = [Face(self.verts[0], mid[0], mid[1]),
                Face(self.verts[1], mid[0], mid[2]),
                Face(self.verts[2], mid[1], mid[2]),
                Face(*mid)]
        return subs
    def get_verts(self):
        return self.verts
    def __repr__(self):
        return '(%d, %d, %d)' %tuple(item.get_index() for item in self.verts)
    def __str__(self):
        return '''
  %s
 /\t\t\t\t   \\
%s -- %s
        ''' %tuple(self.verts)

def gen_icosahedron():
    '''Generate the vertices of icosahedron, the center is set at (0, 0 ,0)
    The initial coordinates for the 12 vertices are:
    (    0, +/-1, +/-p)
    ( +/-1, +/-p,    0)
    ( +/-p,    0, +/-1)
    where p = (1+sqrt(5))/2
    
    Return a list of Vertex objects.
    '''
    p = (1+np.sqrt(5))/2
    verts = []
    for i in [-1, 1]:
        for j in [-p, p]:
            verts.append(Vertex(0, i, j))
            verts.append(Vertex(i, j, 0))
            verts.append(Vertex(j, 0, i))
    return verts

def judge_face(v1, v2, v3, radius = 1):
    '''Judge if the three vertices are valid to form an initial icosahedron triangle
    face on a sphere with the given radius.
    '''
    return Vertex.get_dist(v1, v2) == Vertex.get_dist(v1, v3) == \
            Vertex.get_dist(v2, v3) == 2*radius/np.sqrt((5+np.sqrt(5))/2)
                

def gen_icosphere(subdivide = 0, radius = 1):
    '''First build an icosahedron, which will produce 12 vertices and 20 faces.
    Subdivide each face into 4 subfaces, and continue this subdivision until the
    designated depth is reached. The final grids consist of inital 12 vertices and
    all the midpoints generated in the process.
    
    Return a numpy array of the icosphere grids.
    '''
    verts = gen_icosahedron()
    Vertex.clear_midpoints()
    faces = []
    ## The following block of code might be optimized
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            for k in range(j+1, len(verts)):
                if judge_face(verts[i], verts[j], verts[k], radius = 1):
                    faces.append(Face(verts[i], verts[j], verts[k]))

    for _ in range(subdivide):
        subfaces = []
        for face in faces:
            subfaces.extend(face.subdivide())
        faces = subfaces
    verts.extend(list(Vertex.midpoints.values()))
    return Vertex.get_coords(verts)

def gen_uvsphere(radius = 1, slices = 30):
    '''The UV sphere will be created by looping over sliced theta and phi
    '''
    coords = []
    for theta in np.linspace(0, np.pi, slices+1):
        for phi in np.linspace(0, 2*np.pi, slices*2+1)[:-1]:
            coords.append([radius*np.sin(theta)*np.cos(phi),
                          radius*np.sin(theta)*np.sin(phi),
                          radius*np.cos(phi)])
    return np.array(coords)*radius

########## main ##########
centers, radius, outname, choice, sub = get_args()
if choice == 'ico':
    grids = gen_icosphere(sub, radius)
else:
    grids = gen_uvsphere(radius = radius)
coords = [grids+center for center in centers]
coords = np.concatenate(tuple(coords), axis = 0)
np.savetxt(outname, coords, fmt = ['%12.5f' for _ in range(3)])
