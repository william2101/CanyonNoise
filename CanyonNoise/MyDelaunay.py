# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:47:25 2023

@author: willi
"""

from scipy.spatial import Delaunay
import numpy as np

class MyDelaunay(Delaunay):
    """Add helper functions"""
    
    def __init__(self, *args, pts = np.asarray([0,0,0])):
        """pts is the set of xyz pts"""
        Delaunay.__init__(self, *args)
        self.pts = pts

    def centroid(self):
        """Centroids of delauney triangularization
        """
        s = self.simplices
        p = self.pts
        
        num_tri = len(s)
        
        cen = np.zeros([num_tri,3])
        
        for i in range(0,num_tri):
            cen[i,0] = sum(p[s[i],0])/3
            cen[i,1] = sum(p[s[i],1])/3
            cen[i,2] = sum(p[s[i],2])/3
       
        return cen

    def normal(self):
        """Upward normals of delauney triangularization
        """
        s = self.simplices
        p = self.pts
        
        num_tri = len(s)
        
        normal = np.zeros([num_tri,3])
        
        for i in range(0,num_tri):
            ##Choose any two edges to find normal
            pt0 = p[s[i,0]]
            pt1 = p[s[i,1]]
            pt2 = p[s[i,2]]

            ##Edge vector given by difference between two points
            cp = np.cross(pt0 - pt1, pt1 - pt2) ##cross product
            normal[i] = cp / np.linalg.norm(cp)
       
        return normal
    
    def distance(self, source):
        """Return the distance between a source point
        and the centroid of every point on the mesh"""
        centroids = self.centroid()
        
        dist = np.zeros(len(centroids))
        
        for i in range(0,len(centroids)):
            dist[i] = np.linalg.norm(source - centroids[i])
        
        return dist
    
    def angle(self, source):
        """Return the angle between the vector defined 
        by source - centroid versus the normal vector
        of the triangular mesh element"""
        
        centroids = self.centroid()
        normals = self.normal()
        
        ang = np.zeros(len(centroids))
        
        for i in range(0,len(centroids)):
            s_mag = np.linalg.norm(source - centroids[i])
            #Magnitude of the normal vector is defined to be 1
            
            ang[i] = np.arccos(np.dot(source-centroids[i], normals[i]) / s_mag )
        
        return ang
    
    def area(self):
        """Get the area of the simplices"""
        s = self.simplices
        p = self.pts
        
        num_tri = len(s)
        
        area = np.zeros(num_tri)
        
        for i in range(0,num_tri):
            ##Choose any two edges to find normal
            pt0 = p[s[i,0]]
            pt1 = p[s[i,1]]
            pt2 = p[s[i,2]]

            ##Edge vector given by difference between two points
            cp = np.cross(pt0 - pt1, pt1 - pt2) ##cross product
            area[i] = 0.5 * np.linalg.norm(cp)
       
        return area
        
    
    def plane_height(self):
        """Return a vector, d, whose value is the height 
        of a plane given by n dot point = d. Where n is 
        the normal vector, point is the centroid, and dot 
        is the dot product. I.e., nx * x + ny * y + nz * z = d
        """ 

        return np.dot( self.centroid(), self.normal() )
    
    def check_blocked(self, source, step_factor = 10):
        """Check if the vector between the source and 
        a simplex on the mesh is blocked by any other
        simplices.
        
        To do so, parametrize the vector between source
        and destination simplex. Step through the line 
        and identify from the x,y coordinates on the line
        which simplex it resides in (in 2d space). Then check
        if the z value of the line is less than the z value of the
        identified simplex. If so, the line is under the graph,
        and therefore has crossed a boundary
        """
        centroids = self.centroid()
        
        check = np.ones(len(centroids)) 
        
        ##Set the step size = 5* sqrt(of average triangle area) / distance
        avg_len = np.average(self.area()) ** (1/2) 
        
        #Set source as starting point and the search simplex as end
        
        r0 = source
        for i in range(0, len(centroids)):
            r1 = centroids[i]
            
            step = step_factor * avg_len / np.linalg.norm(r1 - r0)
            
            #Initialize the check point as r0
            check_pt = r0 
            
            while (np.linalg.norm(check_pt - r0) < np.linalg.norm(r1 - r0)):
                
                simplex = self.find_simplex(check_pt[0:2])
                
                ##Check the check point against lowest point in simplex                
                if (check_pt[2] < self.pts[self.simplices[simplex],2].min()):
                    check[i] = 0
                    break
                
                check_pt = check_pt + step * (r1 - r0)

        return check
    
    def element_blocking_matrix(self):
        """Return a matrix which identifies the element 
        to element visibility. I.e., check if a line 
        between two elements is blocked by any other
        elements
        """
        centroids = self.centroid()
        block_mat = np.zeros( [len(centroids), len(centroids)] )
        step_factor = 10
        
        ##Set the step size 
        avg_len = np.average(self.area()) ** (1/2) 
        
        for i in range (0, len(centroids)):
            r0 = centroids[i]
            print(i)
            
            for j in [x for x in reversed(range(i, len(centroids))) if x != i]:
                block_mat[i,j] = 1
                
                r1 = centroids[j]
                step = step_factor * avg_len / np.linalg.norm(r1 - r0)
                
                check_pt = r0
                
                while (np.linalg.norm(check_pt - r0) < np.linalg.norm(r1 - r0)):
                    
                    simplex = self.find_simplex(check_pt[0:2])
                    
                    ##Check the check point against lowest point in simplex                
                    if (check_pt[2] < self.pts[self.simplices[simplex],2].min()):
                        block_mat[i,j] = 0
                        break
                    
                    check_pt = check_pt + step * (r1 - r0)
        
        return block_mat + block_mat.T
        

    def reflection_vector_from_source(self, source):
        """Return the reflected angle from the source
        """
        centroids = self.centroid()
        normals = self.normal()
        
        #r is reflection vector
        r = np.zeros([len(centroids),3]) 
        
        for i in range(0, len(centroids)):
            ##Let d be the vector between source and centroid
            d = centroids[i] - source
            
            r[i] = d - 2 * np.dot(normals[i], d) * normals[i]
            r[i] = r[i] / np.linalg.norm(r[i])
            
        return r
    
    def reflection_intensity(self, source, blocking_matrix, I_incident = 1, dispersion_angle = 0):
        """Return the intensity at every simplex from the
        first reflection from the source
        source = source location
        I_incident = incident intensity
        dispersion_angle = width of the specular reflection in radians
        """
        centroids = self.centroid()
        reflecs = self.reflection_vector_from_source(source)
        cos_dispersion = np.cos(dispersion_angle)
        diffuse_loss = 100
                
        #intensity map
        I = np.zeros(len(centroids))
        
        #Loop through every reflection (secondary source) element
        for i in range(0, len(centroids)):
            #Loop through every destination element
            for j in [x for x in range(0, len(centroids)) if x != i]:
                if blocking_matrix[i,j]: ##Check if elements can see each other
                    v_mag = np.linalg.norm( centroids[i] - centroids[j] )
                    v = (centroids[j] - centroids[i]) / v_mag
                
                    #Specular reflection
                    if (abs(np.dot(v, reflecs[i])) > cos_dispersion):
                        I[j] = I[j] + I_incident[i] / v_mag ** 2
                        
                    #Diffuse reflection
                    else:
                        I[j] = I[j] + I_incident[i] / v_mag ** 2 / diffuse_loss
                    
        return I
                
                
            
        
    
    
    
        
        
    
    