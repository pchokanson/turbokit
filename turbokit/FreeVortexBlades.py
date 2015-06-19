# FreeVortexBlades.py
# Copyright (c) 2014, 2015 Peter Hokanson
# Vertical Limit Labs

import os, sys, shutil
import math
from subprocess import call, check_call, check_output

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

import PyFoam as pf
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.STLFile import STLFile


from Splines import *
from FreeVortex import FreeVortex
import stl_writer

def condense_face(face):
	"""Remove consecutive duplicate vertices from a face, so we don't have extra
	triangles in the STL file."""
	if len(face) == 3:
		return face
	if len(face) == 4:
		for i in range(len(face)):
			if face[i] == face[i-1]:
				del face[i]
				#print("removed vertex")
				return face
		return face
	else:
		raise Exception("Face has incorrect number of vertices")

class FreeVortexBlades(FreeVortex):
	"""Subclass of FreeVortex meant to implement bladed flow shapes"""
	def __init__(self, 
	             Z=7, 
	             Omega=7330.0,
	             thickness_fn_l=lambda m,s: 0 if m == 0 or m == 1 else 0.001, 
	             thickness_fn_t=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             interblade_faces = 6,
	             hub_solid = True,
	             shroud_solid = False,
	             **kwargs):
		"""Create a representation of free-vortex flow through a bladed region.
		
		Keyword arguments:
		(same as FreeVortex)
		Z -- blade count
		Omega -- angular velocity
		thickness_fn_l -- function for leading edge offset from blade centerline
		thickness_fn_t -- function for trailing edge offset from blade centerline
		interblade_faces -- number of faces between blades
		hub_solid -- whether to make a solid region on the hub
		shroud_solid -- whether to make a solid region for the shroud"""
		super(FreeVortexBlades, self).__init__(**kwargs)
		self.Z = Z
		self.Omega = Omega
		self.thickness_fn_l = thickness_fn_l
		self.thickness_fn_t = thickness_fn_t
		self.hub_solid = hub_solid
		self.shroud_solid = shroud_solid
		self.interblade_faces = interblade_faces
		
		assert hub_solid and not shroud_solid, "Hub/shroud surface options not implemented"
		
		self.makeBladeProfile()
		self.makeMesh()
	
	def makeBladeProfile(self):
		"""Calculate the angular position of the blade at each point (m, s).  This
		is done by numerically integrating the relative velocity."""
		# This is a bit of a hack, but we only have midpoint values and we need
		# to interpolate points outside the convex hull
		interp = scipy.interpolate.NearestNDInterpolator(self.rz_points, self.u_rtz_points)
		
		th = np.zeros(self.r.shape)
		beta = np.zeros(self.r.shape)
		for s in range(self.r.shape[1]):
			for m in range(1, self.r.shape[0]):
				midpoint = np.array([[(self.r[m,s]+self.r[m-1,s])/2, (self.z[m,s]+self.z[m-1,s])/2]])
				u_midpoint = interp(midpoint)[0] # Have to pass in an array
				# Relative velocity terms:
				w_m = math.sqrt(u_midpoint[0]**2 + u_midpoint[2]**2)
				w_th = (u_midpoint[1] - self.Omega * self.r[m,s])
				# Linear displacement from previous grid point
				x_m = math.sqrt((self.r[m,s] - self.r[m-1,s])**2 + (self.z[m,s] - self.z[m-1,s])**2)
				# Final blade angular position and angle
				th[m,s] = th[m-1,s] + x_m * w_th / (self.r[m,s] * w_m)
				beta[m,s] = math.atan2(w_m, w_th)
			th[:,s] -= th[-1,s] # readjust the columns so the outlet side is aligned at 0
		self.th = th
		self.beta = beta
		self.beta[0,:] = self.beta[1,:] # slightly better than using zero, still not perfect
	
	def makeMesh(self):
		"""Enumerate all of the faces required to make a mesh."""
		# NOTE: makes solid-centered meshes for now.  Not desirable for stators, typically.
		# NOTE: Probably swaps thickness functions when Omega is negative
		self.faces = []

		
		th_l = np.copy(self.th)
		th_t = np.copy(self.th)
		for s in range(0, self.r.shape[1]):
			for m in range(0, self.r.shape[0]):
				s_n = s / (self.r.shape[1]-1)
				m_n = m / (self.r.shape[0]-1)
				thickness_l = self.thickness_fn_l(m_n,s_n)
				thickness_t = self.thickness_fn_t(m_n,s_n)
				dth_l = (thickness_l) * math.sin(self.beta[m,s]) / self.r[m,s]
				dth_t = -(thickness_t) * math.sin(self.beta[m,s]) / self.r[m,s]
				th_l[m,s] += dth_l
				th_t[m,s] += dth_t
		self.th_l = th_l
		self.th_t = th_t
		
		# Make the actual mesh
		for i in range(self.Z): 
			# For each blade
			th_i = 2 * math.pi * i / self.Z
			th_next = 2 * math.pi * (i+1) / self.Z
			self.makeBlade(th_i, th_next)
	
	def makeBladeLeadingEdge(self, th_i):
		"""Makes the leading edge faces at m=0."""
		# NOTE: This may be eliminated in favor of zero-thickness at the leading and
		# trailing edges, which would simplify the mesh topology a bit.
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[0,s-1], self.th_l[0,s-1]+th_i, self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s-1], self.th_t[0,s-1]+th_i, self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_t[0,s  ]+th_i, self.z[0,s  ]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_l[0,s  ]+th_i, self.z[0,s  ]])])
	
	def makeBladeTrailingEdge(self, th_i):
		"""Makes the leading edge faces at m=1."""
		# NOTE: This may be eliminated in favor of zero-thickness at the leading and
		# trailing edges, which would simplify the mesh topology a bit.
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[-1,s-1], self.th_l[-1,s-1]+th_i, self.z[-1,s-1]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_l[-1,s  ]+th_i, self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_t[-1,s  ]+th_i, self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s-1], self.th_t[-1,s-1]+th_i, self.z[-1,s-1]])])
	
	def makeBladeLeadingSide(self, th_i):
		"""Make the leading side of the blade (pressure side)."""
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_l[m-1,s-1]+th_i, self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_l[m-1,s  ]+th_i, self.z[m-1,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_l[m  ,s  ]+th_i, self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_l[m  ,s-1]+th_i, self.z[m  ,s-1]])])
	
	def makeBladeTrailingSide(self, th_i):
		"""Make the trailing side of the blade (suction side)."""
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_t[m-1,s-1]+th_i, self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_t[m  ,s-1]+th_i, self.z[m  ,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_t[m  ,s  ]+th_i, self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_t[m-1,s  ]+th_i, self.z[m-1,s  ]])])
	
	def makeBladeShroudEdge(self, th_i):
		"""Make the faces on the shroud side of the blade (s=1).  For unshrouded 
		rotors, this is typically called.  Not typically called for shrouded rotors 
		or stators."""
		for m in range(1, self.r.shape[0]):
			# Faces at blade (shroud) ends
			self.faces.append([rtz_to_xyz([self.r[m-1,-1], self.th_t[m-1,-1]+th_i, self.z[m-1,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_t[m  ,-1]+th_i, self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_l[m  ,-1]+th_i, self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m-1,-1], self.th_l[m-1,-1]+th_i, self.z[m-1,-1]])])
	
	def makeBladeHubEdge(self, th_i):
		"""Make the faces on the hub side of the blade (s=0).  For solid-centered 
		rotors, this isn't called, but typically will be for stators."""
		for m in range(1, self.r.shape[0]):
			# Faces at blade hub ends
			self.faces.append([rtz_to_xyz([self.r[m-1,0], self.th_l[m-1,0]+th_i, self.z[m-1,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_l[m  ,0]+th_i, self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_t[m  ,0]+th_i, self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m-1,0], self.th_t[m-1,0]+th_i, self.z[m-1,0]])])
	
	def makeBladeSpan(self, th_l0, th_t1):
		"""Given the leading edge profile of the current blade, and the trailing
		edge of the next blade, create the junction between them."""
		for m in range(1, self.r.shape[0]):
			
			# Faces between blades
			# theta points on upstream/downstream side of each face
			th_ma = np.linspace(th_l0[m,0], th_t1[m,0], num=self.interblade_faces+1) 
			th_mb = np.linspace(th_l0[m-1,0], th_t1[m-1,0], num=self.interblade_faces+1)
			for j in range(self.interblade_faces):
				# We want to produce a mesh that interpolates between th_l[m0,0]+th_i and th_t[m0,0]+th_next
				# and the same for m0
				self.faces.append([rtz_to_xyz([self.r[m-1,0], th_mb[j]  , self.z[m-1,0]]),
				                   rtz_to_xyz([self.r[m  ,0], th_ma[j]  , self.z[m  ,0]]),
				                   rtz_to_xyz([self.r[m  ,0], th_ma[j+1], self.z[m  ,0]]),
				                   rtz_to_xyz([self.r[m-1,0], th_mb[j+1], self.z[m-1,0]])])
	
	def makeBlade(self, th_i, th_next):
		
		# For convenience
		th_l = self.th_l
		th_t = self.th_t
		
		#self.makeBladeLeadingEdge(th_i)
		#self.makeBladeTrailingEdge(th_i)
		self.makeBladeLeadingSide(th_i)
		self.makeBladeTrailingSide(th_i)
		self.makeBladeShroudEdge(th_i)
		#self.makeBladeHubEdge(th_i)
		
		self.makeBladeSpan(th_l + th_i, th_t + th_next)
		
		# Blade sides
		

		for j in range(self.interblade_faces):
			# Cap at inlet
			th_ma = np.linspace(th_l[0,0]+th_i, 
			                    th_t[0,0]+th_next, 
			                    num=self.interblade_faces+1) 
			self.faces.append([rtz_to_xyz([self.r[0,0], th_ma[j]  , self.z[0,0]]),
			                   rtz_to_xyz([self.r[0,0], th_ma[j+1]  , self.z[0,0]]),
			                   rtz_to_xyz([0, 0, self.z[0,0]])])
			# Cap at outlet
			th_mb = np.linspace(th_l[-1,0]+th_i, 
			                    th_t[-1,0]+th_next, 
			                    num=self.interblade_faces+1)
			self.faces.append([rtz_to_xyz([self.r[-1,0], th_mb[j+1]  , self.z[-1,0]]),
			                   rtz_to_xyz([self.r[-1,0], th_mb[j]  , self.z[-1,0]]),
			                   rtz_to_xyz([0, 0, self.z[-1,0]])])
	
	def writeStlMesh(self, outfilename):
		"""Write out an STL file from the face data."""
		stl_f = open(outfilename, "wb")
		stl = stl_writer.Binary_STL_Writer(stl_f)
		print("Writing STL with %d faces" % len(self.faces))
		for quad in self.faces:
			stl.add_face(condense_face(quad))
		stl.close()

if __name__ == "__main__":
	fvb_rotor = FreeVortexBlades()
	fvb_rotor.writeStlMesh("rotormesh.stl")
	