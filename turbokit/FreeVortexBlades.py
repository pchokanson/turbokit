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
from BladeFactoryBase import BladeFactoryBase, BladeBase, BladeCompleterBase,\
                             BladeEdgeCompleter, BladeHubCompleter
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
		raise Exception("Face %s has incorrect number of vertices" % face)

class FreeVortexBlades(FreeVortex):
	"""Subclass of FreeVortex meant to implement bladed flow shapes"""
	def __init__(self,
	             Z=7,
	             bladeFactories = None,
	             Omega=7330.0,
	             thickness_fn_l=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             thickness_fn_t=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             interblade_faces = 6,
	             **kwargs):
		"""Create a representation of free-vortex flow through a bladed region.

		Keyword arguments:
		(same as FreeVortex)
		Z -- blade count
		bladeFactories -- factory objects for making each blade
		Omega -- angular velocity
		thickness_fn_l -- function for leading edge offset from blade centerline
		thickness_fn_t -- function for trailing edge offset from blade centerline
		interblade_faces -- number of faces between blades
		hub_solid -- whether to make a solid region on the hub
		shroud_solid -- whether to make a solid region for the shroud"""
		super(FreeVortexBlades, self).__init__(**kwargs)
		self.Z = Z

		if bladeFactories is not None:
			self.bladeFactories = bladeFactories
		else:
			self.bladeFactories = [BladeFactoryBase() for i in range(self.Z)]

		self.Omega = Omega
		self.thickness_fn_l = thickness_fn_l
		self.thickness_fn_t = thickness_fn_t
		self.interblade_faces = interblade_faces

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
		# NOTE: Probably swaps thickness functions when Omega is negative
		self.faces = []
		self.blades = []

		for i in range(0, self.Z):
			th_i = i * 2 * np.pi / self.Z
			blade = self.bladeFactories[i](self.r, self.z, self.th + th_i, self.beta)
			self.faces.extend(blade.makeBladeFaces())
			self.blades.append(blade)

		self.hubCompleter = BladeHubCompleter(self.blades, self.r, self.z, 0)
		self.faces.extend(self.hubCompleter.faces)
		self.shroudCompleter = BladeEdgeCompleter(self.blades, self.r, self.z, 1)
		self.faces.extend(self.shroudCompleter.faces)

	def writeStlMesh(self, outfilename):
		"""Write out an STL file from the face data."""
		stl_f = open(outfilename, "wb")
		stl = stl_writer.Binary_STL_Writer(stl_f)
		print("Writing STL with %d faces" % len(self.faces))
		for quad in self.faces:
			stl.add_face(condense_face(quad))
		stl.close()

if __name__ == "__main__":
	fvb_rotor = FreeVortexBlades(points_m=15, points_s = 10)
	fvb_rotor.writeStlMesh("rotormesh.stl")
