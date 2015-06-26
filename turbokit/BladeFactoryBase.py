# BladeFactoryBase.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

import math

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

from Splines import *

class BladeFactoryBase(object):
	"""Base class for blades in the turbine flow.  Presents as a factory object
	capable of returning all of the blade faces at a given angular offset.  This
	is then called repeatedly to make all of the matching blades."""

	def __init__(self,
	             thickness_fn_l=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             thickness_fn_t=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             m_min = 0,
	             m_max = 1,
	             hubEdge=False,
	             shroudEdge=True):
		self.thickness_fn_l = thickness_fn_l
		self.thickness_fn_t = thickness_fn_t
		self.m_min = m_min # TODO: not implemented
		self.m_max = m_max # TODO: not implemented
		self.hubEdge = hubEdge
		self.shroudEdge = shroudEdge

	def __call__(self, r, z, th, beta):
		th_l = np.copy(th)
		th_t = np.copy(th)
		for s in range(0, r.shape[1]):
			for m in range(0, r.shape[0]):
				s_n = s / (r.shape[1]-1)
				m_n = m / (r.shape[0]-1)
				thickness_l = self.thickness_fn_l(m_n,s_n)
				thickness_t = self.thickness_fn_t(m_n,s_n)
				dth_l = (thickness_l) * math.sin(beta[m,s]) / r[m,s]
				dth_t = -(thickness_t) * math.sin(beta[m,s]) / r[m,s]
				th_l[m,s] += dth_l
				th_t[m,s] += dth_t
		return BladeBase(r, z, th_l, th_t)

class BladeBase(object):
	"""Class for representing blade geometry.  Calculates faces and stores
	geometry.  Generally meant to be produced by below factory type.  Blade
	may or may not encompass entire flow region (as in a splitter blade in a
	centrifugal pump or compressor).

	Note that r and z may be limited to a subset of the flow, but th_l and th_t
	should be the full size, since they are referenced to create the inter-blade
	geometry, which extends beyond partial blades."""
	def __init__(self, r, z, th_l, th_t,
	             hubEdge = False,
	             shroudEdge = True):
		self.r = r
		self.z = z
		self.th_l = th_l
		self.th_t = th_t
		self.hubEdge = hubEdge
		self.shroudEdge = shroudEdge

		self.makeBladeFaces()

	def makeBladeLeadingEdge(self):
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[0,s-1], self.th_l[0,s-1], self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s-1], self.th_t[0,s-1], self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_t[0,s  ], self.z[0,s  ]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_l[0,s  ], self.z[0,s  ]])])

	def makeBladeTrailingEdge(self):
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[-1,s-1], self.th_l[-1,s-1], self.z[-1,s-1]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_l[-1,s  ], self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_t[-1,s  ], self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s-1], self.th_t[-1,s-1], self.z[-1,s-1]])])

	def makeBladeLeadingSide(self):
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_l[m-1,s-1], self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_l[m-1,s  ], self.z[m-1,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_l[m  ,s  ], self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_l[m  ,s-1], self.z[m  ,s-1]])])

	def makeBladeTrailingSide(self):
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_t[m-1,s-1], self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_t[m  ,s-1], self.z[m  ,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_t[m  ,s  ], self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_t[m-1,s  ], self.z[m-1,s  ]])])

	def makeBladeShroudEdge(self):
		for m in range(1, self.r.shape[0]):
			# Faces at blade (shroud) ends
			self.faces.append([rtz_to_xyz([self.r[m-1,-1], self.th_t[m-1,-1], self.z[m-1,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_t[m  ,-1], self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_l[m  ,-1], self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m-1,-1], self.th_l[m-1,-1], self.z[m-1,-1]])])

	def makeBladeHubEdge(self):
		for m in range(1, self.r.shape[0]):
			# Faces at blade hub ends
			self.faces.append([rtz_to_xyz([self.r[m-1,0], self.th_l[m-1,0], self.z[m-1,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_l[m  ,0], self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_t[m  ,0], self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m-1,0], self.th_t[m-1,0], self.z[m-1,0]])])

	def makeBladeFaces(self):
		self.faces = []
		self.makeBladeLeadingEdge()
		self.makeBladeTrailingEdge()
		self.makeBladeLeadingSide()
		self.makeBladeTrailingSide()
		if self.shroudEdge:
			self.makeBladeShroudEdge()
		if self.hubEdge:
			self.makeBladeHubEdge()
		# for face in self.faces:
		# 	print(face, len(face))
		return self.faces
