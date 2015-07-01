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
	             m_max = 1):
		self.thickness_fn_l = thickness_fn_l
		self.thickness_fn_t = thickness_fn_t
		self.m_min = m_min # TODO: not implemented
		self.m_max = m_max # TODO: not implemented

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
	geometry, which extends beyond partial blades.

	This represents only the geometry of the leading and trailing edges (inlet
	and outlet sides) and the front and back faces of the blade.  The separate
	BladeCompleterBase objects are responsible for the geometry between blades:
	for example a solid hub, or the edge of an unshrouded blade."""
	def __init__(self, r, z, th_l, th_t):
		self.r = r
		self.z = z
		self.th_l = th_l
		self.th_t = th_t

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

	def makeBladeFaces(self):
		self.faces = []
		self.makeBladeLeadingEdge()
		self.makeBladeTrailingEdge()
		self.makeBladeLeadingSide()
		self.makeBladeTrailingSide()
		return self.faces

class BladeCompleterBase(object):
	def __init__(self, blades, r, z, s):
		self.s = s
		self.r = r
		self.z = z
		assert self.s == 0 or self.s == 1, "BladeCompleterBase called for neither hub nor shroud"
		self.blades = blades
		self.makeFaces()

class BladeEdgeCompleter(BladeCompleterBase):
	def makeFaces(self):
		self.faces = []
		for blade in self.blades:
			self.makeBladeEdge(blade)

	def makeBladeEdge(self, b):
		for m in range(1, self.r.shape[0]):
			if self.s == 0:
				# Faces at blade hub ends
				self.faces.append([rtz_to_xyz([self.r[m-1,0], b.th_l[m-1,0], self.z[m-1,0]]),
				                   rtz_to_xyz([self.r[m  ,0], b.th_l[m  ,0], self.z[m  ,0]]),
				                   rtz_to_xyz([self.r[m  ,0], b.th_t[m  ,0], self.z[m  ,0]]),
				                   rtz_to_xyz([self.r[m-1,0], b.th_t[m-1,0], self.z[m-1,0]])])
			elif self.s == 1:
				# Faces at blade shroud ends
				self.faces.append([rtz_to_xyz([self.r[m-1,-1], b.th_t[m-1,-1], self.z[m-1,-1]]),
				                   rtz_to_xyz([self.r[m  ,-1], b.th_t[m  ,-1], self.z[m  ,-1]]),
				                   rtz_to_xyz([self.r[m  ,-1], b.th_l[m  ,-1], self.z[m  ,-1]]),
				                   rtz_to_xyz([self.r[m-1,-1], b.th_l[m-1,-1], self.z[m-1,-1]])])

class BladeHubCompleter(BladeCompleterBase):
	def __init__(self, blades, r, z, s, interblade_faces = 6):
		self.interblade_faces = interblade_faces
		super(BladeHubCompleter, self).__init__(blades, r, z, s)

	def makeFaces(self):
		self.faces = []
		self.makeBladeSpans()
		self.makeInletCap()
		self.makeOutletCap()

	def makeInletCap(self):
		for i in range(1, len(self.blades)):
			self.makeInletCapSingleBlade(self.blades[i-1], self.blades[i])
		self.makeInletCapSingleBlade(self.blades[-1], self.blades[0], wrapBlade=True)

	def makeInletCapSingleBlade(self, b0, b1, wrapBlade = False):
		# wrapBlade indicates whether to subtract 2*pi from the blade pair,
		# otherwise there is one that interpolates the long way around.
		th_ma = np.linspace(b0.th_l[0,0],
		                    b1.th_t[0,0] + wrapBlade * 2 * np.pi,
		                    num=self.interblade_faces+1)
		for j in range(self.interblade_faces):
			# Cap at inlet
			self.faces.append([rtz_to_xyz([self.r[0,0], th_ma[j]  , self.z[0,0]]),
			                   rtz_to_xyz([self.r[0,0], th_ma[j+1]  , self.z[0,0]]),
			                   rtz_to_xyz([0, 0, self.z[0,0]])])

	def makeOutletCap(self):
		for i in range(1, len(self.blades)):
			self.makeOutletCapSingleBlade(self.blades[i-1], self.blades[i])
		self.makeOutletCapSingleBlade(self.blades[-1], self.blades[0], wrapBlade=True)

	def makeOutletCapSingleBlade(self, b0, b1, wrapBlade = False):
		pass
		th_mb = np.linspace(b0.th_l[-1,0],
		                    b1.th_t[-1,0] + wrapBlade * 2 * np.pi,
		                    num=self.interblade_faces+1)
		for j in range(self.interblade_faces):
			# Cap at outlet
			self.faces.append([rtz_to_xyz([self.r[-1,0], th_mb[j+1]  , self.z[-1,0]]),
			                   rtz_to_xyz([self.r[-1,0], th_mb[j]  , self.z[-1,0]]),
			                   rtz_to_xyz([0, 0, self.z[-1,0]])])

	def makeBladeSpans(self):
		for i in range(1, len(self.blades)):
			self.makeBladeSpan(self.blades[i-1].th_l, self.blades[i].th_t)
		# Finally a little black magic to make the last inter-blade span work
		# correctly.  The 2*pi offset ensures that they don't wrap the linear
		# interpolation.
		self.makeBladeSpan(self.blades[-1].th_l - 2 * np.pi, self.blades[0].th_t)

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
