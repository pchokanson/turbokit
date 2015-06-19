# FreeVortexBlades.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

class BladeBase(object):
	"""Base class for blades in the turbine flow.  Presents as a factory object
	capable of returning all of the blade faces at a given angular offset.  This
	is then called repeatedly to make all of the matching blades."""
	
	def __init__(self, r, z, th, beta,
	             thickness_fn_l=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             thickness_fn_t=lambda m,s: 0 if m == 0 or m == 1 else 0.001,
	             hubEdge=False, 
	             shroudEdge=False):
		self.r = r
		self.z = z
		self.th = th
		self.beta = beta
		self.hubEdge = hubEdge
		self.shroudEdge = shroudEdge
		
		self.calcVertices()
		self.faces = []

	def calcVertices(self):
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

	def makeBladeLeadingEdge(self, th_i):
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[0,s-1], self.th_l[0,s-1]+th_i, self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s-1], self.th_t[0,s-1]+th_i, self.z[0,s-1]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_t[0,s  ]+th_i, self.z[0,s  ]]),
			                   rtz_to_xyz([self.r[0,s  ], self.th_l[0,s  ]+th_i, self.z[0,s  ]])])
	
	def makeBladeTrailingEdge(self, th_i):
		for s in range(1, self.r.shape[1]):
			self.faces.append([rtz_to_xyz([self.r[-1,s-1], self.th_l[-1,s-1]+th_i, self.z[-1,s-1]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_l[-1,s  ]+th_i, self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s  ], self.th_t[-1,s  ]+th_i, self.z[-1,s  ]]),
			                   rtz_to_xyz([self.r[-1,s-1], self.th_t[-1,s-1]+th_i, self.z[-1,s-1]])])
	
	def makeBladeLeadingSide(self, th_i):
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_l[m-1,s-1]+th_i, self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_l[m-1,s  ]+th_i, self.z[m-1,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_l[m  ,s  ]+th_i, self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_l[m  ,s-1]+th_i, self.z[m  ,s-1]])])
	
	def makeBladeTrailingSide(self, th_i):
		for m in range(1, self.r.shape[0]):
			for s in range(1, self.r.shape[1]):
				self.faces.append([rtz_to_xyz([self.r[m-1,s-1], self.th_t[m-1,s-1]+th_i, self.z[m-1,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s-1], self.th_t[m  ,s-1]+th_i, self.z[m  ,s-1]]),
				                   rtz_to_xyz([self.r[m  ,s  ], self.th_t[m  ,s  ]+th_i, self.z[m  ,s  ]]),
				                   rtz_to_xyz([self.r[m-1,s  ], self.th_t[m-1,s  ]+th_i, self.z[m-1,s  ]])])
	
	def makeBladeShroudEdge(self, th_i):
		for m in range(1, self.r.shape[0]):
			# Faces at blade (shroud) ends
			self.faces.append([rtz_to_xyz([self.r[m-1,-1], self.th_t[m-1,-1]+th_i, self.z[m-1,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_t[m  ,-1]+th_i, self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m  ,-1], self.th_l[m  ,-1]+th_i, self.z[m  ,-1]]),
			                   rtz_to_xyz([self.r[m-1,-1], self.th_l[m-1,-1]+th_i, self.z[m-1,-1]])])
	
	def makeBladeHubEdge(self, th_i):
		for m in range(1, self.r.shape[0]):
			# Faces at blade hub ends
			self.faces.append([rtz_to_xyz([self.r[m-1,0], self.th_l[m-1,0]+th_i, self.z[m-1,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_l[m  ,0]+th_i, self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m  ,0], self.th_t[m  ,0]+th_i, self.z[m  ,0]]),
			                   rtz_to_xyz([self.r[m-1,0], self.th_t[m-1,0]+th_i, self.z[m-1,0]])])

	def makeBlade(self, th_i, th_next):
		# For convenience
		th_l = self.th_l
		th_t = self.th_t
		
		self.makeBladeLeadingEdge(th_i)
		self.makeBladeTrailingEdge(th_i)
		self.makeBladeLeadingSide(th_i)
		self.makeBladeTrailingSide(th_i)
		self.makeBladeShroudEdge(th_i)




