
import numpy as np
import math

# The top two functions don't belong here, but don't justify their own files
def intersection_2d(a1, a2, b1, b2):
	"""Find the intersection of the lines defined a1-a2 and b1-b2 or returns false"""
	D = (a1[0]-a2[0])*(b1[1]-b2[1]) - (a1[1]-a2[1])*(b1[0]-b2[0])
	if D == 0:
		return False
	else:
		c0 = ((b1[0]-b2[0])*(a1[0]*a2[1]-a1[1]*a2[0]) - 
		      (a1[0]-a2[0])*(b1[0]*b2[1]-b1[1]*b2[0])) / D
		c1 = ((b1[1]-b2[1])*(a1[0]*a2[1]-a1[1]*a2[0]) - 
		      (a1[1]-a2[1])*(b1[0]*b2[1]-b1[1]*b2[0])) / D
		return np.array([[c0],[c1]])
	
def rtz_to_xyz(rtz):
	return [rtz[0] * math.cos(rtz[1]), rtz[0] * math.sin(rtz[1]), rtz[2]]


class BezierCurve(object):
	"""Represents a 1-D Bezier curve with the given control points."""
	def __init__(self, ctrlpoints):
		self.ctrlpoints = ctrlpoints
		self.order = ctrlpoints.shape[0]
		#print(self.ctrlpoints[0])
	
	def __call__(self, u):
		assert 0 <= u <= 1, "Bezier curve parameter u=%f out of range." % u
		if self.order == 2: # linear
			return (1-u) * self.ctrlpoints[0] + \
			       u * self.ctrlpoints[1]
		elif self.order == 3: # quadratic
			return (1-u)**2 * self.ctrlpoints[0] + \
			       2*u*(1-u) * self.ctrlpoints[1] + \
			       u**2 * self.ctrlpoints[2]
		elif self.order == 4: # cubic
			return (1-u)**3 * self.ctrlpoints[0] + \
			       3*u*((1-u)**2) * self.ctrlpoints[1] + \
			       3*(u**2)*(1-u) * self.ctrlpoints[2] + \
			       u**3 * self.ctrlpoints[3]
		else:
			raise NotImplemented("BezierCurve only supports up to cubic curves currently")

class BezierSurface(object):
	"""Represents a 2-D Bezier curve with the given control points."""
	def __init__(self, ctrlpoints):
		"""Construct from a 2D rectangular array of control points."""
		self.ctrlpoints = ctrlpoints
		self.order_u = ctrlpoints.shape[0]
		self.order_v = ctrlpoints.shape[1]
		#print("BezierSurface(u=%d, v=%d)" % (self.order_u, self.order_v))
	
	def __call__(self, u, v):
		# This function doesn't claim to be a terribly efficient solution
		# Graphics applications typically use De Casteljau's algorithm.
		
		# Reduce the problem to the 1D case: solve along v for a 1D Bezier curve in u
		
		k_u = []
		for j in range(self.order_u):
			b_v = BezierCurve(self.ctrlpoints[j])
			k_u.append(b_v(v))
		b_u = BezierCurve(np.array(k_u))
		return b_u(u)

if __name__ == "__main__":
	a1 = np.array([-1,1])
	a2 = np.array([-1,-1])
	b1 = np.array([0,2])
	b2 = np.array([2,2])
	print(intersection_2d(a1, a2, b1, b2))