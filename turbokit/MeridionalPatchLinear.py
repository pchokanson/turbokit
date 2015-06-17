# MeridionalPatchLinear.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

from MeridionalPatch import MeridionalPatch
from Splines import *

class MeridionalPatchLinear(MeridionalPatch):
	
	def __init__(self, m0_s0, m0_s1, m1_s0, m1_s1):
		"""Construct a meridional patch between 2D (r,z) points.
		
		Arguments:
		m0_s0 -- vector at m=0, s=0
		m0_s1 -- vector at m=0, s=1
		m1_s0 -- vector at m=1, s=0
		m1_s1 -- vector at m=1, s=1"""
		self.m0_s0 = m0_s0
		self.m0_s1 = m1_s1
		self.m1_s0 = m1_s0
		self.m1_s1 = m1_s1
		self.k_array = np.array([[m0_s0, m0_s1],
		                         [m1_s0, m1_s1]])
		self.b = BezierPatch(self.k_array)
		
	def __call__(self, m, s):
		return self.b(m, s)
