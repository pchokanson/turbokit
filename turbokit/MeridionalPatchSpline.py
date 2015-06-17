# MeridionalPatchSpline.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

from MeridionalPatch import MeridionalPatch
from Splines import *

class MeridionalPatchSpline(MeridionalPatch):
	
	def __init__(self, m0_s0, m0_s1, m1_s0, m1_s1, v_m0, v_m1):
		"""Construct a meridional patch interpolating between 2D (r,z) points with
		given inlet and outlet velocities (also 2D).
		
		Arguments:
		m0_s0 -- vector at m=0, s=0
		m0_s1 -- vector at m=0, s=1
		m1_s0 -- vector at m=1, s=0
		m1_s1 -- vector at m=1, s=1
		v_m0 -- 2D velocity vector at m=0 (inlet)
		v_m1 -- 2D velocity vector at m=1 (outlet)"""
		self.m0_s0 = m0_s0
		self.m0_s1 = m1_s1
		self.m1_s0 = m1_s0
		self.m1_s1 = m1_s1
		
		s0_ctrlpoint = intersection_2d(m0_s0, m0_s0+v_m0, m1_s0, m1_s0+v_m1)
		s1_ctrlpoint = intersection_2d(m0_s1, m0_s1+v_m0, m1_s1, m1_s1+v_m1)
		
		self.k_array = np.array([[m0_s0, m0_s1],
		                         [s0_ctrlpoint, s1_ctrlpoint],
		                         [m1_s0, m1_s1]])
		self.b = BezierPatch(self.k_array)
		
	def __call__(self, m, s):
		return self.b(m, s)
