# MeridionalPatch.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

class MeridionalPatch(object):
	"""Abstract base class for defining meridional flowfield.  When called as a 
	function, this takes two arguments: m and s, for the inlet-to-outlet and 
	hub-to-shroud parameters, both varying from 0 to 1.  This returns a numpy 
	vector (r,z) for the relevant """
	
	def __init__(self):
		assert False
		
	def __call__(self, m, s):
		assert False