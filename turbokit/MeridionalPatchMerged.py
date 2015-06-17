# MeridionalPatchMerged.py
# Copyright (c) 2015 Peter Hokanson
# Vertical Limit Labs

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

from MeridionalPatch import MeridionalPatch
from Splines import *

class MeridionalPatchMerged(MeridionalPatch):
	"""Meridional patch composed of other meridional patch shapes merged on the 
	meridional (inlet-to-outlet) axis.  Splits given patches evenly in meridional
	coordinate. Possibly useful for combined pump-inducer shapes."""
	
	def __init__(self, patch_list):
		"""Construct a meridional patch merged from a list of patches."""
		self.patch_list = patch_list
		self.patch_count = len(patch_list)
		
		# NOTE: Not currently asserting that the patches actually align or anything.
		
	def __call__(self, m, s):
		
		# Figure out which subpatch to call with what meridional parameter
		patch_idx = int(math.floor(m * self.patch_count))
		patch_m = m * self.patch_count - patch_idx
		patch = self.patch_list[patch_idx]
		return patch(patch_m, s)
