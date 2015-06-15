# FreeVortex.py
# Copyright (c) 2014, 2015 Peter Hokanson
# Vertical Limit Labs

import os, sys, shutil
import math
import csv
import io
from subprocess import call, check_call, check_output

import numpy as np
import scipy
try:
	import matplotlib.pyplot as plt
	MATPLOTLIB_DISABLED = False
except:
	MATPLOTLIB_DISABLED = True

import PyFoam as pf
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


from Splines import *

def loadPatchVectorSamples(sampleFile):
	with open(sampleFile) as csvfile:
		datareader = csv.reader(csvfile, delimiter=' ')
		next(datareader)
		next(datareader)
		
		data = [row for row in datareader]
		points = np.array([[float(r[0]), float(r[1]), float(r[2])] for r in data])
		field_vals = np.array([[float(r[3]), float(r[4]), float(r[5])] for r in data])
		return points, field_vals

class FreeVortex(object):
	"""Representation of a free-vortex region of flow.  This is used as a base
	for both rotor and stator segments."""
	def __init__(self, 
	              casename="cases/freevortex",
	              compressible=False,
	              rho=1e3,
	              inlet_s0=np.array([3e-3, 7e-3]),
	              inlet_s1=np.array([7.8e-3, 7e-3]),
	              inlet_v=np.array([0.0, 0.0, -39.6]),
	              outlet_s0=np.array([12.8e-3, 0.0]),
	              outlet_s1=np.array([12.8e-3, 2.0e-3]),
	              outlet_v=np.array([39.63, -19.15, 0.0]),
	              points_m = 40,
	              points_s = 20):
		"""Create a representation of free-vortex flow through a region.
		
		Keyword arguments:
		casename -- directory to create for OpenFOAM files
		rho -- fluid density in kg/m^3
		inlet_s0 -- numpy array specifying (r, z) coordinates at inner edge of inlet
		inlet_s1 -- numpy array specifying (r, z) coordinates at outer edge of inlet
		inlet_v -- numpy array specifying (r, th, z) velocity (uniform) at inlet
		outlet_s0 -- numpy array specifying (r, z) coordinates at inner edge of outlet
		outlet_s1 -- numpy array specifying (r, z) coordinates at outer edge of outlet
		outlet_v -- numpy array specifying (r, th, z) velocity (uniform) at outlet
		points_m -- number of vertices in the meridional direction (inlet to outlet)
		points_s -- number of vertices in the shroud direction (hub to shroud)
		"""
		
		# Case directory
		self.casename = casename
		self.case_template = "case_templates/freevortex"
		
		# Simulation properties
		self.points_m = points_m
		self.points_s = points_s
		
		# Flow properties
		self.rho = rho
		self.compressible = compressible
		assert not self.compressible, "Compressible flow not yet supported"
		
		self.inlet_s0 = inlet_s0
		self.inlet_s1 = inlet_s1
		self.inlet_v = inlet_v
		self.outlet_s1 = outlet_s1
		self.outlet_s0 = outlet_s0
		self.outlet_v = outlet_v
		
		# set up folder structure
		self.makeOFCase()
		self.makeMeridionalPatch()
		self.makeOFMesh()
		self.setOFBoundaries()
		self.solve()
	
	def makeOFCase(self):
		"""Remove target case directory and copy OpenFOAM case template."""
		print("Copying OpenFOAM case %s from template %s" % 
		       (self.casename, self.case_template))
		shutil.rmtree(self.casename)
		shutil.copytree(self.case_template, self.casename)
	
	def makeMeridionalPatch(self):
		"""Produces a Bezier patch connecting the inlet and the outlet.  This is
		second order, with three control points.  The center control point is chosen
		such that inlet and outlet velocities in the meridional plane match the
		stated values.
		
		NOTE: This will likely break for a purely axial flowfield, and will need
		to go to fourth order somehow."""
		k_11 = self.inlet_s0
		k_12 = self.inlet_s1
		
		k_31 = self.outlet_s0
		k_32 = self.outlet_s1
		
		# Let's find the intersection:
		# Our velocities are 3D, so we'll need to strip out the tangential component.
		v_in_rz = np.array([self.inlet_v[0], self.inlet_v[2]]).T
		v_out_rz = np.array([self.outlet_v[0], self.outlet_v[2]]).T
		k_21 = intersection_2d(k_11, k_11+v_in_rz, k_31, k_31+v_out_rz)
		k_22 = intersection_2d(k_12, k_12+v_in_rz, k_32, k_32+v_out_rz)
		
		self.b = BezierSurface(np.array([[k_11, k_12],
		                            [k_21, k_22],
		                            [k_31, k_32]]))
		
		# Grid points
		self.r = np.zeros((self.points_m, self.points_s))
		self.z = np.zeros((self.points_m, self.points_s))
		for m in range(self.points_m):
			for s in range(self.points_s):
				pt_rz = self.b(m/(self.points_m-1), s/(self.points_s-1))
				self.r[m][s] = pt_rz[0]
				self.z[m][s] = pt_rz[1]
		return (self.r, self.z)
	
	def makeOFMesh(self, runBlockMesh=True):
		"""Update the blockMeshDict file in the OpenFOAM case to represent our new
		mesh.  Optionally runs blockMesh."""
		
		# The points are numbered in a single array, front and back.  These
		# functions index the different sides' points.
		idx1 = lambda m, s: s * self.points_m + m
		idx2 = lambda m, s: self.points_m * self.points_s + s * self.points_m + m
		
		vertices = {}
		for m in range(self.points_m):
			for s in range(self.points_s):
				vertex = [self.r[m][s], self.z[m][s], self.r[m][s] * 0.01]
				vertices[idx1(m, s)] = vertex
				#vertex[2] *= -1
				vertices[idx2(m, s)] = [vertex[0], vertex[1], -vertex[2]]
		
		blocks = []
		for m in range(self.points_m-1):
			for s in range(self.points_s-1):
				corners = [idx2(m,s), idx2(m+1,s), idx2(m+1,s+1), idx2(m,s+1),
				           idx1(m,s), idx1(m+1,s), idx1(m+1,s+1), idx1(m,s+1)]
				cells = [1,1,1]
				simpleGrading = [1,1,1]
				blocks.append({"corners" : corners, "cells" : cells, "simpleGrading" : simpleGrading})
		
		edges = {}
		# We don't do anything interesting with the edges in this case
		
		# boundary faces are listed CW from within the block
		front = {"type":"wedge"}
		faces = []
		for m in range(self.points_m-1):
			for s in range(self.points_s-1):
				faces.append([idx1(m,s), idx1(m+1,s), idx1(m+1,s+1), idx1(m,s+1)])
		front["faces"] = faces
		
		back = {"type":"wedge"}
		faces = []
		for m in range(self.points_m-1):
			for s in range(self.points_s-1):
				faces.append([idx2(m,s), idx2(m,s+1), idx2(m+1,s+1), idx2(m+1,s)])
		back["faces"] = faces
		
		inlet = {"type":"patch"}
		faces = []
		m = 0
		for s in range(self.points_s-1):
			faces.append([idx1(m,s), idx1(m,s+1), idx2(m,s+1), idx2(m,s)])
		inlet["faces"] = faces
		
		outlet = {"type":"patch"}
		faces = []
		m = self.points_m-1
		for s in range(self.points_s-1):
			faces.append([idx1(m,s), idx2(m,s), idx2(m,s+1), idx1(m,s+1)])
		outlet["faces"] = faces
		
		wallShroud = {"type":"wall"}
		faces = []
		s = self.points_s-1
		for m in range(self.points_m-1):
			faces.append([idx1(m,s), idx1(m+1,s), idx2(m+1,s), idx2(m,s)])
		wallShroud["faces"] = faces
		
		wallHub = {"type":"wall"}
		faces = []
		s = 0
		for m in range(self.points_m-1):
			faces.append([idx1(m,s), idx2(m,s), idx2(m+1,s), idx1(m+1,s)])
		wallHub["faces"] = faces
		
		boundary = {"front":front, 
		            "back":back, 
		            "inlet":inlet, 
		            "wallShroud":wallShroud, 
		            "wallHub":wallHub,
		            "outlet":outlet}
		
		self.blockmesh_data = {
			"vertices" : vertices,
			"blocks" : blocks,
			"edges" : edges,
			"boundary" : boundary
		}
		self.writeOFMesh(runBlockMesh)
		
	def writeOFMesh(self, runBlockMesh=True):
		"""Write OpenFOAM mesh data to the blockMeshDict file.  Optionally run blockMesh."""
		bm = ParsedParameterFile(os.path.join(self.casename, 
		                                      "constant/polyMesh/blockMeshDict"),
		                         longListOutputThreshold=0)
		
		#TODO: Fix PyFoam bug that writes incorrect list length
		
		vertex_keys = sorted(self.blockmesh_data["vertices"].keys(), key=lambda x: int(x))
		vertices = [self.blockmesh_data["vertices"][key] for key in vertex_keys]
		bm["vertices"] = vertices
		
		blocks = []
		for b in self.blockmesh_data["blocks"]:
			#print(b)
			blocks += ["hex"]
			blocks += [b["corners"]]
			blocks += [b["cells"]]
			blocks += ["simpleGrading"]
			blocks += [b["simpleGrading"]]
		#print(blocks)
		bm["blocks"] = blocks
		
		edges = []
		for e in self.blockmesh_data["edges"]:
			pass # TODO: implement when necessary
		
		boundaries = []
		#print(bm["boundary"])
		for key in self.blockmesh_data["boundary"].keys():
			boundaries.append(key)
			boundaries.append(self.blockmesh_data["boundary"][key])
		#print(boundaries)
		bm["boundary"] = boundaries
		
		print(bm.content)
		bm.writeFile()
		
		if runBlockMesh:
			check_call(["blockMesh"], cwd=self.casename)
			check_call(["checkMesh"], cwd=self.casename)
		
	def setOFBoundaries(self):
		"""Set up boundary conditions"""
		self.boundaries = {
			"U": {
				"inlet": {
					"type": "fixedValue",
					"value": "uniform (%f %f %f)" % (self.inlet_v[0], 
					                                 self.inlet_v[1],
					                                 self.inlet_v[2])},
				"outlet": {
					"type": "fixedValue",
					"value": "uniform (%f %f %f)" % (self.outlet_v[0], 
					                                 self.outlet_v[1],
					                                 self.outlet_v[2])},
				"wallShroud": {"type":"slip"},
				"wallHub":{"type":"slip"},
				"front":{"type":"wedge"},
				"back":{"type":"wedge"}
				},
			"p": {}
		}
		self.writeOFBoundaries()
	
	def writeOFBoundaries(self):
		"""Write boundary conditions to 0/<field> file in case directory."""
		for field in self.boundaries:
			f = ParsedParameterFile(os.path.join(self.casename, "0/" + field))
			for boundary in self.boundaries[field]:
				f["boundaryField"][boundary] = self.boundaries[field][boundary]
			f.writeFile()
	
	def solve(self):
		"""Call OpenFOAM solver for case, then read back solved data and convert
		it to cylindrical coordinates."""
		check_call(["simpleFoam"], cwd=self.casename)
		
		# Get velocity figures at grid points:
		check_call(["sample", "-latestTime"], cwd=self.casename)
		end_time = check_output(["foamListTimes", "-latestTime"], cwd=self.casename)
		end_time = end_time.decode('utf-8')[:-1]
		xyz_points, u_xyz_points = loadPatchVectorSamples(os.path.join(self.casename, "postProcessing/surfaces", end_time, "U_frontWall.raw"))
		
		self.xyz_points = xyz_points
		self.u_xyz_points = u_xyz_points
		
		self.rz_points = np.array([[np.sqrt(xyz[0]**2 + xyz[2]**2), xyz[1]] for xyz in xyz_points])
		self.th_points = np.array([math.atan2(xyz[2], xyz[0]) for xyz in xyz_points])
		self.u_rtz_points = np.array([[math.cos(th) * xyz[0] + math.sin(th) * xyz[2], 
		                               -math.sin(th) * xyz[0] + math.cos(th) * xyz[2], 
		                               xyz[1]] for (xyz, th) in zip(u_xyz_points, self.th_points)])


if __name__=="__main__":
	fv = FreeVortex()
	fv.solve()
	