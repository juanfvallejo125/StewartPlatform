from vpython import *
import numpy as np
from scipy.spatial.transform import Rotation as R
from math import *

class Simulation:
	def __init__(self, 		   #Parameters needed to make the simulation. Probalby pass in kinematics
				r: np.float64, #Radius of inscribed circle in inner triangle
				s: np.float64, #Radius of inscribed circle in outer triangle
				): 
		self.r = r
		self.s = s
		self.origin = vec(0,0,0)
		self.ek = np.zeros((6,3))# Corners of the plate initial position
		self.bk = np.zeros((6,3))# Points of servo axis
		self.hk = np.zeros((6,3))
		self.beta_k = np.zeros(6)# Angles of the servo horns in XY plane
		self.alpha_k = np.zeros(6)# Angles of the servo axis
		self.z0 = 50 # Base height of plate
		self.a = (2*self.r-self.s)/sqrt(3) #Side of the pyramidion
		self.d = 70 # Length of the rods
		self.h = self.a/2# Length of servo horns
		self.calculateInitialPlatePoints()
		self.drawPlate()
		self.calculateServoAxisPoints()
		self.calculateBetaK()
		self.calculateHK()
		self.drawHorns()
		self.drawRods()
		ring(pos=self.origin, axis=vec(0,0,1), radius=s, thickness=0.1)
		ring(pos=self.origin, axis=vec(0,0,1), radius=r, thickness=0.1)

	def calculateInitialPlatePoints(self):
		self.vertices = [vertex(pos=vec(0,0,self.z0))]
		self.triangles = []
		for i,e in enumerate(self.ek):
			k = i+1
			ang = 2*pi/3*(k)/2
			self.ek[i,0] = self.s*cos(ang)+self.a/2*((-1)**k)*sin(ang)
			self.ek[i,1] = self.s*sin(ang)+self.a/2*((-1)**k)*-cos(ang)
			self.ek[i,2] = self.z0
			self.vertices.append(vertex(pos=vec(self.ek[i,0], self.ek[i,1], self.ek[i,2])))
		
	def drawPlate(self):
		# Plot the platform with triangles
		for i in range(len(self.vertices)-2):
			self.triangles.append(triangle(vs = [self.vertices[0], self.vertices[i+1], self.vertices[i+2]]))
		self.triangles.append(triangle(vs = [self.vertices[0], self.vertices[-1], self.vertices[1]]))
		# self.extrusionPath = [self.vertices[0].pos,self.triangles[0].v0.normal.norm()*self.thickness]
		# self.plate = extrusion(path=self.extrusionPath, shape=self.triangles[0], color=color.black)

	def updatePlate(self):
		for i in range(len(self.vertices)):
			self.vertices[i].pos = vec(self.ek[i,0], self.ek[i,1], self.ek[i,2])

	def calculateServoAxisPoints(self):
		color_list = [color.black, color.blue, color.cyan, color.orange, color.green, color.magenta]
		for i,e in enumerate(self.bk):
			k = i+1
			ang = 2*pi/3*(k)/2 + pi/3
			self.bk[i,0] = (self.s+10)*cos(ang)+self.a/4*((-1)**k)*sin(ang)
			self.bk[i,1] = (self.s+10)*sin(ang)+self.a/4*((-1)**k)*-cos(ang)
			sphere(pos=vec(self.bk[i,0], self.bk[i,1], self.bk[i,2]), radius=2, color=color_list[i])

	def calculateBetaK(self):
		i=0
		while i < len(self.ek):
			curr_idx = i % len(self.bk)
			next_idx = (i+1) % len(self.bk)
			v_side = self.ek[next_idx] - self.ek[curr_idx]
			self.beta_k[i] = atan2(v_side[1], v_side[0])+pi/3
			self.beta_k[i+1] = atan2(-v_side[1], -v_side[0])+pi/3
			print(self.beta_k[i]*180/pi)
			print(self.beta_k[i+1]*180/pi)
			i += 2

	def calculateHK(self):
		for i in range(len(self.alpha_k)):
			h_hat = np.array([cos(self.alpha_k[i])*cos(self.beta_k[i]),
						      cos(self.alpha_k[i])*sin(self.beta_k[i]), 
						      sin(self.alpha_k[i])])
			h_vec = h_hat*self.h
			self.hk[i] = self.bk[i] + h_vec

	def drawHorns(self):
		self.horns = []
		for i in range(len(self.bk)):
			self.horns.append(curve(vec(self.bk[i,0], self.bk[i,1], self.bk[i,2]), 
								   vec(self.hk[i,0], self.hk[i,1], self.hk[i,2])))

	def updateHorns(self):
		for i in range(len(self.bk)):
			self.horns[i].pos = [vec(self.bk[i,0], self.bk[i,1], self.bk[i,2]), 
								vec(self.hk[i,0], self.hk[i,1], self.hk[i,2])]

	def drawRods(self):
		self.rods = []
		for i in range(len(self.hk)):
			next_idx = (i+1) % len(self.hk)
			self.rods.append(curve(vec(self.hk[i,0], self.hk[i,1], self.hk[i,2]), 
								   vec(self.ek[next_idx,0], self.ek[next_idx,1], self.ek[next_idx,2])))

	def updateRods(self):
		for i in range(len(self.hk)):
			next_idx = (i+1) % len(self.hk)
			self.rods[i].pos = [vec(self.hk[i,0], self.hk[i,1], self.hk[i,2]), 
								vec(self.ek[next_idx,0], self.ek[next_idx,1], self.ek[next_idx,2])]

	def updatePlatform(self, new_alpha_k, new_ek):
		#Update the platform
		self.alpha_k = new_alpha_k
		self.ek = new_ek
		self.calculateHK()
		self.updatePlate()
		self.updateHorns()
		self.updateRods()


class Kinematics:
	def __init__(self, 		   #Parameters needed for kinematics
				r: np.float64, #Radius of inscribed circle in inner triangle
				s: np.float64, #Radius of inscribed circle in outer triangle
				d: np.float64,
				h: np.float64,

				):
		self.r = r
		self.s = s
		self.origin = vec(0,0,0)
		self.pk = np.zeros((6,3))# Corners of the plate in plates frame
		self.bk = np.zeros((6,3))# Points of servo axis
		self.beta_k = np.zeros(6)# Angles of the servo horns in XY plane
		self.alpha_k = np.zeros(6)# Angles of the servo axis
		self.z0 = 0 # Base height of plate, needs to be calculated
		self.T0 = np.array([0,0,self.z0])
		self.T = np.zeros(3)
		self.a = (2*self.r-self.s)/sqrt(3) #Side of the pyramidion
		self.d = 70 # Length of the rods
		self.h = self.a/3# Length of servo horns
		self.calculateInitialPlatePoints()
		self.drawPlate()
		self.calculateServoAxisPoints()
		self.calculateBetaK()
		self.calculateHK()
		self.drawHorns()
		self.drawRods()

	def calculateInitialPlatePoints(self):
		for i,e in enumerate(self.pk):
			k = i+1
			ang = 2*pi/3*(k)/2
			self.pk[i,0] = self.s*cos(ang)+self.a/2*((-1)**k)*sin(ang)
			self.pk[i,1] = self.s*sin(ang)+self.a/2*((-1)**k)*-cos(ang)
			self.pk[i,2] = 0

	def solveIK(T_rel, R):
		ik=1
		#Solve inverse kinematics, solves for alpha k

if __name__ == "__main__":
	sim = Simulation(50,60)
	# print(sim.ek)
	#Do something