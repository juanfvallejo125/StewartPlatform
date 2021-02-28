from vpython import *
import numpy as np
from scipy.spatial.transform import Rotation as R
from math import *

class Simulation:
	def __init__(self, model): 
		self.r = model.r
		self.s = model.s
		self.origin = vec(0,0,0)
		self.ek = np.zeros((6,3))# Corners of the plate in base frame
		self.hk = np.zeros((6,3))
		self.kine = model
		scene.width=1920*0.8
		scene.height=1080*0.8
		self.generateVertices()
		self.drawPlate()
		self.drawServoAxisPoints()
		self.calculateHK()
		self.drawHorns()
		self.drawRods()
		ring(pos=self.origin, axis=vec(0,0,1), radius=self.s, thickness=0.1)
		ring(pos=self.origin, axis=vec(0,0,1), radius=self.r, thickness=0.1)

	def generateVertices(self):
		self.vertices = [vertex(pos=vec(self.kine.T[0],self.kine.T[1],self.kine.T[2]))]
		self.triangles = []
		self.ek = self.kine.Pk_solved
		for i,e in enumerate(self.ek):
			self.vertices.append(vertex(pos=vec(self.ek[i,0], self.ek[i,1], self.ek[i,2])))
		
	def drawPlate(self):
		# Plot the platform with triangles
		for i in range(len(self.vertices)-2):
			self.triangles.append(triangle(vs = [self.vertices[0], self.vertices[i+1], self.vertices[i+2]]))
		self.triangles.append(triangle(vs = [self.vertices[0], self.vertices[-1], self.vertices[1]]))
		# self.extrusionPath = [self.vertices[0].pos,self.triangles[0].v0.normal.norm()*self.thickness]
		# self.plate = extrusion(path=self.extrusionPath, shape=self.triangles[0], color=color.black)

	def updatePlate(self):
		self.ek = self.kine.Pk_solved
		self.vertices[0].pos = vec(self.kine.T[0],self.kine.T[1],self.kine.T[2])
		for i in range(len(self.ek)):
			self.vertices[i+1].pos = vec(self.ek[i,0], self.ek[i,1], self.ek[i,2])

	def drawServoAxisPoints(self):
		color_list = [color.black, color.blue, color.cyan, color.orange, color.green, color.magenta]
		for i,e in enumerate(self.kine.bk):
			sphere(pos=vec(self.kine.bk[i,0], self.kine.bk[i,1], self.kine.bk[i,2]), radius=2, color=color_list[i])

	def calculateHK(self):
		h_hat_x = np.reshape(np.cos(self.kine.alpha_k)*np.cos(self.kine.beta_k), (len(self.kine.alpha_k),1))
		# print("h_hat_x")
		# print(h_hat_x)
		h_hat_y = np.reshape(np.cos(self.kine.alpha_k)*np.sin(self.kine.beta_k), (len(self.kine.alpha_k),1))
		# print("h_hat_y")
		# print(h_hat_y)
		h_hat_z = np.reshape(np.sin(self.kine.alpha_k),(len(self.kine.alpha_k),1))
		# print("h_hat_z")
		# print(h_hat_z)
		h_hat = np.concatenate((h_hat_x, h_hat_y, h_hat_z), axis=1)
		# print("h_hat")
		# print(h_hat)
		h_vec = h_hat*self.kine.h
		self.hk = self.kine.bk + h_vec

	def drawHorns(self):
		self.horns = []
		for i in range(len(self.kine.bk)):
			self.horns.append(curve(vec(self.kine.bk[i,0], self.kine.bk[i,1], self.kine.bk[i,2]), 
								   vec(self.hk[i,0], self.hk[i,1], self.hk[i,2])))

	def updateHorns(self):
		for i in range(len(self.kine.bk)):
			self.horns[i].modify(0,vec(self.kine.bk[i,0], self.kine.bk[i,1], self.kine.bk[i,2]))
			self.horns[i].modify(1,vec(self.hk[i,0], self.hk[i,1], self.hk[i,2]))

	def drawRods(self):# Once I pass everything in as the kinematics object. Indexes should match in HK and ek
		self.rods = []
		for i in range(len(self.hk)):
			self.rods.append(curve(vec(self.hk[i,0], self.hk[i,1], self.hk[i,2]), 
								   vec(self.ek[i,0], self.ek[i,1], self.ek[i,2])))

	def updateRods(self):# Same as note above
		for i in range(len(self.hk)):
			self.rods[i].modify(0,vec(self.hk[i,0], self.hk[i,1], self.hk[i,2]))
			self.rods[i].modify(1,vec(self.ek[i,0], self.ek[i,1], self.ek[i,2]))

	def updatePlatform(self):
		#Update the platform
		self.updatePlate()
		self.calculateHK()
		self.updateHorns()
		self.updateRods()


class Kinematics:
	def __init__(self, 		   #Parameters needed for kinematics
				r: np.float64, #Radius of inscribed circle in inner triangle
				s: np.float64, #Radius of inscribed circle in outer triangle
				d: np.float64, #Rod length
				h: np.float64, #Servo Horn Length 
				):
		self.r = r
		self.s = s
		self.origin = vec(0,0,0)
		self.T0 = np.array([0,0,0])# Home position translation from base to plate origin
		self.pk = np.zeros((6,3))# Corners of the plate in plate's frame
		self.Pk_solved = np.zeros((6,3))# Corners of the plate in the base frame with desired pose
		self.bk = np.zeros((6,3))# Points of servo axis
		self.beta_k = np.zeros(6)# Angles of the servo horns in XY plane
		self.alpha_k = np.zeros(6)# Angles of the servo axis
		self.T = np.zeros(3)
		self.Rot = R.from_euler('x', 0, degrees=True)
		self.a = (2*self.r-self.s)/sqrt(3) #Side of the pyramidion
		self.d = d # Length of the rods
		self.h = h# Length of servo horns
		self.lk = np.array((6,3))# Vector from bk to ek in base frame
		self.calculateInitialPlatePoints()
		self.calculateServoAxisPoints()
		self.calculateBetaK()
		self.lineUpIndices()
		self.calculateZ0()
		self.solveIK(np.array([0,0,0]), self.Rot)

	def calculateInitialPlatePoints(self):
		for i,e in enumerate(self.pk):
			k = i+1
			ang = 2*pi/3*(k)/2
			self.pk[i,0] = self.s*cos(ang)+self.a/2*((-1)**k)*sin(ang)
			self.pk[i,1] = self.s*sin(ang)+self.a/2*((-1)**k)*-cos(ang)
			self.pk[i,2] = 0

	def calculateZ0(self):
		self.T0[2] = sqrt(d**2 + h**2 - (self.pk[0,0] - self.bk[0,0])**2 - (self.pk[0,1] - self.bk[0,1])**2)

	def calculateServoAxisPoints(self):
		for i,e in enumerate(self.bk):
			k = i+1
			ang = 2*pi/3*(k)/2 + pi/3
			self.bk[i,0] = (self.s+10)*cos(ang)+self.a/4*((-1)**k)*sin(ang)
			self.bk[i,1] = (self.s+10)*sin(ang)+self.a/4*((-1)**k)*-cos(ang)

	def calculateBetaK(self):
		i=0
		while i < len(self.bk):
			curr_idx = i % len(self.bk)
			next_idx = (i+1) % len(self.bk)
			v_side = self.pk[next_idx] - self.pk[curr_idx]
			self.beta_k[i] = atan2(v_side[1], v_side[0])+pi/3
			self.beta_k[i+1] = atan2(-v_side[1], -v_side[0])+pi/3
			i += 2
		# print("Initial Beta K")
		# print(self.beta_k)
		# print("Initial bk")
		# print(self.bk)
		# print("")

	def lineUpIndices(self):
		temp_bk = np.copy(self.bk)
		temp_beta_k = np.copy(self.beta_k)
		for i in range(len(self.bk)):
			idx = (i+1) % len(self.bk)
			self.bk[idx] = temp_bk[i]
			self.beta_k[idx] = temp_beta_k[i]
		# print("Corrected Indices Beta K")
		# print(self.beta_k)
		# print("Corrected Indices bk")
		# print(self.bk)
		# print("")


	def solveIK(self, T_rel, rot):
		T_prev = np.copy(self.T)
		Rot_prev = np.copy(self.Rot)
		Pk_solved_prev = np.copy(self.Pk_solved)
		self.T = self.T0 + T_rel
		self.Rot = rot
		self.lk = np.tile(self.T, (6,1)) + self.Rot.apply(self.pk) - self.bk
		self.Pk_solved = np.tile(self.T, (6,1)) + self.Rot.apply(self.pk)
		# print("Kinematics: lk")
		# print(self.lk)
		self.e_k = 2*self.h*self.lk[:,2] # shape is 6
		# print("Kinematics: e_k")
		# print(self.e_k)
		self.fk = 2*self.h*(np.cos(self.beta_k)*self.lk[:,0] + np.sin(self.beta_k)*self.lk[:,1])
		# print("Kinematics: f_k")
		# print(self.fk)
		self.gk = np.linalg.norm(self.lk, axis=1)**2 - (self.d**2 - self.h**2)
		# print("Kinematics: g_k")
		# print(self.gk)
		# print("Input to Arcsin")
		# print(self.gk/np.sqrt(self.e_k**2 + self.fk**2))
		try:
			self.alpha_k = np.arcsin(self.gk/np.sqrt(self.e_k**2 + self.fk**2)) - np.arctan2(self.fk, self.e_k)
			successful = True
		except:
			print("Kinematically unfeasible")
			self.T = T_prev
			self.Rot = Rot_prev
			self.Pk_solved = Pk_solved_prev
			successful = False

		return successful
		# print("Kinematics: Alpha K")
		# print(self.alpha_k)
		#Solve inverse kinematics, solves for alpha k

if __name__ == "__main__":
	ro = 60
	ri = 50
	d = 50
	h = 10
	np.seterr(all='raise')
	kine = Kinematics(ri, ro, d, h)
	sim = Simulation(kine)
	T = np.array([0,0,0])
	r = R.from_euler('x', 0, degrees=True)
	print("T0")
	print(kine.T0)
	print("T")
	print(T)
	print("R")
	print(r.as_euler('xyz'))
	x_rot=0
	y_rot=0
	z_rot=0
	x=0
	y=0
	z=0
	x_rot_prev=0
	y_rot_prev=0
	z_rot_prev=0
	x_prev=0
	y_prev=0
	z_prev=0
	dpos = 0.1
	drot = 0.1
	dlength = 0.1
	increment = True
	decrement = False
	while(True):
		rate(50)
		k = keysdown() # a list of keys that are down
		if 'left' in k: x -= dpos
		if 'right' in k: x += dpos
		if 'down' in k: y -= dpos
		if 'up' in k: y += dpos
		if '.' in k: z -= dpos
		if ';' in k: z += dpos
		if 'a' in k: x_rot -= drot
		if 'q' in k: x_rot += drot
		if 's' in k: y_rot -= drot
		if 'w' in k: y_rot += drot
		if 'd' in k: z_rot -= drot
		if 'e' in k: z_rot += drot

		if 't' in k: kine.d += dlength
		if 'g' in k: kine.d -= dlength
		if 'h' in k: kine.h += dlength
		if 'y' in k: kine.h -= dlength


		T = np.array([x,y,z])
		r = R.from_euler('xyz', [x_rot,y_rot,z_rot], degrees=True)
		if 'z' in k:
			T = np.array([0,0,0])
			r = R.from_euler('xyz', [0,0,0], degrees=True)
			x = 0
			y = 0
			z = 0
			x_rot = 0
			y_rot = 0
			z_rot = 0
		success = kine.solveIK(T,r)
		if(not success):
			x = x_prev
			y = y_prev
			z = z_prev
			x_rot = x_rot_prev
			y_rot = y_rot_prev
			z_rot = z_rot_prev
		else:
			x_prev = x
			y_prev = y
			z_prev = z
			x_rot_prev = x_rot
			y_rot_prev = y_rot
			z_rot_prev = z_rot

		sim.updatePlatform()

		if('p' in k):
			print("Pose")
			print("T")
			print(T)
			print("R")
			print(r.as_euler('xyz'))
			print("Servo Angles")
			print(kine.alpha_k)
