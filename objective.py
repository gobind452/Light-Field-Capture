class Objective: # Simulate the effect of an objective lens
	def __init__(self,NA,flen):
		self.NA = NA
		self.flen = flen

	def distanceMagnification(self,u):
		f_inverse = 1/float(self.flen)
		u_inverse = 1/float(u)
		v_inverse = f_inverse - u_inverse
		v = 1/v_inverse
		return (round(v,5),round(v*u_inverse,5))

	def blurPsf(self,image,index): # Simulate the diffraction effects (To be implemented)
		return 0

	def propagateImage(self,image,u):
		return 0

