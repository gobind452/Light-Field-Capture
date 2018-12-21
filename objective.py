import numpy as np
import cv2

class Objective: # Simulate the effect of an objective lens
    def __init__(self,NA,flen):
        self.NA = NA
        self.flen = flen

    def distanceMagnification(self,u): # Simple lens equation
        f_inverse = 1/float(self.flen)
        u_inverse = 1/float(u)
        v_inverse = f_inverse - u_inverse
        v = 1/v_inverse
        return (round(v,5),round(v*u_inverse,5))

    def blurPsf(self,image,index): # Simulate the diffraction effects (To be implemented)
        return 0

    def propagateImage(self,image,u):
        return 0

    def imageFormation(self,image, u, x_final, y_final):
        abs_mag = round(self.flen/(u - self.flen), 5) # absolute value of magnification
        x_image, y_image = image.shape[0], image.shape[1]
        X = int(x_image*abs_mag) #size in pixels of the image formed by objective
        Y = int(y_image*abs_mag)
        IMG = np.zeros((X, Y))

        for x in range(x_image):
            for y in range(y_image):
                IMG[int(X/2 - abs_mag*(x - x_image/2)) - 1][int(Y/2 - abs_mag*(y - y_image/2)) - 1] = image[x][y]

        res =  cv2.resize(IMG, (x_final, y_final), interpolation = cv2.INTER_AREA) #resizing IMG to fit the desired size
	#interpolation chosen assuming that resize will correspond to shrinking
        return res.astype(np.uint8)
