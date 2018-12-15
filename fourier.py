import numpy as np 
import cv2
from matplotlib import pyplot as plt 

def defocusImage(image,z): # Defocussing algorithm 
	otf = np.zeros([image.shape[0],image.shape[1]])
	for x in range(image.shape[0]):
		for y in range(image.shape[1]):
			otf[x][y] = np.exp((-z*z)*((x-image.shape[0]/2)*(x-image.shape[0]/2)+(y-image.shape[1]/2)*(y-image.shape[1]/2))) ## Need a more accurate OTF
	fft = np.fft.fftshift(np.fft.fft2(image))
	fft = np.multiply(fft,otf)
	fft = np.fft.ifft2(fft)
	fft = np.abs(fft)
	return fft

# Different file created for exploring other algorithms in Fourier domain (If of no use, will be removed)


