from matplotlib import pyplot as plt 
import numpy as np 
import cv2

def createImage(image,image_index):
	new_image = np.zeros([image.shape[0],image.shape[1]])
	x_start = 100 + (image_index-1)*10
	y_start = 100 - (image_index-1)*10
	y_end = 100 + (image_index-1)*10
	for x in range(x_start,x_start+10):
		for y in range(y_start-10,y_end+10):
			new_image[x][y] = image[x][y]
	x_start = 100 - (image_index-1)*10
	for x in range(x_start-10,x_start):
		for y in range(y_start-10,y_end+10):
			new_image[x][y] = image[x][y]
	x_end = 100 + (image_index-1)*10
	for y in range(y_start-10,y_start):
		for x in range(x_start,x_end):
			new_image[x][y] = image[x][y]
	for y in range(y_end,y_end+10):
		for x in range(x_start,x_end):
			new_image[x][y] = image[x][y]
	return new_image

image = cv2.imread("test.png",0)

for index in range(10):
	section = createImage(image,index) # Creates an section according to some rule (can be changed in stack.py)
	cv2.imwrite("Stack/image"+str(index)+".png",section)



