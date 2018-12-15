import objective,stack,fourier
import numpy as np 
import cv2
from matplotlib import pyplot as plt 

# All distances in meters

f_objective = 0.01 # Objective focal lens
specimen_distance = 0.0105 # Specimen to objective distance
NA = 0.4 # Numerical aperture 
d_objective = 0 # Diameter of objective lens
field_of_view_x = 0 # Field of view in the image plane (Magnification*Actual specimen size)
field_of_view_y = 0 

lens = objective.Objective(NA,f_objective) # Lens object (implemented in objective.py)
array_distance,magnification = lens.distanceMagnification(specimen_distance) # Distance for the microlens_array such that there is a one to one mapping

microlens_pitch = 0 #Microlens pitch
f_microlens = 0 # Microlens focal length
pixel_ratio = 0 # Number of sensor_pixels per m
x_lenslet = 0 # Number of lenslets in x
y_lenslet = 0 # Number of lenslets in y
pixels = 10 #Pixels behind each microlens (in 1D) 

''' Points to note
1. Pixel ratio should be fixed, as the sensor is a given, rather than pixels behind each microlens as it depends on the microlens diameter. However for convenience calibrating the other way round is better.
2. The microlens are assumed to be spherical in nature, so pixels behind each microlens can be given as 10 in 1D (Actual pixels would be 100 in 2D)
3. The actual fixed quantities before an experiment would only be the specimen dimensions. Then one would choose an objective, which fixes the objective focal
length, diameter and NA. This then fixes the field_of_view, which would then help to determine the number of lenslets to include in the array and where to place the array. This would 
ultimately fix the pixels behind each microlens.
'''
## All 0 initialised values will be calibrated later

# Create initial stack
image = cv2.imread("test.png",0) # 2D image for creating the stack

# Assumption - The number of pixels in the image = number of lenslets (Allows modelling of the image as collection of sources)
x_lenslet = image.shape[0] 
y_lenslet = image.shape[1]


## Uncomment for creating a stack

#for index in range(10):  # 10 is the number of slices
#	section = stack.createImage(image,index) # Creates an section according to some rule (can be changed in stack.py)
#	cv2.imwrite("Stack/image"+str(index)+".png",section)


## End 

# Calibration of the experimental details
def matching():
	global d_objective  
	global field_of_view_x
	global field_of_view_y
	global pixel_ratio 
	global microlens_pitch 
	global f_microlens 
	d_objective = NA*2*f_objective # Calculated using basic formulae for NA
	field_of_view_x = field_of_view_y =  d_objective # Assumption (Can be set)
	pixel_ratio = x_lenslet*pixels/field_of_view_x #Pixels per m  
	microlens_pitch = field_of_view_x/x_lenslet # Microlens array covers the entire image plane field of view
	f_microlens = microlens_pitch*array_distance*2/d_objective 
	# Calculated by using the fact that the cone of rays from the objective should converge on the image plane and diverge to cover only the pixels behind each microlens

## Convert into image stack
 
## Nothing here 
## Would contain the code for modelling the effects of the objective lens (objective.py)

## End

z_spacing = 0.000002 # Stack spacing ( Assumption ) (Lower stack spacing allows faster processing)

sensor = np.zeros([x_lenslet*pixels,y_lenslet*pixels]) # Initialise the sensor 
matching() # Calibrate

def findIntersection(x1,y1,z1,x2,y2,z2,z_plane): # Helper function for finding intersection of line passing through (x1,y1,z1) and (x2,y2,z2) with the plane z = z_plane
	lamda = float(z_plane-z1)/(z1-z2)
	return(int(round(x1+lamda*(x1-x2),0)),int(round(y1+lamda*(y1-y2),0)))

def calculateRadius(image_distance,distance_to_array): # Calculates the radius of the cone of rays incident on the sensor (In pixels)
	return int(round(d_objective*(distance_to_array+f_microlens)*pixel_ratio/float(2*image_distance),0))

def sensorMapping(x_lens,y_lens,x,y,x1,y1,distance_to_array): # Given a microlens coordinate and the angle of the ray incident, gets the sensor pixel to be activated
	x_ratio = float(x-x1)/(distance_to_array+f_microlens)
	y_ratio = float(y-y1)/(distance_to_array+f_microlens)
	return (int(x_lens*pixels+pixels/2+x_ratio*f_microlens),int(y_lens*pixels+pixels/2+y_ratio*f_microlens))

def realCoord(x,y): # Helper function 
 	return (x*pixels+float(pixels-1)/2,y*pixels+float(pixels-1)/2)

def propagateToSensor(image,image_distance): # Main function for getting the images
	distance_to_array = round(array_distance - image_distance,5) # Distance to the microlens array
	radius = calculateRadius(image_distance,distance_to_array) # Sensors activated 
	z_add = (f_microlens+distance_to_array)*pixel_ratio # Used for calculating the cosine
	normalise = round(1/float(radius*radius),5) # Intensity has to be conserved
	print(radius,distance_to_array,normalise)
	for x in range(x_lenslet):
		for y in range(y_lenslet):
			if image[x][y] == 0: 
				continue
			x_real,y_real = realCoord(x,y) # Get the corresponding sensor coordinates 
			center_x,center_y = findIntersection(x_real,y_real,distance_to_array,x_lenslet*pixels/2,y_lenslet*pixels/2,array_distance,-f_microlens)
			# Get the center of the cone of rays (Formed by extending the line between image point and center of objective lens)
			for x1 in range(center_x-radius,center_x+radius):
				extent = int(np.sqrt(radius*radius-(x1-center_x)*(x1-center_x)))
				for y1 in range(center_y-extent,center_y+extent): # For sensor pixels within that circle 
					lens_x,lens_y = findIntersection(x1,y1,-f_microlens,x_real,y_real,distance_to_array,0) # Find which lens it passes through
					lens_x = int((lens_x-pixels/2)/pixels)
					lens_y = int((lens_y-pixels/2)/pixels)
					sensor_x,sensor_y = sensorMapping(lens_x,lens_y,x_real,y_real,x1,y1,distance_to_array) # Get the sensor pixel it would hit
					cosine = np.power(float(z_add)/((x_real-x1)*(x_real-x1)+(y_real-y1)*(y_real-y1)+z_add*z_add),1.5)
					sensor[sensor_x][sensor_y] = sensor[sensor_x][sensor_y] + image[x][y]*cosine*normalise # Add the normalised intensity

# Driving loop
for index in range(2,5): # 2,5 arbitrary. Can range between the number of slices
	section = cv2.imread("Stack/image"+str(index)+".png",0) 
	distance = specimen_distance +(index-2)*z_spacing # 2 here determines that slice in image2 is in focus
	image_distance,mag = lens.distanceMagnification(distance)
	print(image_distance,mag)
	propagateToSensor(section,image_distance)

image = cv2.normalize(sensor,image, 0, 255, cv2.NORM_MINMAX) # Normalise the entire image
cv2.imwrite("sensor.png",sensor)



