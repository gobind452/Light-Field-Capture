import objective,stack
import numpy as np 
import cv2
from matplotlib import pyplot as plt 
import json

# All distances in meters

para = {}

def readParameters(para): # Reads the experimental data from a json file
	with open('data.txt') as infile:
		temp = json.load(infile)
	for key in temp.keys():
		para[key] = temp[key]

readParameters(para)

'''Parameters dictionary
1. specimen_distance = Distance of specimen to objective
2. specimen_x,specimen_y = Dimensions of specimen 
3. NA,f_objective = Numerical aperture and focal length of objective 
4. x_lenslet,y_lenslet = Number of lenslets in the microlens array 
5. x_sensor,y_sensor = Dimensions of sensor 
6. array_spacing = Radius of a microlens in microlens array 
7. z_spacing = Spacing between two different sections in the stack 
8. x_tilt,y_tilt = cosine of the tilt of the array plane about origin in x and y direction (1,1 corresponds to no tilt)
'''

lens = objective.Objective(para['NA'],para['f_objective']) 

specimen = cv2.imread("specimen.png",0)

x_matrix = np.zeros([3,3]) #Rotation about x_axis
y_matrix = np.zeros([3,3]) # Rotation about y_axis
x_matrix[0][0] = x_matrix[1][1] = x_matrix[2][2] = y_matrix[0][0] = y_matrix[1][1] = y_matrix[2][2] = 1 #Identity Matrix

transform = np.zeros([3,3]) # Matrix to transform from non tilted plane to tilted plane coordinates
inverse_transform = np.zeros([3,3]) # Inverse


# Calibration of the experimental details
def calibration(para):
	para['d_objective'] = para['NA']*2*para['f_objective'] # Calculated using basic formulae for NA
	para['image_distance'],para['magnification'] = lens.distanceMagnification(para['specimen_distance']) 
	para['x_image'] = para['magnification']*para['specimen_x'] # Image height = Magnification*Object Height
	para['y_image'] = para['magnification']*para['specimen_y']
	para['image_spacing'] = para['x_image']/specimen.shape[0] # Spacing between the pixels of the image
	para['x_array'] = para['x_lenslet']*para['array_spacing'] # Dimensions of microlens array
	para['y_array'] = para['y_lenslet']*para['array_spacing'] 
	para['sensor_xpixels'] = int(para['x_sensor']/para['sensor_spacing']) # Number of sensors in x,y
	para['sensor_ypixels'] = int(para['y_sensor']/para['sensor_spacing'])
	if para['x_tilt'] is not 1: # Create rotation matrix 
		global x_matrix
		x_matrix[1][1] = para['x_tilt'] #cosine
		x_matrix[1][2] = np.sqrt(1-para['x_tilt']*para['x_tilt']) #sin
		x_matrix[2][1] = -1*np.sqrt(1-para['x_tilt']*para['x_tilt'])
		x_matrix[2][2] = para['x_tilt']
	if para['y_tilt'] is not 1:
		global y_matrix
		y_matrix[0][0] = para['y_tilt']
		y_matrix[0][2] = np.sqrt(1-para['y_tilt']*para['y_tilt'])
		y_matrix[2][0] = -1*np.sqrt(1-para['y_tilt']*para['y_tilt'])
		y_matrix[2][2] = para['y_tilt']
	global transform 
	transform = np.matmul(y_matrix,x_matrix)
	global inverse_transform
	inverse_transform = np.matmul(np.transpose(x_matrix),np.transpose(y_matrix))
	vector = [0,0,1] # Normal to the plane in the plane coordinates
	vector = np.matmul(inverse_transform,vector) # Convert to coordinates in xyz plane
	para['a'] = vector[0] #ax+by+cz = 0 (Equation of microlens plane (Assumed passes through origin)
	para['b'] = vector[1] #[a,b,c] is normal to the plane with the equation ax+by+cz = 0
	para['c'] = vector[2]
	para['imagex_shift'] = para['imagey_shift'] = 0
	para['sensorx_shift'] = para['sensory_shift'] = 0


calibration(para)


## Uncomment for creating a stack

#for index in range(10):  # 10 is the number of slices
#	section = stack.createImage(image,index) # Creates an section according to some rule (can be changed in stack.py)
#	cv2.imwrite("Stack/image"+str(index)+".png",section)


## End 


## Convert into image stack
 
## Nothing here 
## Would contain the code for modelling the effects of the objective lens (objective.py)

## End

para['z_spacing'] = 0.000002 # Stack spacing ( Assumption ) (Lower stack spacing allows faster processing)
epsilon = 0.01 #Approximation


sensor = np.zeros([para['sensor_xpixels'],para['sensor_ypixels']]) # Initialise the sensor 

def findIntersection(x1,y1,z1,x2,y2,z2,z_plane): # Helper function for finding intersection of line passing through (x1,y1,z1) and (x2,y2,z2) with the plane z = z_plane
	lamda = float(z_plane-z1)/(z1-z2)
	return(x1+lamda*(x1-x2),y1+lamda*(y1-y2))

def findIntersectionAngle(x1,y1,z1,x2,y2,z2): # Used for calculated intersection with the plane ax+by+cz = 0
	lamda = -(x1+y1+z1)/(para['a']*(x1-x2)+para['b']*(y1-y2)+para['c']*(z1-z2))
	return (x1+lamda*(x1-x2),y1+lamda*(y1-y2),z1+lamda*(z1-z2))

def calculateRadius(image_distance,distance_to_array): # Calculates the radius of the cone of rays incident on the sensor (In pixels)
	return int(round(para['d_objective']*(distance_to_array+para['f_microlens'])/float(2*image_distance*para['sensor_spacing']),0))

def sensorMapping(x_lens,y_lens,x1,y1,x,y,distance_to_array): # Given a microlens coordinate and the angle of the ray incident, gets the sensor pixel to be activated
	lenslet_x,lenslet_y = realToPixel(x_lens,y_lens,'array')
	if lenslet_x >= para['x_lenslet'] or lenslet_x < 0: # If lens is outside the array
		return (x1,y1) # No deviation
	if lenslet_y >= para['y_lenslet'] or lenslet_y < 0:
		return (x1,y1) # No deviation
	lens_x,lens_y = pixelToReal(lenslet_x,lenslet_y,'array')
	if (lens_x-x_lens)*(lens_x-x_lens) + (lens_y-y_lens)*(lens_y-y_lens) > para['array_spacing']*para['array_spacing']/4:
		return (x1,y1) # Doesnt fall in the lens
	x_ratio = float(x1-x)/(distance_to_array+para['f_microlens'])
	y_ratio = float(y1-y)/(distance_to_array+para['f_microlens'])
	return (lens_x+x_ratio*para['f_microlens'],lens_y+y_ratio*para['f_microlens'])
	
def sensorMappingAngle(x_lens,y_lens,z_lens,x1,y1,x,y,distance_to_array): # For tilted arrays
	lenslet_x,lenslet_y = realToPixelAngle(x_lens,y_lens,z_lens,'array')
	if lenslet_x >= para['x_lenslet'] or lenslet_x < 0: # If lens is outside the array
		return (x1,y1)
	if lenslet_y >= para['y_lenslet'] or lenslet_y < 0:
		return (x1,y1)
	lens_x,lens_y,lens_z = pixelToRealAngle(lenslet_x,lenslet_y,'array')
	x_ratio = float(x1-x)/(distance_to_array+para['f_microlens'])
	y_ratio = float(y1-y)/(distance_to_array+para['f_microlens'])
	return (lens_x+x_ratio*(para['f_microlens']+lens_z),lens_y+y_ratio*(para['f_microlens']+lens_z))

def pixelToReal(x,y,flag): # Helper function (Maps pixels to the real coordinates in the setup). Flag here corresponds to which object we are talking about (Sensor, Microlesn array, Image) 
 	return(y*para[flag+'_spacing']-0.5*(para['x_'+flag]-para[flag+'_spacing'])+para[flag+'x_shift'],-x*para[flag+'_spacing']+0.5*(para['y_'+flag]-para[flag+'_spacing'])+para[flag+'y_shift'])
 	#Observe that (1,1) pixel in the image,array and sensor could be at different places depending on the spacing and dimensions of the object
	
def realToPixel(x,y,flag): # Inverse
	return(int(round(-100*(y-para[flag+'y_shift']-0.5*(para['y_'+flag]-para[flag+'_spacing']))/(100*para[flag+'_spacing']),0)),int(round(100*(x-para[flag+'x_shift']+0.5*(para['x_'+flag]-para[flag+'_spacing']))/(100*para[flag+'_spacing']),0)))

def pixelToRealAngle(x,y,flag): #For tilted angles
	x_real,y_real = pixelToReal(x,y,flag) # Get pixel values
	vector = [x_real,y_real,0] # Represent vector
	vector = np.matmul(inverse_transform,vector)
	return(vector[0],vector[1],vector[2])

def realToPixelAngle(x,y,z,flag):
	vector = [x,y,z] # Represent the real vector
	vector = np.matmul(transform,vector)
	return(realToPixel(x,y,flag)) # Return the pixel values

def propagateToSensor(image,image_distance): # Main function for getting the images
	distance_to_array = round(para['array_distance']-image_distance,5) # Distance to microlens array
	radius = calculateRadius(image_distance,distance_to_array) # Radius of the cone
	normalize = 1/float(radius*radius) # Intensity to be conserved
	print(radius)
	for x in range(image.shape[0]):
		for y in range(image.shape[1]):
			if image[x][y] == 0:
				continue
			x_real,y_real = pixelToReal(x,y,'image') # Get the coordinates of this pixel
			center_x,center_y = findIntersection(x_real,y_real,distance_to_array,0,0,para['array_distance'],-para['f_microlens']) # Sensor coordinates
			center_x,center_y = realToPixel(center_x,center_y,'sensor') # Get sensor pixel
			for x1 in range(center_x-radius,center_x+radius): # Within that circle
				extent = int(np.sqrt(radius*radius-(x1-center_x)*(x1-center_x)))
				for y1 in range(center_y-extent,center_y+extent):
					x1_real,y1_real = pixelToReal(x1,y1,'sensor')
					lensx_real,lensy_real = findIntersection(x1_real,y1_real,-para['f_microlens'],x_real,y_real,distance_to_array,0) #Find the lens coordinates
					sensor_x,sensor_y = sensorMapping(lensx_real,lensy_real,x1_real,y1_real,x_real,y_real,distance_to_array) # Find the sensor mapping
					if np.abs(sensor_x) >= para['x_sensor']/2: # Outside the sensor
						continue
					if np.abs(sensor_y) >= para['y_sensor']/2:
						continue
					sensor_x,sensor_y = realToPixel(sensor_x,sensor_y,'sensor')
					sensor[sensor_x][sensor_y] = sensor[sensor_x][sensor_y] + image[x][y]*normalize

def propagateToSensorAngle(image,image_distance): # Same as above expect everything involves tilted planes
	distance_to_array = round(para['array_distance']-image_distance,5)
	radius = calculateRadius(image_distance,distance_to_array)
	normalize = 1/float(radius*radius)
	for x in range(image.shape[0]):
		for y in range(image.shape[1]):
			if image[x][y] == 0:
				continue
			x_real,y_real = pixelToReal(x,y,'image')
			center_x,center_y = findIntersection(x_real,y_real,distance_to_array,0,0,para['array_distance'],-para['f_microlens']) # Sensor coordinates
			center_x,center_y = realToPixel(center_x,center_y,'sensor') # Get sensor pixel
			for x1 in range(center_x-radius,center_x+radius):
				extent = int(np.sqrt(radius*radius-(x1-center_x)*(x1-center_x)))
				for y1 in range(center_y-extent,center_y+extent):
					x1_real,y1_real = pixelToReal(x1,y1,'sensor')
					lensx_real,lensy_real,lensz_real = findIntersectionAngle(x1_real,y1_real,-para['f_microlens'],x_real,y_real,distance_to_array)
					sensor_x,sensor_y = sensorMappingAngle(lensx_real,lensy_real,lensz_real,x1_real,y1_real,x_real,y_real,distance_to_array)
					if np.abs(sensor_x) >= para['x_sensor']/2:
						continue
					if np.abs(sensor_y) >= para['y_sensor']/2:
						continue
					sensor_x,sensor_y = realToPixel(sensor_x,sensor_y,'sensor')
					sensor[sensor_x][sensor_y] = sensor[sensor_x][sensor_y] + image[x][y]*normalize

# Driving loop
for index in range(3,4): # 2,3 arbitrary. Can range between the number of slices
	section = cv2.imread("Stack/image"+str(index)+".png",0) 
	distance = para['specimen_distance'] +(index-3)*para['z_spacing'] # 2 here determines that slice in image2 is in focus
	image_distance,mag = lens.distanceMagnification(distance)
	print(image_distance,mag)
	if para['x_tilt'] == 1 and para['y_tilt'] == 1:
		propagateToSensor(section,image_distance)
	else:
		propagateToSensorAngle(section,image_distance)

section = cv2.normalize(sensor,section, 0, 255, cv2.NORM_MINMAX) # Normalise the entire image
cv2.imwrite("sensor_less.png",section)

