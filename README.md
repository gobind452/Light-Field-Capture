# Light Field Capture

Ray tracing code for simulating the capture of a light field in a light-field-microscopy setup

## Working

Light field microscopy involves using a microlens array and sensors to capture different light rays coming from an specimen, allowing refocussing and perspective changes. 

For simulating such a capture, we model a 3D object using a stack of slices in which each pixel is considered as a spherical source emitting rays in all directions proportional to its value. Each slice is then propagated through the objective lens, and from the image plane ray tracing occurs. For ray tracing, a cone of rays from each pixel is incident on the sensor. This cone has the objective lens as its base, and is extended to the sensor. Then the cone is discretized into rays by the number of sensor pixels it is hitting, and the microlens array effects are individually applied to these rays.

## Schematics

An image test.png of 200X200 pixels is supplied, along with various files used.

1. fourier.py - Fourier domain algorithms used (Not used in the current version)
2. stack.py - Given an image, it creates a stack of slices. This can be customised to create stacks of different nature.
3. objective.py - Simulates the effect of the objective lens (Not implemented)
4. drive.py - Contains the main code for simulating the capture

## Issues 

1. Only the convex hull of a 3D object can be modelled as a stack. For instance consider an cylinder whose thickness is the lowest at the center. The light coming from the back portion would be indistinguishable.
2. For accuracy, each point in the object should be considered a spherical source whose intensity decreases as 1/r*r where r is the distance from the soruce. As some rays travel more distance than others, these effect should be taken into account for better results.
3. Even if diffraction due to objective lens can be modelled easily using PSFs, one needs to take into account further diffraction as well which cant be done by PSFs as ray tracing is required.
4. Ray tracing is painfully slow, and optimisation is needed. (Maybe use C++? Multithreading?)
5. For image planes formed very close to the sensor, only a few rays are considered (as the radius of the cone incident from a pixel on the sensor is very less). This needs some changes in the algorithm.
6. No general way to create stacks

## Future Work

1. Support for RGB images
2. Explore algorithms in the Fourier space
3. Allow experimental errors to be encoded


