# Light Field Capture

Ray tracing code for simulating the capture of a light field in a light-field-microscopy setup

## Working

Light field microscopy involves using a microlens array and sensors to capture different light rays coming from an specimen, allowing refocussing and perspective changes. 

For simulating such a capture, we model a 3D object using a stack of slices. Each slice is then propagated through the objective lens, and from the image plane ray tracing occurs. For ray tracing, a cone of rays from each pixel is incident on the sensor. This cone has the objective lens as its base, and is extended to the sensor. Then the cone is discretized into rays by the number of sensor pixels it is hitting, and the microlens array effects are individually applied to these rays.

## Schematics

An image test.png of 200X200 pixels is supplied, along with various files used.

1. fourier.py - Fourier domain algorithms used (Not used in the current version)
2. stack.py - Given an image, it creates a stack of slices. This can be customised to create stacks of different nature.
3. objective.py - Simulates the effect of the objective lens (Not implemented)
4. drive.py - Contains the main code for simulating the capture

## Issues 

To be added

## Future Work

To be added
