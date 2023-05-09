misc.py has some miscellaneous functions and as well as constants including control parameters and physical parameters of the quadrotor

quad.py has the dynamic model of the quadrotor in the Quad class. The Quad class has several functions involving transformations of frames and also functions which return the state of the quadrotor.

hover.py sets a constant desired position 

triangle.py changes the desired position in each step according to a triangular trajectory

Running hover.py or triangle.py simulates the motion of the quadrotor

Execution times 
hover.py : < 1 second
triangle.py : Around 5 seconds