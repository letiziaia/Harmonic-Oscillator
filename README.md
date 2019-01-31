# Harmonic-Oscillator
Studying 2nd order differential equation with Runge-Kutta and Verlet algorithms

_p1.c file approximates the solution of x'' = - alpha*x + A*x^2 - B*x^3
with parameters alfa=10.0, A=30, B=20, x0=-0.1, v=0.0
printing the answer in a text file which can be read by GnuPlot

_p2.c file studies the same equation in order to see how the trajectory changes
when x0 changes. Using the same algorithms, the critical value of x0 is 
approximated with 1% precision.

This project was part of the course in Computational Physics. 
