"""
Created on Sun Jul 20 15:10:43 2014

@author: Eryn Cangi
        Week 3
        PHYS 410
        
A program to demonstrate the inaccuracy of the Euler method over long data
sets for calculating the path of a curve of a differential equation. Figures
demonstrating the graphs given by the Euler method and reality are shown.
"""

import numpy as np
import matplotlib.pyplot as mp
from timeit import default_timer as timer

def runge_kutta(y, h):
      
    k1 = y
    k2 = y + (h/2)*k1 
    k3 = y + (h/2)*k2
    k4 = y + h*k3
    
    return (h/6)*(k1 + 2*k2 + 2*k3 + k4)

# Let x and v be simple lists of the positions and velocities over time, given 
# in units of meters and m/s. m = mass in kg, k = N/m. Initially the mass is 
# not moving and is resting at some location away from its equilibrium location.
# Lists are also defined for the exact solutions and the percent error of the 
# Euler method vs. the exact solutions.
x = [1]
v = [0]
exactx = [1]
exactv = [0]
exactE = [1]
percentErrorInX = [0]
percentErrorInV = [0]
E = [0.5]
m = k = 1

# The resonant frequency of any oscillating system.
w0 = np.sqrt(k/m)

# dt = an arbitrary change in time which determines how frequently we would be
# taking measurements. time_interval is the interval over which we will take 
# measurements.
dt = 0.0010 # tenths of a second
time_interval = np.arange(0, 2001, dt, dtype=float)
start = timer()

# Calculate the exact solutions, doesn't need to be done in a loop
exactx = np.cos(w0*time_interval)
exactv = -1*np.sin(w0*time_interval)
exactE = (0.5*m*exactv**2 + 0.5*k*exactx**2)

# For each defined time in the interval, the loop will increment or decrement
# the position or velocity using the Euler method. It then appends the 
# resulting value into a new list for plotting. Also, this loop fills in a
# list containing data about the energy over time.
for i in range(len(time_interval) -1):  
    x.append(x[i] + runge_kutta(v[i], dt))
    v.append(v[i] + (w0**2)*runge_kutta(-x[i], dt))
    E.append(0.5*m*(v[i]**2 + (w0**2)*(x[i]**2)))
    
# Just populates the percent error lists.    
for i in range(len(time_interval) -1):
    percentErrorInX.append((x[i] - exactx[i])/(exactx[i]))
    percentErrorInV.append((v[i] - exactv[i])/(exactv[i]))
    
dt = timer() - start

print("Elapsed time: %f s" % dt)

# This part just lets us pick a particular time to check the percent errors.
index = input("Enter a time index, representing  tenths of a second between 0"
            " and 1001 tenths of a second, at which you'd like to calculate the" 
            " error between the Euler method and reality: ")
            
print("Percent error for position: %0.2f%%" % (percentErrorInX[index]*100))
print("Percent error for velocity: %0.2f%%" % (percentErrorInV[index]*100))

# Shows a graph of the percent error lists, but it is ugly.
showErrorGraph = raw_input("Would you like to show the graph of error over "
                            "time? warning: it is hard to read. (Please enter "
                            "Y/N): """)

if showErrorGraph.lower() == 'y':
    showErrorGraph = True
else: 
    showErrorGraph = False

# Plots Figure 1, position over time.
mp.figure(1)
mp.plot(time_interval, x, 'g-')
mp.plot(time_interval, E, 'r-')
mp.legend(['Position', 'Energy'], loc='upper right')
mp.xlabel("Time")
mp.ylabel("Position")
mp.title("Position of harmonic oscillator using Runge-Kutta method")
mp.xlim(0, 287001)
mp.ylim(-10, 10)
mp.show()
  
# Plots Figure 2, velocity over time.  
mp.figure(2)
mp.plot(time_interval, v, 'b-')
mp.plot(time_interval, E, 'r-')
mp.legend(['Velocity', 'Energy'], loc='upper right')
mp.xlabel("Time")
mp.ylabel("Velocity")
mp.title("Velocity of harmonic oscillator using Runge-Kutta method")
mp.xlim(0, 2001)
mp.ylim(-10, 10)
mp.show()

# Plots Figure 3, exact solution for position, velocity and energy 
mp.figure(3)
mp.plot(time_interval, exactx, 'k-')
mp.plot(time_interval, exactv, 'b-')
mp.plot(time_interval, exactE, 'r-')
mp.legend(['Position', 'Velocity', 'Energy'], loc='upper right')
mp.xlabel("Time")
mp.ylabel("Position, Velocity and Energy")
mp.title("Exact solution for harmonic oscillator")
mp.xlim(0, 1001)
mp.ylim(-10, 10)
mp.show()

# Plots Figure 4, graph of percent error difference, based on user input
if showErrorGraph == True:
    mp.figure(4)
    mp.plot(time_interval, percentErrorInX, 'k-')
    mp.plot(time_interval, percentErrorInV, 'g-')
    mp.legend(['Error in Position', 'Error in Velocity'], loc='upper right')
    mp.xlabel("Time")
    mp.ylabel("Error")
    mp.title("Percent Error in Runge-Kutta method vs. reality")
    mp.xlim(0, 1001)
    mp.ylim(-10, 1000)
    mp.show()
else:
    pass
    
