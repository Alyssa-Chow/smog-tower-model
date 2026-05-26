# Project 2: Smog Tower Analysis
# PROGRAM DESCRIPTION:
# This program prompts the user for their weighted importance of cost and 
# efficiency and uses this to calculate the resulting voltage powering the 
# tower. Then, using physics equations it calculates the net force on the 
# particulate matter and uses two rounds of Euler's method calculations to 
# estimate first velocity then position values over 1 second. Then, the program 
# graphs Position vs. Time and Velocity vs. Time of this particulate matter. 
# Finally, the program calculates cost and change in concentration of 
# pollutants. The program outputs the initial conditions, changes in pollutant 
#concentrationm, and the yearly power cost of running the smog tower.
#
#------------------------------------------------------------------------------
# IMPORT LIBRARIES
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# CONSTANTS
#------------------------------------------------------------------------------

# Initial pollutant values as percents of air
NO2Init = 5 # initial percent of NO2 in air
PM10Init = 5 # initial percent of PM10 in air
PM25Init = 5 # initial percent of PM2.5 in air
ozoneInit = 5 # initial percent of ozone in air
o2 = 20 # percent of oxygen in air
totParticulate = NO2Init + PM25Init + PM10Init + ozoneInit # total particulate as a percent

# X-EFIELD CONSTANTS
q = -1 * totParticulate #test charge values (C)
df = 4 # length between electric charged plates (m)

# X-SURROUNDING PARTICLES CONSTANTS
hPlate = 5 # height of the plate (m)
e = 8.85 * 10**-12 # constant e naught (F/m)
infPlate = 100 # plate's range of charge (m)

#EULERS METHOD CONSTANTS
mPart = 1.3497 * 10**-25 # mass of 1 particle (kg)

#------------------------------------------------------------------------------
# UDFs - GET/CHECK INPUTS
#------------------------------------------------------------------------------

# accepts two numbers as user input for cost and efficiency weights and checks they are int.s
def getInput(prompt):
    while True:
        weight = input(prompt)
        try: 
            weight = int(weight)
            # check if input is in proper range
            if (weight < 0  or weight > 100): # if input is outside accepted range, output error message
                print("Please enter an integer in the range 0 - 100")
            else: # return input if it is in the correct range
                return int(weight)
        except ValueError: # if value entered is not an int, output error message
            print("Invalid number. Please enter an integer.")

# checks that both integers add up to 100
def checkInput():
    while True:
        efficiencyWeight = getInput("Please enter your efficiency weighting: ")
        costWeight = getInput("Please enter your cost weighting: ")
        
        # check if weights add up to 100
        if (efficiencyWeight + costWeight == 100):
            return efficiencyWeight, costWeight
        else:
            print("Your integers do not add up to 100. Please try again")           
#------------------------------------------------------------------------------
# UDFs - NET FORCE, VELOCITY, AND POSITION CALCULATIONS
#------------------------------------------------------------------------------

# Calculates the net force in the x-direction on the particle
def xNetForce(volt, q, df, hPlate, m, e):
   
    elecFieldForce = (volt * q) / df
    surrPartForce = 2 * (q / (volt * hPlate * e) + q / (4 * volt * hPlate)) - ((q * infPlate) / (4 * e * hPlate * volt) + (q * infPlate) / (hPlate * volt))
    xNetF = elecFieldForce - surrPartForce * 10**-12
   
    return xNetF

# Uses Euler's Method to calculate a velocity estimation of the particle
def EulerVel(vel0, t0Vel, tfVel, h, volt, q, df, hPlate, mPart, e):
  
    tVals = np.arange(t0Vel, tfVel + h, h) # set a list of time values from initial to final time in increments of "h"
    vVals = [vel0] # initialize matrix for velocity values
    
    # Calculates the velocity for a range of time values 
    for t in tVals[:-1]:
        xForce = xNetForce(volt, q, df, hPlate, mPart, e) # calculate net force
        accel = xForce / mPart # calculate acceleration
        v = vVals[-1] + h * accel # update velocity
        vVals.append(v) # add the most recent velocity to the list of velocity values
        
    return tVals, np.array(vVals)  

# Uses Euler's Method to calculate a position estimation of the particle
# Returns a list of time values, position values, and velocity values
def EulerPos(vel0, t0Vel, tfVel, h, volt, q, df, hPlate, mPart, e):
    
    tVals, vVals = EulerVel(vel0, t0Vel, tfVel, h, volt, q, df, hPlate, mPart, e) # gets list of times and velocity values using Euler's estimation of velocity
    
    xVals = [0] # initialize position
    x0 = 1 # previous x value is initialized to initial x position
    
    # Use velocity to calculate position
    for v in vVals[:-1]: # for every velocity value, calculate corresponding position values
        x = x0 + h * v # calculate a new x value based on previous x value and velocity
        xVals.append(x) # add the most recent position to the list of position values
        x0 = x # update x0 (previous x value) to current position for next iteration
        
    return tVals, np.array(xVals), vVals

#------------------------------------------------------------------------------
# FUNCTION CALLS
#------------------------------------------------------------------------------

# Get inputs from user for efficiency and cost weights
print("For the following prompts, please enter 2 positive integers that add up to 100 (0 being least important and 100 being most important)") # user prompt
efficiencyWeight, costWeight = checkInput() # assign and check user input values

# Calculates a voltage(V) based on the weighted efficiency for range of voltages [1.5 V,100 V]
volt = 0.985 * efficiencyWeight + 1.5 

# Initialize paramaters for Euler's Method Calculations
vel0 = -6.5 * 10**27 # initial velocity (m/s)
t0Vel = 0  # start time (s)
tfVel = 1  # end time (s)
h = 0.1 # step size for Euler's method

# Calculates The lists of time, position, and velocity values using Euler's Method
tVals, xVals, vVals = EulerPos(vel0, t0Vel, tfVel, h, volt, q, df, hPlate, mPart, e)
vVals *= 10**-27
xVals *= 10**-27

#------------------------------------------------------------------------------
# PLOTTING
#------------------------------------------------------------------------------

# Create subplots for Velocity vs. Time and Pos vs. Time
fig, axs = plt.subplots(2,1) # creates a figure for 2 sub plots with 2 rows and 1 column

# Plot Velocity vs Time on first subplot (top)
axs[0].plot(tVals, vVals, 'r') # plot time on the x-axis and velocity on the y-axis in red
axs[0].set_title("Euler's Approx. for Velocity(m/s) vs. Time(s)") # assign graph title
axs[0].set_xlabel("Time (s)")  # x-axis label
axs[0].set_ylabel("Velocity (m/s)") # y-axis label

# Plot Position vs Time on second subplot (bottom)
axs[1].plot(tVals, xVals, 'b') # plot time on the x-axis and position on the y-axis in blue
axs[1].set_title("Euler's Approx. for Position(m) vs. Time(s)") # assign graph title
axs[1].set_xlabel("Time (s)") # x-axis label
axs[1].set_ylabel("Position (m)") # y-axis label

# Automatically adjust subplot layout for optimal spacing
plt.tight_layout()

# Show plot
plt.show()

#------------------------------------------------------------------------------
# CONCENTRATION CALCULATIONS
#------------------------------------------------------------------------------

# Calculate the initial concentration of pollutant
initialPartConcentration = PM25Init + PM10Init + ozoneInit + NO2Init

# Calculate the final concentration of particulate
NO2Final = NO2Init # NO2 concentration stays the same from start to finish
# if the particle hits the plate at x = 7 m, then the final PM will be completely filtered and the final ozone will be only what the smog tower produces
if xVals[-1] <= -7: 
    PM25Final = 0 # all PM2.5 is filtered out
    PM10Final = 0 # all PM2.5 is filtered out
    ozoneFinal = 0.05*o2 # all initial ozone is filtered out, but ozone is created from running the smog tower
# if the particle falls short of the plate at x = 7m, then no particulate matter gets cleared out and ozone will increase from running the smog tower
else:
    PM25Final = PM25Init # no PM2.5 is filtered out
    PM10Final = PM10Init # no PM10 is filtered out
    ozoneFinal = 0.05*o2 + ozoneInit # ozone is created from running the smog tower in addition to the initial ozone in the air

finalPartConcentration = PM25Final + PM10Final + ozoneFinal + NO2Final # final concentration of the particles as a percent

# positive percent change indicates an increase in pollutant, and a negative percent change indicates a decrease in pollutant
percentChange = 100 * (finalPartConcentration - initialPartConcentration) / initialPartConcentration # percent change between final and initial pollutant concentration

#------------------------------------------------------------------------------
# YEARLY POWER COST CALCULATIONS
#------------------------------------------------------------------------------

# Cost of precipiator plates
elecFieldForce = (volt * q) / df # force due to the electric field (N), needed to calculated the cost of powering the plates
costPlates = abs(elecFieldForce) * 4 * 2.7778 * 10**-7 * 0.19 * 60**2 * 24 *365 # cost to power plates $ per year

#Cost of Fan 
costFan = 199.93 # cost to power fan in $ per year

# Total yearly power cost
totPwrCost = costPlates + costFan # total cost to power the smog tower in $ per year

#------------------------------------------------------------------------------
# OUTPUTS
#------------------------------------------------------------------------------
print("\n---------------------------------------")

# Initial conditions for inputted weights and pollutant concentration:
print(f"Efficiency Weight: {efficiencyWeight}")
print(f"Cost Weight: {costWeight}") 
print("---------------------------------------")

# Concentration Data:
print("Concentration Data:")
# Initial concentration pollutants
print(f"Initial Concentration of Pollutants: {initialPartConcentration}%") 
# Final concentration pollutants
print(f"Final Concentration of Pollutants: {finalPartConcentration}%")
# Changes in concentration
print(f"Total pollutant has changed by {percentChange:.2f}% from initial readings.")
print("---------------------------------------")

#Yearly Power Cost:
print(f"Required Voltage: {volt: .2f} V")
print(f"Yearly Power Cost: ${totPwrCost:.2f} per year")
print("---------------------------------------")
print(f"x last {xVals[-1]}")
