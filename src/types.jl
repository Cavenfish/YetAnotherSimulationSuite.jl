"""
Abstract types used in JMD
"""

#MD image (a time snapshot)
abstract type MyImage end

#MD trajectory
abstract type MyTraj end

#Potential specific variables
abstract type PotVars end

#Atoms struct
abstract type MyAtoms end

#Cell struct
abstract type MyCell end

# Variables for thermostats
abstract type ThermoVars end

# Variables for barostats
abstract type BaroVars end

# Calculator
abstract type MyCalc end

# Thermostat
abstract type MyThermostat end

# Constraints
abstract type MyConstraint end