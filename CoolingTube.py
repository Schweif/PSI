#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#Calculates Reynolds Number and thermal conductivity for a simple tube mode. The model assumes, that the power load is equally distributed
#Updates: 
#added possibility to caluclate system with parallel tubes (if they have the same diameter). 
#added outside wall temperature to output 
#added calculation for thermal radiation
#added pressure loss due to bends
#Version 1.4, Jan 2019
#David Just, Paul Scherrer Institut
#david.just@psi.ch

from __future__ import division
from termcolor import colored
import math

########################################################################
#CONFIGURATION SECTION

#Boundary conditions, set according to machine
P0= 315                 #J/s = W heating power
T_i = 20                #°C water inlet temperature

#Constants Water
roh_water = 1000        #kg/m**3 density at 25 °C
Cp_water = 4180         #J/(kg*K) thermal capacity
lambda_water = 0.6      #W/(m*K) thermal conductivity
nue_water = 890e-6      #kg/(m*s) dynamic viscosity at 25 °C

#Constants Solid
lambda_solid = 390      #W/(m*K) thermal conductivity of the wall material CuCr1Zr= 320, Glidcop =365, Cu =390 W/m*K)
epsylon_solid = 0.6     #w/o unit, emission number e.g. Cu polished = 0.04, Cu oxidized = 0.6, black colored 0.9

#Setup of the cooled device
thickness = 0.015       #m thickness of the material between fluid and power in
width =  0.021          #m irradiation width on which the power is applied
length = 0.007          #m irradiation length on which the power is applied

#Setup of the cooling channel
d= 0.004                #m diameter of water tube
l= 0.252                #m length of water tube in series e.g. if you have two parallel tube with a length l each, enter l
n= 1                    #number of tubes (with the same diameter) in parallel configuration flowing in the same direction

#choose which parameter should be calculated (set to False if it should be calculated)
delta_T = False         # K difference between inlet and outlet temperature, set to False to calculate delta T
v_flow_l_n = 2          #Volume flow  of the water in l/min through all parallel tubes, set to False to calculate volume current in l/min
model= 'Wagner'         #Select a calculation model for the Nusselt Number. Available models are:
#'Wagner'         = DEFAULT Formula 3.78, from Walter Wagner, Waermeuebertagung, Vogelfachbuch, 5. Ausgabe, 1998
#'Gnielinski'     = Gnielinski correlation, from: https://en.wikipedia.org/wiki/Nusselt_number, 9.11.2018
#'Dittus_Boelter' = Dittus-Boelter equation, from: https://en.wikipedia.org/wiki/Nusselt_number, 9.11.2018

#Experimental configurations (Pressure differences and counterflow)
k_tube = 5e-6           #m surface roughness of the cooling tube
n_bends = 0             #number of bend inside the cooling channel, set to 0 to ignore or if straight channel only
r_bend = 0.018          #m bending radius of bends inside the cooling channel
bend_angle = 180        #° bending angle of bends inside the cooling channel
counterflow= False      #boolean, specifiy if the cooling tubes run in counteflow mode or not; If in counterflow the average temperture will be used to qualify solid temperatures, default False
counterflow_factor= 0.5 #0..1, a factor defining the fraction between the efficinecy of the counterflow. if set to 1 the water inlet temperature wil be choosen as reference temperature for the maximum solid temperature calculations if set to 0 the maximum water temperature is used (i.e. without counterflow). Default is 0.5

#CONFIGURATION ENDS HERE
############################################################################

r= d/2
A= r**2*math.pi                     #m**2 cross section of the water tube
area = width*length                 #m**2 are on which the power is applied
P= P0/n                             #power divided into parallel tubes
v_flow_l= v_flow_l_n/n              #Volume flow  of the water in l/min through one single tube
mue_water = nue_water / roh_water   #m**2/s kinetic viscosity
BoltzmannConst = 5.6704e-8          #W/(m**2*K**4), Stefan Boltzmann Constant


def calc_water_flow_from_deltaT(P,Cp,roh,delta_T):
    m_flow = P/(Cp*delta_T)         #kg/s mass flow needed to cool heating power P
    v_flow = m_flow/roh             #m**3/s volume flow needed to cool heating power P
    v_flow_l = v_flow*1000.0*60     #l/min volume flow in liter
    v_flow_l_n = v_flow_l*n
    print "Flow needed is:                          " +str(round(v_flow_l_n,1)) +" l/min"
    return v_flow

def calc_deltaT_from_volumeFlow(P,Cp,roh,v_flow):
    m_flow = v_flow*roh
    delta_T = P/(m_flow*Cp)
    print "delta T needed is:                       " +str(round(delta_T,1)) +" K"
    return delta_T

def calc_pressure_loss_due_to_bends(Re,d,omega,roh,k_tube=0,n_bends=0,r_bend=1,bend_angle=0):
    #www.uni-magdeburg.de/isut/LSS/Lehre/Arbeitsheft/VIII.pdf p.7  
    Cang = 0.0                      #Constant adapting for bending angle
    Cre = 0.0                       #Constant adapting for Reynolds number
    Crou = 0.0                      #Constant abapting for surface roughtness
    if bend_angle <= 30.0 :
        Cang = 0.1
    elif bend_angle > 30.0 and bend_angle <= 45.0:
        Cang = 0.135
    elif bend_angle > 45.0 and bend_angle <= 60.0:
        Cang = 0.17
    elif bend_angle > 60.0 and bend_angle <= 90.0:
        Cang = 0.21
    elif bend_angle > 90.0 and bend_angle <= 200.0:
        Cang = 0.24
    elif bend_angle > 200.0:
        Cang = 0.26
    if Re > 3e3 and Re < 1e5:
        Cre = 20.2*Re**(-0.25)
    if Re > 1e5:
        Cre = 1
    if Re > 4e4:
        if 0 < k_tube/d and k_tube/d < 0.47*Re**(-0.75):
            Crou = 1.0
        if 0.47*Re**(-0.75) < k_tube/d and k_tube/d < 1e-3:
            Crou = 1+k_tube/d*1e3
        if k_tube/d > 1000:
            Crou = 2
    if Re < 4e4 and Re > 4e3:
            Crou = 0.8
    zeta = (Cang*Cre*Crou)/((r_bend/d)**(-0.5))
    if zeta == 0:
        print colored("The pressure loss due to bends could not be calculated","blue")
    delta_p_bends = zeta*roh*0.5*omega**2*n_bends
    return delta_p_bends

def calc_pressure_loss(Re,d,L,omega,roh,k_tube=0,n_bends=0,r_bend=1,bend_angle=0):
#calculate Dracy friction factor lambda_tube (Rohrreibungszahl)
#Formulas form: HTBL-Kapfenberg Druckverlust in Rohrleitungen; http://formeln.technis.ch/Formelsammlungen/FORMELNSAMMLUNG%20STROMUNGSLEHRE1.pdf
    lambda_tube = 0.0                                   #Darcy friction factor
    if Re*(k_tube/d) <= 65:                             #hydraulic flat surface
        if Re >= 2320 and Re <= 1e5:
            lambda_tube = 0.3164/(Re**(0.25))
            printOut = colored('Darcy friction factor for pressure loss calculation after Blasius','blue')
        elif Re > 1e5 and Re <= 5*e6:
            print '3'
            lambda_tube = 0.0032*+0.221*Re**(-0.237)
            printOut = colored('Darcy friction factor for pressure loss calculation after Nikuradse','blue')
        elif Re >= 1e6:
            print '4'
            #1/lambda_tube_Srt = 1/(2*math.log10(Re*lambda_tube_Srt) - 0.8
            printOut = colored('Function to calculate Darcy friction factor for pressure loss calculation is not implemented. Try to lower Re or implement iterative solver','red')
            return none
    elif Re*(k_tube/d) > 65 and Re*(k_tube/d) <= 1300:   #Uebergangsbereich
        printOut = colored('function to calculate Darcy friction factor for pressure loss calculation is not implemented. Try to lower Re or implement iterative solver','red')
    elif Re*(k_tube/d) > 1300:                           #hydraulic rough surface
        lambda_tube = 0.0055+0.15(k_tube/d)**(1/3)
        printOut = colored('Darcy friction factor for pressure loss calculation after Moody','blue')
    delta_p = lambda_tube*L/d*roh/2*omega**2             #Pa pressure loss in one tube
    delta_p_bends = calc_pressure_loss_due_to_bends(Re,d,omega,roh,k_tube,n_bends,r_bend,bend_angle)
    if delta_p_bends > 0.0:
        printOut = printOut +colored('\npressure loss due to bends: ' +str(delta_p_bends/100000/n),'blue')
    delta_p = delta_p + delta_p_bends
    return delta_p, printOut;

def calc_outside_wall_temp(d,A,lamda_solid,P,T_innerWall):
    T_outherWall= P*d/(lamda_solid*A)+T_innerWall
    return T_outherWall

def calc_thermal_radiation(A,T,epsylon):
    P_rad= epsylon*BoltzmannConst*(T**4)*A
    return P_rad

def evaluate_and_print_Prandt_Number():
    print "Prandt Number Pr is:                     " +str(round(Pr,1))
    if Pr>0.6 and Pr< 160:
        print colored("Prand Number ideal",'green')
    if Pr<0.6:
        print colored("Prand Number too small",'magenta','on_yellow')
    if Pr>160: 
        print colored("Prand Number too large", 'red')

def evaluate_and_print_velocity():
    if omega <= 1.5:
        print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nIdeal velocity in terms of vibrations','green')
    if omega <= 3.0 and omega >= 1.5:
        print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is ok','green')
    if omega <= 4.0 and omega >= 3.0:
        print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is too high','grey','on_yellow')
    if omega >= 4.0:
        print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is dangerously high, cavities might build up','white','on_red')


if v_flow_l_n:
    print colored('volume flow ('+str(v_flow_l_n)+' l/min) is given, delta T will be calculated','blue')
    v_flow = v_flow_l/(1000.0*60)   #volume flow of the water in m**3/s
    delta_T = calc_deltaT_from_volumeFlow(P,Cp_water,roh_water,v_flow)
else:
    print colored('delta T ('+str(delta_T)+' °C) is given, volume flow will be calculated','blue')
    v_flow = calc_water_flow_from_deltaT(P,Cp_water,roh_water,delta_T)



T_o =  T_i + delta_T                #°C water outlet Temperature
T_av = (T_i+T_o)/2                  #°C average water Temp
omega= v_flow/A                     #m/s flow speed of water


#calculate Reynolds  and Prandt Number
Re= omega*d/mue_water               #dimension less
Pr = nue_water*Cp_water/lambda_water

#calculate Nusselt Number
if Re>3000 and Re<5e6 and model == 'Gnielinski':
    print colored("Flow conditions is turbulent; Nu will be calculated after Gnieliski. See https://en.wikipedia.org/wiki/Nusselt_number",'blue') 
    Re_staus= "Gnieliski"
    f= (0.79*math.log(Re)-1.64)**(-2)
    Nu= (f/8)*(Re-1000)*Pr/(1+12.7*(f/8)**(0.5)*(Pr**(2/3)-1))

elif Re>10000 and model == 'Dittus_Boelter':
    print colored("Flow conditions is turbulent; Nu will be calculated after Dittus-Boelter. See https://en.wikipedia.org/wiki/Nusselt_number",'blue') 
    Re_staus= "Dittus_Boelter"
    Nu= 0.023*Re**(4/5)*Pr**(0.4)

elif Re>2300 and Re<1e6:
    print colored("Flow conditions is mixed turbulent; Nu will be calculated after Wagner 3.78",'blue')
    Re_staus= "Wagner"
    Nu= 0.0235*(Re**0.8-230)*Pr**0.48
    f6= 1+(d/l)**(2/3)

elif Re<2300:
    print colored("Flow condition is laminar, Nu is static after Wagner 3.28",'blue')
    Re_status= "laminar"
    Nu = 3.66


alpha= Nu*lambda_water/d                                                        #W/(m**2*K), heat transfer coefficient Alpha
T_chan_av = P/(d*math.pi*l*alpha) +T_av                                         #°C average wall temperature on the surface of the fluid channel
T_chan_max = P/(d*math.pi*l*alpha) +T_o                                         #°C maximum wall temperature on the surface of the fluid channel
T_chan_in = P/(d*math.pi*l*alpha) +T_i                                          #°C minimum wall temperature on the surface of the fluid channel
T_chan_counter = T_chan_in+(T_chan_max-T_chan_in)*counterflow_factor            #°C weigthed wall temperature on the surface of the fluid channel under counterflow condition
delta_p, printOutPressure= calc_pressure_loss(Re,d,l,omega,roh_water,k_tube,n_bends,r_bend,bend_angle)
delta_p = delta_p/100000/n                                                      #bar, pressure loss in tube in bar for all parallel tubes

if counterflow==False:
    T_outherWall_max = calc_outside_wall_temp(thickness,area,lambda_solid,P0,T_chan_max)
else:
    T_outherWall_max = calc_outside_wall_temp(thickness,area,lambda_solid,P0,T_chan_counter) 
P_rad = calc_thermal_radiation(area*7,T_outherWall_max+273.15,epsylon_solid)    #factor area*7 accounts for cubic volume from area plus some estimated factor

#generate output
evaluate_and_print_velocity()

if counterflow==False:
    if T_chan_max < 100:
        print colored("The max. channel wall temperature is:    " +str(round(T_chan_max,1)) +" °C",'green')
    else:
        print colored("The max. channel wall temperature is:    " +str(round(T_chan_max,1)) +" °C",'red')
else:
    if T_chan_max < 100:
        print colored("The max. channel wall temperature is:    " +str(round(T_chan_counter,1)) +" °C",'green')
    else:
        print colored("The max. channel wall temperature is:    " +str(round(T_chan_counter,1)) +" °C",'red')



print "The pressure loss in the tube is:        " +str(round(delta_p,2)) +" bar"
if delta_p > 4:
    print colored("Warning: Pressure difference is too high, should be redesigned if possible",'red')

print "The power loss due to radiation is:      " +str(round(P_rad,1)) +" W"

print  colored("\nAdditional vallues and information:", "blue")
evaluate_and_print_Prandt_Number()
print "The Reynolds Number Re is:               " +str(round(Re,1))
print "Nusselt Number Nu is:                    " +str(round(Nu,1))
print "The heat transfer coefficient Alpha is:  " +str(round(alpha,1))+" W/(m**2*K)"
print "The inlet channel wall temperature is:   " +str(round(T_chan_in,1)) +" °C"
print "The average channel wall temperature is: " +str(round(T_chan_av,1)) +" °C"
print "The max. outside wall temperature is:    " +str(round(T_outherWall_max,1)) +" °C"
print printOutPressure

print colored("\nThe values for delta_t and for the volume flow with accounting for thermal radiation are:","blue")
calc_deltaT_from_volumeFlow(P-P_rad,Cp_water,roh_water,v_flow)
calc_water_flow_from_deltaT(P-P_rad,Cp_water,roh_water,delta_T)


