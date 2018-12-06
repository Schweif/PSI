#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#Calculates Reynolds Number and thermal conductivity for a simple tube model
#Version 1, Nov 2018
#David Just, Paul Scherrer Institut
#david.just@psi.ch

from __future__ import division
from termcolor import colored
import math

#Constants Water
roh_water = 997.05 # kg/m**3 density at 25 °C
Cp_water = 4179.0 #J/(kg*K) thermal capacity
lambda_water = 0.598 #W/(m*K) thermal conductivity
nue_water = 890.45e-6 #kg/(m*s) dynamic viscosity at 25 °C
mue_water = nue_water / roh_water #m**2/s kinetic viscosity


#Boundary conditions if 0 it will be calculated, set according to machine
P = 3900.0 #J/s = W heating power 
d= 0.004 #m diameter of water tube
r= d/2
l= 12*0.220 #m length of water tube
A= r**2*math.pi #m**2 cross section
T_i = 25.0 #°C water inlet temperature

#choose which parameter should be calculated (set to False if it should be calculated)
delta_T = 40.0 # K difference between inlet and outlet temperature, set to False to calculate delta T
v_flow_l = False #Volume flow  of the water in l/min, set to False to calculate volume current in l/min
model= 'Wagner' #Select a calculation model. Available models are:
#'Wagner'         = DEFAULT Formula 3.78, from Walter Wagner, Waermeuebertagung, Vogelfachbuch, 5. Ausgabe, 1998
#'Gnielinski'     = Gnielinski correlation, from: https://en.wikipedia.org/wiki/Nusselt_number, 9.11.2018
#'Dittus_Boelter' = Dittus-Boelter equation, from: https://en.wikipedia.org/wiki/Nusselt_number, 9.11.2018


def calc_water_flow_from_deltaT(P,Cp,roh,delta_T):
    #calculate water flow
    m_flow = P/(Cp*delta_T) #kg/s mass flow needed to cool heating power P
    v_flow = m_flow/roh #m**3/s volume flow needed to cool heating power P
    v_flow_l = v_flow*1000.0*60 #l/min volume flow in liter
    print "Flow needed is:                          " +str(round(v_flow_l,1)) +" l/min"
    return v_flow

def calc_deltaT_from_volumeFlow(P,Cp,roh,v_flow):
    m_flow = v_flow*roh
    delta_T = P/(m_flow*Cp)
    print "delta T needed is:                       " +str(round(delta_T,1)) +" K"
    return delta_T

def calc_pressure_loss(Re,d,L,omega,roh,k=0,n_bends=0,r_bend=0):
#calculate dracy frition factor lambda_tube (Rohrreibungszahl)
#Formulas form: HTBL-Kapfenberg Druckverlust in Rohrleitungen; https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=6&ved=2ahUKEwit1rHWqvzeAhWBl4sKHbMRBooQFjAFegQICRAC&url=http%3A%2F%2Fformeln.technis.ch%2FFormelsammlungen%2FFORMELNSAMMLUNG%2520STROMUNGSLEHRE1.pdf&usg=AOaw230lwH_PbpJvZ5VC3nbCQJ
    lambda_tube = 0.0
    if Re*(k/d) <= 65: #hydraulic flat surface
        if Re >= 2320 and Re <= 1e5:
            lambda_tube = 0.3164*Re**(-0.25) #darcy friction factor
            print colored('darcy friction factor calculated after Blasius','blue')
        elif Re > 1e5 and Re <= 5*e6:
            lambda_tube = 0.0032*+0.221*Re**(-0.237) #darcy friction factor
            print colored('darcy friction factor calculated after Nikuradse','blue')
        elif Re >= 1e6:
            #1/lambda_tube_Srt = 1/(2*math.log10(Re*lambda_tube_Srt) - 0.8 #darcy friction factor
            print colored('function to calculate darcy friction factor is not implemented. Try to lower Re or implement iterative solver','red')
            return none
    elif Re*(k/d) > 65 and Re*(k/d) <= 1300: #Uebergangsbereich
        print colored('function to calculate darcy friction factor is not implemented. Try to lower Re or implement iterative solver','red')
    elif Re*(k/d) > 1300: #hydraulic rough surface
        lambda_tube = 0.0055+0.15(k/d)**(1/3)#darcy friction factor
        print colored('darcy friction factor calculated after Moody','blue')
    delta_p = lambda_tube*L/d*roh/2*omega**2 # Pa pressure loss in tube
    return delta_p


if v_flow_l:
    print colored('volume flow ('+str(v_flow_l)+' l/min) is given, delta T will be calculated','blue')
    v_flow = v_flow_l/(1000.0*60) #volume flow of the water in m**3/s
    delta_T = calc_deltaT_from_volumeFlow(P,Cp_water,roh_water,v_flow)
else:
    print colored('delta T ('+str(delta_T)+' °C) is given, volume flow will be calculated','blue')
    v_flow = calc_water_flow_from_deltaT(P,Cp_water,roh_water,delta_T)



T_o =  T_i + delta_T #°C water outlet Temperature
T_1 = 0.0 #°C wall temperature on water side
T_2 = 0.0 #°C wall temperature warm side
T_av = (T_i+T_o)/2 #°C average water Temp




omega= v_flow/A #m/s flow speed of water
if omega <= 1.5:
    print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nIdeal velocity in terms of vibrations','green')
if omega <= 3.0 and omega >= 1.5:
    print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is ok','green')
if omega <= 4.0 and omega >= 3.0:
    print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is too high','grey','on_yellow')
if omega >= 4.0:
    print colored('Velocity omega is:                       '+str(round(omega,1))+' m/s\nVelocity is dangerously high, cavities might build up','white','on_red')

#calculate Reynolds Number
Re= omega*d/mue_water #dimension less
Pr = nue_water*Cp_water/lambda_water

print "Prandt Number Pr is:                     " +str(round(Pr,1))
if Pr>0.6 and Pr< 160:
    print colored("Prand Number ideal",'green')
if Pr<0.6:
    print colored("Prand Number too small",'mangenta','on_yellow')
if Pr>160: 
    print colored("Prand Number too large", 'red')

print "The Reynolds Number Re is:               "+str(round(Re,1))

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
    print colored("Flow condition is laminar, Nu will be calculated after Wagner 3.28",'blue')
    Re_status= "laminar"
    Nu =3.66



print "Nusselt Number Nu is: Nu =               " +str(round(Nu,1))

alpha= Nu*lambda_water/d
print "The heat transfer coefficient Alpha is:  " +str(round(alpha,1))+" W/(m**2*K)"
T_1 = P/(d*math.pi*l*alpha) +T_av
print "The average wall temperature is          " +str(round(T_1,1)) +" °C"
print "The outlet wall temperature is           " +str(round(T_1-T_av+T_o,1)) +" °C"
print "The inlet wall temperature is            " +str(round(T_1-T_av+T_i,1)) +" °C"

delta_p = calc_pressure_loss(Re,d,l,omega,roh_water)/100000 #pressure loss in tube in bar
print "The pressure loss in the tube is         " +str(round(delta_p,1)) +" bar"
if delta_p > 4:
    print colored("Warning: Pressure difference is too high, should be redesigned if possible",'red')

