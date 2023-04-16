#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 19:17:09 2021

@author: Mizaan Katchera / Romain Joue CUPGE MP L2 CYU Cergy Paris Université
"""

#Il est conseillé d'ouvrir en grand écran la fenêtre d'affichage graphique.

#Enlève les erreurs.

import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

#Importation des bibliothèques utiles

import numpy as np
import math as m
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Début du programme 

#Définiton des Variables

g = 9.81 #constante de pesanteur
L = 1.0 #longueur
h=0.02 #pas
w=5 #pulsation

#CondtionsInitiales

angle0=45 #angle initial
theta0=np.radians(angle0) #angle initial
x0=2 #Vitesse
y0=np.array([theta0,x0]) #liste de la position verticale au cours du temps 
t0=0 #temps initial 
tps=np.arange(0,10,0.02) #liste de temps 
theta3=np.linspace(0,np.pi,50) #valeurs angles pi pour le diagramme de phase

#FonctionPendule

def pendule(y,t): #Fonction pendule
    theta = y[0] #liste position 
    theta_point = y[1] #liste vitesse
    return np.array([theta_point, -g/L*np.sin(theta)]) #equation du pendule

def pendule2(y,t): #Fonction pendule
    theta,theta_point = y #liste position 
    return np.array([theta_point, -np.sin(theta)]) #equation du pendule approchée pour diagramme de phase

#MéthodeEuler

def euler_bis (F, t0 , tf , y0 , n): #EulerAméliorée
    """ Fonction qui résout numériquement l’équation y ’=F(t,y)
    avec y(t0) = y0 par méthode d’Euler """
    pas = (tf - t0) / n #pas
    temps = np.linspace(t0 , tf , n + 1) #temps
    pos = [0 for i in range(len(temps))] #liste position
    pos[0] = y0  #liste position
    #pos = np.array(pos)
    for i in range(1 ,len(pos)):
        pos[i] = pos[i - 1] + pas * F(pos[i - 1], temps [i - 1] ) #position en fonction du temps et autres variables.
    return temps , np.array(pos)

"""x,y=euler(F,t0,4,0.1,20,tps,pos,vit)"""
x ,y = euler_bis(pendule, t0 , 10 , y0 , 10**6) #attribution variables pour tracé sur graphe.

#MéthodeEquaDiff

theta1=odeint(pendule,[theta0,x0],tps) #utilisation odeint et attribution variable pour tracé avec équa diff.

# Génération des subplots

fig, axs = plt.subplots(2,2, figsize=(10,8))

# Premier graphique (haut gauche)

axs[0,0].plot(x,y[:,0],lw=5.0,label='euler') #tracé Euler
axs[0,0].plot(tps,theta1[:,0],lw=1.5,label='equa diff') #tracé equa diff

axs[0,0].set_title("Trajectoire du pendule") #titre du graphe
axs[0,0].set_xlabel("Temps") #légende Ox
axs[0,0].set_ylabel("Position")  #légende Oy
axs[0,0].legend(loc='lower right') #déplacement légende
axs[0,0].grid() #grille

# Deuxième graphique (haut droit)

for i in theta3: #tracé pour différentes valeurs de theta
    tps1 = np.arange(0,25,0.02) #Liste de temps
    theta, theta_point = odeint(pendule2, (i, 0), tps1).T #équa diff
    axs[0,1].plot(theta, theta_point) #tracé différentes valeurs.
    
axs[0,1].set_title("Portrait de phase du pendule") #titre
axs[0,1].set_xlabel("Position") #légende Ox
axs[0,1].set_ylabel("Vitesse")  #légende Oy
axs[0,1].set_xlim(-np.pi, np.pi) #limitation axe Ox
axs[0,1].set_ylim(-2, 2) #limitation axe Oy
axs[0,1].grid() #grille

# Troisième graphique

legende=["\u03B8 = π/2","\u03B8 = π/3","\u03B8 = π/4","\u03B8 = π/6"] #légende pour différents angles
theta_arr=[np.pi/2,np.pi/3,np.pi/4,np.pi/6]

for k in range(0,4): #tracé pour différentes valeurs de theta
    theta = theta_arr[k]
    theta1=odeint(pendule,[theta,x0],tps) #équa diff
    theta,theta_point=theta1[:,0],theta1[:,1]#attribution variables
    axs[1,0].plot(tps,theta,label=legende[k])#tracé trajectoire pendule différents angles
    
axs[1,0].set_title("Trajectoire du pendule différents angles") #titre
axs[1,0].set_xlabel("Temps") #légende Ox
axs[1,0].set_ylabel("Position") #légende Oy
axs[1,0].legend(ncol=1, bbox_to_anchor= (1.02, 1.01), borderaxespad=0, frameon=False) #déplcement de la légende
axs[1,0].grid() #grille

# Supprimer le dernier subplot
gs = axs[1, 1].get_gridspec()
for ax in axs[1:, -1]:
    ax.remove()
    
fig.tight_layout()

plt.show()
