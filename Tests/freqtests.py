# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:11:02 2016

@author: Scott
"""
import numpy as np
import matplotlib.pyplot as plt

times = np.linspace(0,400*np.pi, 5000)
E1 = 10*np.sin(times)
E2 = 2*np.sin(times/2 - 0.5)
#E3 = 5*np.sin(1.5*times)
#E3 = 0.01
E3 = 10

E = E1 + E2 + E3

plt.figure(1)
plt.clf()
plt.plot(times, E)
plt.title("Electric field vs. time")

# Calculate energy
J1 = np.mean(np.abs(E1)**2)
J2 = np.mean(np.abs(E2)**2)
J3 = np.mean(np.abs(E3)**2)
print "J1,J2,J3:", J1, J2, J3, "Jsum:", J1 + J2 + J3
J = np.mean(np.abs(E)**2)
print "Jtot:", J

# Break down by frequency
print J1


dt = np.mean(np.diff(times))
freqHz = np.fft.rfftfreq(len(times), d = dt) # Makes equally spaced frequency bins
fr_fund = 1/(2*np.pi)
freq = freqHz/fr_fund

Eft = np.fft.rfft(E) # SERIAL FFT OPTION
pwr = np.absolute(Eft)**2

plt.figure(2)
plt.clf()
plt.plot(freq, pwr)
plt.title('Frequency vs Eft')
plt.xlim(-0.5, 2.5)
plt.xlabel('Frequency')
plt.ylabel('A.U.')

plt.figure(3)
plt.clf()
plt.plot(np.diff(freq))
plt.title("Difference between frequency bins (negligible)")

# Make power integrate to total energy, over all fundamental freqs
df = np.mean(np.diff(freq))

pwr = J * (pwr/(np.sum(pwr)*df)) # Normalized value = Total energy * Value / (Integral over all frequencies)

print "Integral is now:", np.sum(pwr)*df


plt.figure(3)
plt.clf()
plt.plot(freq, pwr)
plt.title('Frequency vs Eft')
plt.xlim(-0.5, 2.5)
plt.xlabel('Angular frequency / $\omega_{laser}$')
plt.ylabel('Joules / $\Delta\omega$ / cm')


cdt = freq < 0.2
Jlow = np.sum(pwr[cdt])*df

cdt = np.logical_and(freq > 0.9, freq < 1.1)
Jfund = np.sum(pwr[cdt])*df

cdt = np.logical_and(freq > 0.3, freq < 0.7)
Jhalf = np.sum(pwr[cdt])*df

cdt = np.logical_and(freq > 1.8, freq < 2.3)
Jtwo = np.sum(pwr[cdt])*df

cdt = np.logical_and(freq > 1.3, freq < 1.7)
J3havs = np.sum(pwr[cdt])*df

print "Jlow, Jfund, Jhalf, Jtwo, J3havs:", Jlow, Jfund, Jhalf, Jtwo, J3havs

#plt.figure(4)
#plt.clf()
#plt.plot(np.diff(freq))
#plt.title("Difference between frequency bins (negligible)")

#
#freqHz = np.fft.rfftfreq(len(times), d = dt/1e15) # Frequencies in Hz, for upcoming real FFT
#freq = freqHz/fr_fund # Frequencies in units of the fundamental
#
#Eft = np.fft.rfft(Ez, axis = 0) # SERIAL FFT OPTION
#pwr = np.absolute(Eft)**2
#
#Intensity = np.sum(Ez**2,0)
#pwr = np.absolute(Eft)**2
#FTmap_0_5 = np.sum(pwr[np.logical_and(freq > 0.3, freq < 0.7)],0)
#FTmap_1 = np.sum(pwr[np.logical_and(freq > 0.9, freq < 1.1)],0)
#FTmap_1_5 = np.sum(pwr[np.logical_and(freq > 1.3, freq < 1.7)],0)
#FTmap_2 = np.sum(pwr[np.logical_and(freq > 1.8, freq < 2.3)],0)
#FTmap_0_85 = np.sum(pwr[np.logical_and(freq > 0.7, freq < 0.9)],0)
#pwr_sum = np.sum(pwr, (1,2))