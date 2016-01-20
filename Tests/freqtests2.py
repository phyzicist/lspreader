# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:11:02 2016

@author: Scott
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

times = np.linspace(0,400*np.pi, 5000)
E1 = 10*np.sin(times)
E2 = 2*np.sin(times/2 - 0.5)
#E3 = 5*np.sin(1.5*times)
#E3 = 0.01
E3 = 10

# Simulate the data file
eps = sc.epsilon_0 # Call the dielectric constant equal to vacuum permittivity, which is wrong within a plasma
mu = sc.mu_0 # Call the magnetic permeability the vacuum permeability, which is wrong within a plasma

Emag = (E1 + E2 + E3)/eps
Bmag = (B1 + B2 + B3)*mu
data = {}
data['times'] = times


# Begin real analysis
eps = sc.epsilon_0 # Call the dielectric constant equal to vacuum permittivity, which is wrong within a plasma
mu = sc.mu_0 # Call the magnetic permeability the vacuum permeability, which is wrong within a plasma

JE = (0.5*eps)*Emag**2 # Electric field energy density, Joule / m^3, for each frame
JB = (0.5/mu)*Bmag**2 # Magnetic field energy density, Joule / m^3, for each frames

J = JE + JB # Energy density, J / m^3
# Calculate total mean energy density (J/m^3) across frames (will be used to normalize)

Evec #(time x space x space x vector component)

def emFFT(Evec, data):
    """ Given the electric field vector (time x space x space x vector component), outputs the energy breakdowns"""
    nframes = len(data['times']) # Number of frames (times) in this data dict
    
    # Assume equal spacing in time and space, and get the deltas
    dt = np.mean(np.diff(data['times']))*1e-9 # dt in seconds
    dx = np.mean(np.diff(data['xgv']))*1e-2 # dx in m
    dz = np.mean(np.diff(data['zgv']))*1e-2 # dz in m

    # Calculate energy content of fields
    eps = sc.epsilon_0 # Call the dielectric constant equal to vacuum permittivity, which is wrong within a plasma
    mu = sc.mu_0 # Call the magnetic permeability the vacuum permeability, which is wrong within a plasma
    Jvecmean = (0.5*eps)*np.mean(Evec**2, 0) # Average the field energy density, for each component, along the time axis (axis 0). Now (space x space x vector component)
    Jvectot = np.sum(Jvecmean, axis=(0,1))*dx*dz # 1D array of total field energy per component, in Joules/m per component; Array dims: (component)

    
    # Calculate the frequency of the laser from its wavelength (800 nm)
    wl = 0.8e-6 # Wavelength of laser in m
    fr_fund = sc.c/wl # Laser frequency (the 'fundamental'), in Hz (Freq = Speed of light / wavelength)
    
    ## Break down by energy presence by frequency
    freqHz = np.fft.rfftfreq(nframes, d = dt) # Makes equally spaced frequency bins
    freq = freqHz/fr_fund # Re-define frequency from Hz to units of the laser frequency
    df = np.mean(np.diff(freq)) # Get the frequency step df

    Evecft = np.fft.rfft(Evec, axis=0) # Perform real frequency Fourier transform along time axis (axis 0)
    pwrvec = np.absolute(Evecft)**2 # Non-normalized value (freq x space x space x vector component)

    pwrvectot = np.sum(pwrvec, axis=(0,1,2))*df*dx*dz # Integral across space and time, for each vector component. 1D Arr dims: (component)
    for i in range(pwrvec.shape[3]): # Iterate over vector components
        pwrvec[:,:,:,i] = pwrvec[:,:,:,i] * (Jvectot[i] / pwrvectot[i]) # Normalized-to-J power. Each component integrates, over all frequencies and space, to the total energy (J/m) of that component
    # Now, each vector component of pwrvec integrates, over all frequencies and space, to the energy of that component
    # In other words, pwrvec is still (freq x space x space x vector), but an integral over all space, frequencies, and components gives the total field energy (J/m)
    
    return pwrvec
    
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