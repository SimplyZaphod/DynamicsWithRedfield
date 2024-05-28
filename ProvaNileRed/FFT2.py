import sys

import numpy as np
from scipy.fft import fftfreq, fft, fftshift
from scipy.fft import ifft
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.gridspec as GS
from scipy import signal
import os

def DFT(fx, t, e, n_t=1):
    dt = t[1]-t[0]
    y = np.zeros((len(fx)), dtype='complex')
    for i in range(0, len(fx), n_t):
        y[i] = integrate.simpson(fx*np.exp(-1j*t*e[i]), t, dx=dt)
    return y

def rescale_mean(y):
    y += -np.mean(np.real(y)) - 1j*np.mean(np.imag(y))
    return y
def rescale_minmax(y):
    min = np.min(np.real(y))
    y -= min
    min = np.min(np.imag(y))
    y -= 1j*min
    max = np.max(np.real(y))
    y -= max/2
    max= np.max(np.imag(y))
    y -= 1j*max/2
    return y
def square_ffty(y):
    return np.sqrt(y*np.conj(y))
def mirroring(x,y):
    x1 = -x[1:]
    y1 = np.conjugate(y[1:])
    x = np.hstack((x1,x))
    y = np.hstack((y1,y))
    y = y[np.argsort(x)]
    x = x[np.argsort(x)]
    return x,y
def reorder(x,y):
    y = y[np.argsort(x)]
    x = x[np.argsort(x)]
    return x,y

gs = GS.GridSpec(2,2)

for i in range(1, len(sys.argv)):
    name = sys.argv[i]

    a = np.loadtxt(name)
    x = a[:,0]
    y = a[:,1] + 1j*a[:,2]

    y = rescale_minmax(y)
    x,y = reorder(x, y)
    y *= np.exp(-0.05 * np.abs(x))

    # x, y = mirroring(x, y)
    # BH = signal.blackmanharris(len(y))
    # y*= BH
    
    plt.subplot(gs[0,:])
    plt.plot(x, np.real(y), label=f'Re({name})')
    plt.scatter(x, np.imag(y),s=0.1, label=f'Im({name})')
    # plt.xlim(0,50)
    e = np.linspace(-10., 10., len(y))
    ffty = DFT(y, x, e,n_t=20)
    e = e*0.658

    # ffty = fft(y)
    # dt = x[1]-x[0]
    # e = fftfreq(len(x), dt)
    # ffty = fftshift(ffty)
    # e = fftshift(e)

    if(False):
        ffty = square_ffty(ffty)
        plt.subplot(gs[1,:])
        plt.plot(e,np.real(ffty), label=f'Re({name})')
        plt.xlim(-1,1)
    else:
        plt.subplot(gs[1, 0])
        plt.plot(e, np.real(ffty), label=f'Re({name})')
        plt.xlim(-4, 4)
        # plt.xlim(0,1.5)
        plt.subplot(gs[1, 1])
        plt.plot(e, np.imag(ffty), label=f'Im({name})')
        plt.xlim(-4, 4)

plt.subplot(gs[0,:])
plt.legend()
# plt.subplot(gs[1, 0])
# plt.legend()
# plt.subplot(gs[1, 1])
# plt.legend()
plt.show()
