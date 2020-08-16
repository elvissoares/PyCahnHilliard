import numpy as np
from matplotlib import pyplot as plt
# from scipy.fft import fft2, ifft2
import pyfftw
import multiprocessing
from pyfftw.interfaces.scipy_fftpack import fft2, ifft2
# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2020-08-16

pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
print('Number of cpu cores:',multiprocessing.cpu_count())

"""
 The python script to solve the Cahn-Hilliard equation using
 an implicit pseudospectral algorithm
"""

N = 512
c_hat = np.empty((N,N), dtype=np.complex64)
dfdc_hat = np.empty((N,N), dtype=np.complex64)
c = np.empty((N,N), dtype=np.float32)

dx = 0.1
L = N*dx

noise = 0.1
c0 = 0.3
c[:] = c0 + noise*np.random.standard_normal(c.shape)

# plt.imshow(c)
# plt.colorbar(cmap='viridis')
# # plt.title('$c_0=%.1f$'% c0)
# plt.savefig('cahn-hilliard-input.png')
# plt.show()

print('c0 = ',c.sum()*dx**2/L**2)

Nsteps = 2000
dt = 1.e-2

W = 2.0
M = 1.0 # mobility
kappa = 0.5 #gradient coeficient

kx = ky = np.fft.fftfreq(N, d=dx)*2*np.pi
K = np.array(np.meshgrid(kx , ky ,indexing ='ij'), dtype=np.float32)
K2 = np.sum(K*K,axis=0, dtype=np.float32)

# The anti-aliasing factor  
kmax_dealias = kx.max()*2.0/3.0 # The Nyquist mode
dealias = np.array((np.abs(K[0]) < kmax_dealias )*(np.abs(K[1]) < kmax_dealias ),dtype =bool)

"""
 The interfacial free energy density f(c) = Wc^2(1-c)^2
"""
def finterf(c_hat):
    return kappa*ifft2(K2*c_hat**2).real 

"""
 The bulk free energy density f(c) = Wc^2(1-c)^2
"""
def fbulk(c):
    return W*c**2*(1-c)*c**2

"""
 The derivative of bulk free energy density f(c) = Wc^2(1-c)^2
"""
def dfdc(c):
    return 2*W*(c*(1-c)**2-(1-c)*c**2)

for i in range(Nsteps):
    c_hat[:] = fft2(c)
    dfdc_hat[:] = fft2(dfdc(c)) # the FT of the derivative
    dfdc_hat *= dealias
    c_hat[:] = (c_hat-dt*K2*M*dfdc_hat)/(1+dt*M*kappa*K2**2) # updating in time
    c[:] = ifft2(c_hat).real # inverse fourier transform
    
print('c = ',c.sum()*dx**2/L**2)

plt.imshow(c)
plt.colorbar(cmap='viridis')
plt.title('$c_0=%.1f$'% c0)
plt.savefig('cahn-hilliard-c0-%.1f.png'% c0)
plt.show()