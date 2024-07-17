import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import cm
import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2024-07-10
# Updated: 2024-07-17

"""
 The python script to solve the Cahn-Hilliard equation using
 an implicit pseudospectral algorithm
"""

Nsteps = 10000
dt = 0.1

N = 512
c_hat = torch.empty((N,N), dtype=torch.complex64,device=device)
dfdc_hat = torch.empty_like(c_hat)

c = torch.empty((Nsteps,N,N), dtype=torch.float32,device=device)

L = 64*np.pi
dx = L/N

noise = 0.1
c0 = 0.7

rng = np.random.default_rng(12345) # the seed of random numbers generator

c[0] = c0 + torch.tensor(noise*rng.standard_normal(c[0].shape),dtype=torch.float32, device=device)

plt.imshow(c[0].cpu().numpy())
plt.colorbar(cmap='RdBu_r')
plt.title('$c_0=%.1f$'% c0)
plt.savefig('cahn-hilliard-input.png')
plt.show()

print('c0 = ',c[0].mean().cpu().numpy())

W = 2.0
M = 1.0 # mobility
kappa = 0.5 #gradient coeficient

kx = ky = torch.fft.fftfreq(N, d=dx)*2*np.pi
Kx,Ky = torch.meshgrid(kx,kx,indexing ='ij')
K = torch.stack((Kx,Ky)).to(device)
K2 = torch.sum(K*K,dim=0)

# The anti-aliasing factor  
kcut = kx.max()*2.0/3.0 # The Nyquist mode
dealias = (torch.abs(K[0]) < kcut )*(torch.abs(K[1]) < kcut )

"""
 The interfacial free energy density f(c) = Wc^2(1-c)^2
"""
def finterf(c_hat):
    return kappa*torch.fft.ifftn(K2*c_hat**2).real 

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

cint = c[0].sum()

c_hat[:] = torch.fft.fftn(c[0])
for i in tqdm(range(1,Nsteps)):
    dfdc_hat[:] = torch.fft.fftn(dfdc(c[i-1])) # the FT of the derivative
    dfdc_hat *= dealias # dealising
    c_hat[:] = (c_hat-dt*K2*M*dfdc_hat)/(1+dt*M*kappa*K2**2) # updating in time
    c[i] = torch.fft.ifftn(c_hat).real # inverse fourier transform
    error = torch.abs(c[i].sum()-cint)/cint
    
print('c = ',c[-1].mean().cpu().numpy())

print('relative_error = ',error.cpu().numpy())

plt.imshow(c[-1].cpu().numpy(),cmap='RdBu_r', vmin=0.0, vmax=1.0)
plt.title('$c_0=%.1f$'% c0)
plt.savefig('cahn-hilliard-c0-%.1f.png'% c0)
plt.show()

from matplotlib import animation
from matplotlib.animation import PillowWriter

# generate the GIF animation

fig, ax = plt.subplots(1,1,figsize=(4,4))
im = ax.imshow(c[0].cpu().numpy(),cmap='RdBu_r', vmin=0.0, vmax=1.0)
cb = fig.colorbar(im,ax=ax, label=r'$c(x,y)$', shrink=0.8)
tx = ax.text(400,50,f't={(25*0*dt):.0f}',
         bbox=dict(boxstyle="round",ec='white',fc='white'))
ax.set_title(r'$c_0=%.1f$'% c0)

def animate(i):
    im.set_data(c[25*i].cpu().numpy())
    im.set_clim(0.0, 1.0)
    tx.set_text(f't={(25*i*dt):.0f}')
    return fig,

ani = animation.FuncAnimation(fig, animate, frames= Nsteps//25,
                               interval = 50)
ani.save('ch-c0='+str(c0)+'.gif',writer='pillow',fps=24,dpi=100)