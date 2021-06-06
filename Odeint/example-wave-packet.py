from __future__ import division
import numpy as np
import sys, os
if 'save' in sys.argv:
        import matplotlib
        matplotlib.use('Agg')
import matplotlib.pyplot as plt

name = 'wave-packet'

usebarrier = 'barrier' in sys.argv
if usebarrier:
        name = 'tunneling'

usestep = 'step' in sys.argv
if usestep:
        name = 'step'

usespring = 'spring' in sys.argv
if usespring:
        name = 'spring'

usemovie = 'movie' in sys.argv

savemovie = 'save' in sys.argv and 'movie' in sys.argv

def newfig():
        if not savemovie and 'save' in sys.argv:
                plt.figure(figsize=(6,2.5))
        else:
                plt.figure(figsize=(12.8,7.2))

if usemovie and not savemovie:
        plt.ion() # this turns on interaction---needed for animation

m=1#9.11*10**-31
hbar=1#1.05*10**-34

L = 10
if usestep:
        L = 15

k=20
sig=0.5
if usestep:
        sig = 1.0
x0 = -L/4
p = k*hbar
K = p**2/(2*m)

tf = 1.5*L/(p/m)
if usespring:
        tf *= 3
        x0 = 0
dx=0.02
nx= int(L/dx)+1

x = np.linspace(-L/2,L/2,nx)            # x-array
dx = x[1]-x[0]
#potential
V = np.zeros_like(x)

if usebarrier:
        width = 0.1
        V[(0 < x) & (x < width)] = K*1.5
elif usespring:
        V = 2*K*x**2/(L/2)**2
elif usestep:
        V[x>0] = K*0.8

E = K + V[int((x0+L/2)/dx)]

kmax = np.pi/dx
Kmax = hbar**2*kmax**2/2/m
dt=np.pi*hbar/(Kmax + abs(max(V)))
nt = int(tf/dt)+1

t = np.linspace(0, tf,nt)
R=np.zeros((nt,nx))
I=np.zeros((nt,nx))

R[0]= np.exp(-(x-x0)**2/(2*sig**2))*np.cos(k*x)
I[0]=np.exp(-(x-x0 + p/m*dt/2)**2/(2*sig**2))*np.sin(k*x + K*dt/hbar/2)

#makes imaginary and real arrays
j=0
while j<nt-1:
	i=1
	while i<nx-1:
		I[j+1][i]=I[j][i]+((dt*hbar/(2*m))*(R[j][i+1]+R[j][i-1]-2*R[j][i])/dx**2)-((dt/hbar)*(V[i]*R[j][i]))
		i+=1
	i=1
	while i<nx-1:
		R[j+1][i]=R[j][i]-((dt*hbar/(2*m))*(I[j+1][i+1]+I[j+1][i-1]-2*I[j+1][i])/(dx**2))+((dt/hbar)*(V[i]*I[j+1][i]))
		i+=1
	j+=1

print('Done computing!')

Rt = np.matrix.transpose(R)
It = np.matrix.transpose(I)

if not usemovie:
        newfig()
        pmax = max(Rt.max(), -Rt.min())
        plt.pcolormesh(t,x,Rt, vmin=-pmax, vmax=pmax, cmap='bwr')
        plt.xlabel("$t$")
        plt.ylabel("$x$")
        plt.xlim([min(t), max(t)])
        plt.ylim([min(x), max(x)])
        plt.title("Real part")
        cbar=plt.colorbar()
        cbar.set_label(r'$\Re\psi(x,t)$')

        plt.savefig('%s-real-part.png' % name)

if 'imaginary' in sys.argv:
        newfig()
        pmax = max(It.max(), -It.min())
        plt.pcolormesh(t,x,It, vmin=-pmax, vmax=pmax, cmap='bwr')
        plt.xlim([min(t), max(t)])
        plt.ylim([min(x), max(x)])
        plt.xlabel("$t$")
        plt.ylabel("$x$")
        plt.title("Imaginary part")
        cbar=plt.colorbar()
        cbar.set_label(r'$\Im\psi(x,t)$')

if not usemovie:
        newfig()
        plt.pcolormesh(t,x,Rt**2 + It**2, cmap='hot')
        plt.xlim([min(t), max(t)])
        plt.ylim([min(x), max(x)])
        plt.xlabel("$t$")
        plt.ylabel("$x$")
        plt.title("Particle Probability Density Distribution")
        cbar=plt.colorbar()
        cbar.set_label('Probabilty')

if savemovie:
        os.system('rm -f /tmp/%s-*.png' % name)

fig = newfig()
#animation
anim_dt = 0.001
for this_t in np.arange(0,tf, anim_dt):  # run for one period
    j = int(this_t/dt)
    plt.clf() # erase the existing plot to start over
    plt.plot(x, R[j], label='Real') # plot the data
    Inow = 0.5*(I[j] + I[j+1])
    if usebarrier:
            plt.axvline(0)
            plt.axvline(width)
    if usestep:
            plt.axvline(0)
    plt.plot(x,Inow, label='Imaginary') #green
    plt.plot(x,np.sqrt(Inow**2 + R[j]**2), label='Magnitude')
    plt.xlabel('$x$') # label the axes
    plt.ylim([-2,2])
    plt.xlim([min(x), max(x)])
    plt.ylabel('$\psi(x)$')
    plt.title("$\psi(x)$ at $t=%g$" % this_t)
    plt.legend()
    if savemovie:
            plt.savefig('/tmp/%s-%06d.png' % (name, int(round(this_t/anim_dt))))
            if j % (nt//100) == 0:
                    print('%d%% done' % (100*j//nt))
    else:
            plt.pause(1e-9)
    if this_t == 0:
            plt.savefig('%s-figure.svg' % name)
    if int(this_t/anim_dt) == int(L/4*m/p):
            plt.savefig('%s-figure-later.svg' % name)
    if not usemovie:
            break

if min(V) != max(V) and not usemovie:
        newfig()
        plt.plot(x,V+ min(V), label='Potential')
        plt.plot(x,E +0*x, label='Energy')
        plt.xlim([min(x), max(x)])
        plt.xlabel('$x$')
        plt.ylabel("Energy")
        plt.legend(loc='best')
        plt.savefig('%s-potential.svg' % name)

if savemovie:
        avconv = "ffmpeg -y -r 60 -i /tmp/%s-%%06d.png -b 1000k %s-movie.mp4" % (name, name)
        print(avconv)
        os.system(avconv) # make the movie
else:
        plt.ioff()
        plt.show()
