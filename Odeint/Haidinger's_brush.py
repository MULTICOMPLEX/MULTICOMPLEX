import numpy as np
np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

E = 5 #3

x = np.linspace(-E, E, 1000)
y = np.linspace(-E, E, 1000)

xx, yy = np.meshgrid(x, y)

def DoRotation(xspan, yspan, RotRad=0):
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)], [-np.sin(RotRad), np.cos(RotRad)]])

    x, y = np.meshgrid(xspan, yspan)
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))
    
#zz = np.exp(-(xx**2 + yy**2)) * (yy**2 + yy**2)
zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2)
#zz = np.exp(-(xx**2 + yy**2)) * np.cos(xx**2 - yy**2) * np.sin(xx**2 - yy**2)

#plt.imshow(zz, cmap = 'prism', extent=[-E, E, -E, E])

plt.imshow(zz, cmap = 'viridis', extent=[-E, E, -E, E])

#plt.quiver(xx, yy, zz*xx, zz*yy)
#plt.imshow(zz, cmap = 'gray', extent=[-E, E, -E, E])
#plt.contour(zz, extent=[-E, E, -E, E])

#ax = plt.axes(projection='3d')
#ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

#ax.view_init(60, 35)


plt.title("Haidinger's Brush")
plt.show()

fig  = plt.figure();
frames = [] # for storing the generated images

norm = plt.Normalize(np.min(zz), np.max(zz))
xx, yy = np.meshgrid(x, y)

for i in range(90):
    xx, yy = DoRotation(x, y, math.radians(i*4))
    zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2) * np.sin(math.radians(i*4))
    frames.append([plt.imshow(zz, cmap = 'viridis', norm = norm, extent=[-E, E, -E, E], animated=True)])


ani = animation.ArtistAnimation(fig, frames, interval=20, blit=True,
                                repeat_delay=0)


ani.save('c:/ffmpeg/movies/movie1.gif')
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2)

def data(i):
    xx, yy = np.meshgrid(x, y)
    zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2) * np.sin(math.radians(i*4))
    xx, yy = DoRotation(x, y, math.radians(i*4))
    surf[0].remove()
    surf[0] = ax1.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap='viridis', norm = norm, edgecolor='none')


x = np.linspace(-E, E, 90)
y = np.linspace(-E, E, 90)

xx, yy = np.meshgrid(x, y)
zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2)

ax1 = plt.axes(projection='3d')

surf = [ax1.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap='viridis', norm = norm, edgecolor='none')]

ani = animation.FuncAnimation(fig, data, 90, interval=20, repeat_delay=0)

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.set_xlim(-E,E)
ax1.set_ylim(-E,E)
ax1.set_zlim(np.min(zz), np.max(zz))

plt.title("Photon")

ani.save('c:/ffmpeg/movies/movie2.gif')
plt.show()


