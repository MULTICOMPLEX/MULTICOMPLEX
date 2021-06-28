import numpy as np
import matplotlib.pyplot as plt

E = 5 #3

x = np.linspace(-E, E, 90)
y = np.linspace(-E, E, 90)

xx, yy = np.meshgrid(x, y)
np.set_printoptions(threshold=np.inf)

#zz = np.exp(-(xx**2 + yy**2)) * (yy**2 + yy**2)
zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2)
#zz = np.exp(-(xx**2 + yy**2)) * np.cos(xx**2 - yy**2) * np.sin(xx**2 - yy**2)

#plt.imshow(zz, cmap = 'prism', extent=[-E, E, -E, E])
#plt.imshow(zz, cmap = 'viridis', extent=[-E, E, -E, E])

#plt.quiver(xx, yy, zz*xx, zz*yy)
#plt.imshow(zz, cmap = 'gray', extent=[-E, E, -E, E])
#plt.contour(zz, extent=[-E, E, -E, E])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(xx, yy, zz, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

ax.view_init(60, 35)

plt.title("Haidinger's Brush")
plt.show()