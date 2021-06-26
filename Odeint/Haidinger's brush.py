import numpy as np
import matplotlib.pyplot as plt

E = 3;

x = np.linspace(-E, E, 1000)
y = np.linspace(-E, E, 1000)

xx, yy = np.meshgrid(x, y)
zz = np.exp(-(xx**2 + yy**2)) * (xx**2 - yy**2)

plt.imshow(zz, cmap = 'viridis', extent=[-E, E, -E, E])
#plt.imshow(zz, cmap = 'gray', extent=[-E, E, -E, E])

plt.title("Haidinger's Brush")
plt.show()