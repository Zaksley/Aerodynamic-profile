import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt

data = np.random.rand(4, 6)

for i in range(4):
    for j in range(6):
        data[i][j] *= 255

print(data)

heat_map = sb.heatmap(data, cmap="YlOrRd")
plt.xlabel("Values on X axis")
plt.ylabel('Values on Y axis')

plt.show()











"""
from PIL import Image 
from colour import Color
white = Color("white")
colors = list(white.range_to(Color("red"), 100))
print( tuple(int(255 * i) for i in (colors[98].rgb) ) )

# ==== Creating our image

pressure_map = Image.new('RGB', (10, 10), (255, 255, 255) )
pixels = pressure_map.load()


# =====
for i in range(10):
    for j in range(10):
        pixels[i, j] = tuple(int(255*i) for i in colors[i*10+j].rgb) 
        #print(pixels[i, j])
"""
"""
pressure_map = Image.new('RGB', (10, 10), (255, 255, 255) )
pixels = pressure_map.load()

print(pixels)
pressure_map.save("Pressure_Map", "PNG")
"""


