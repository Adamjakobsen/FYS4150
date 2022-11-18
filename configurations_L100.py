# %%
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

# %%
#import configuration for T=1, MC = 0
config_1_random0 = read_csv("./100/s1_cfg_L100_A0_mc1000000_burn0_t1.000.csv")

#import configurations for T = 1, MC = 10^3
config_1_random3 = read_csv("./100/s2_cfg_L100_A0_mc1000000_burn0_t1.000.csv")

#import configurations for T = 1, MC = 10^6
config_1_random6 = read_csv("./100/s3_cfg_L100_A0_mc1000000_burn0_t1.000.csv")


#import configuration for T=2.4, MC = 0
config_24_random0 = read_csv("./100/s1_cfg_L100_A0_mc1000000_burn0_t2.400.csv")

#import configurations for T = 2.4, MC = 10^3
config_24_random3 = read_csv("./100/s2_cfg_L100_A0_mc1000000_burn0_t2.400.csv")

#import configurations for T = 2.4, MC = 10^6
config_24_random6 = read_csv("./100/s3_cfg_L100_A0_mc1000000_burn0_t2.400.csv")




# %%
#import configuration for T = 2.4, MC = 0
config_1_random0 = config_1_random0.to_numpy()
config_1_random3 = config_1_random3.to_numpy()
config_1_random6 = config_1_random6.to_numpy()
#import configuration for T = 2.4, MC = 10^3


#import configuration for T = 2.4, MC = 10^6
config_24_random0 = config_24_random0.to_numpy()
config_24_random3 = config_24_random3.to_numpy()
config_24_random6 = config_24_random6.to_numpy()


# %%
#plot configurations as 3 subplots 

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))
#plot each configuration using imshow 
ax1.imshow(config_1_random0, cmap='BuPu')
ax2.imshow(config_1_random3, cmap = 'BuPu')
ax3.imshow(config_1_random6, cmap = 'BuPu')
# ax.imshow(config_1_random0, interpolation='nearest', cmap="BuPu")
ax1.text(10, 5, 'L = 100, T = 1 $[J/k_{B}]$, MC = 0', bbox={'facecolor': 'white', 'pad': 10})
ax2.text(10, 5, 'L = 100, T = 1 $[J/k_{B}]$, MC = $10^3$', bbox={'facecolor': 'white', 'pad': 10})
ax3.text(10, 5, 'L = 100, T = 1 $[J/k_{B}]$, MC = $10^6$', bbox={'facecolor': 'white', 'pad': 10})
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
fig.savefig('configurations_T1.pdf', dpi=300, bbox_inches='tight')

#add colorbar

# %%
#plot configurations as 3 subplots 

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))
#plot each configuration using imshow 
ax1.imshow(config_24_random0, cmap = "BuPu")
ax2.imshow(config_24_random3, cmap = "BuPu")
ax3.imshow(config_24_random6, cmap = "BuPu")
# ax.imshow(config_1_random0, interpolation='nearest', cmap="BuPu")
ax1.text(10, 5, 'L = 100, T = 2.4 $[J/k_{B}]$, MC = 0', bbox={'facecolor': 'white', 'pad': 10})
ax2.text(10, 5, 'L = 100, T = 2.4 $[J/k_{B}]$, MC = $10^3$', bbox={'facecolor': 'white', 'pad': 10})
ax3.text(10, 5, 'L = 100, T = 2.4 $[J/k_{B}]$, MC = $10^6$', bbox={'facecolor': 'white', 'pad': 10})
ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
fig.savefig('configurations_T24.pdf', dpi=300, bbox_inches='tight')

# plt.show()


