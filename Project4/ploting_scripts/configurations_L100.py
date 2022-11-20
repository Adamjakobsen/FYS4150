# %%
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import os 

#plotting configurations with MC=10^2 as the intermediate plot 

#import configuration for T=1, MC = 0
config_1_random0 = read_csv("../data/100/s1_cfg_L100_A0_mc1000000_burn0_t1.000.csv")

#import configurations for T = 1, MC = 10^2
config_1_random2 = read_csv("../data/100/s3_cfg_L100_A0_mc100_burn0_t1.000.csv")

#import configurations for T = 1, MC = 10^6
config_1_random6 = read_csv("../data/100/s3_cfg_L100_A0_mc1000000_burn0_t1.000.csv")


#import configuration for T=2.4, MC = 0
config_24_random0 = read_csv("../data/100/s1_cfg_L100_A0_mc1000000_burn0_t2.400.csv")

#import configurations for T = 2.4, MC = 10^2
config_24_random2 = read_csv("../data/100/s3_cfg_L100_A0_mc100_burn0_t2.400.csv")

#import configurations for T = 2.4, MC = 10^6
config_24_random6 = read_csv("../data/100/s3_cfg_L100_A0_mc1000000_burn0_t2.400.csv")


#import configuration for T = 2.4, MC = 0
config_1_random0 = config_1_random0.to_numpy()
config_1_random2 = config_1_random2.to_numpy()
config_1_random6 = config_1_random6.to_numpy()
#import configuration for T = 2.4, MC = 10^3


#import configuration for T = 2.4, MC = 10^6
config_24_random0 = config_24_random0.to_numpy()
config_24_random2 = config_24_random2.to_numpy()
config_24_random6 = config_24_random6.to_numpy()


#plot configurations as 3 subplots 
#T = 1
fig, ax = plt.subplots(1, 3, figsize=(25,10))
configs_1 = [config_1_random0, config_1_random2, config_1_random6]

MC = [0, 100, 100000]
for i in range(3):
    ax[i].imshow(configs_1[i], cmap = "Purples")
    ax[i].axis("off")
    ax[i].text(8, 8, 'T = 1 $[J/k_{B}]$, MC = %d' %MC[i],fontsize=20, bbox={'facecolor': 'white', 'pad': 10})

fig.savefig('../figs/configurations_L100_T1.pdf',dpi=300, bbox_inches='tight')

#T = 2.4
fig, ax = plt.subplots(1, 3, figsize=(25,10))
configs_24 = [config_24_random0, config_24_random2, config_24_random6]

for i in range(3):
    ax[i].imshow(configs_24[i], cmap = "Purples")
    ax[i].axis("off")
    ax[i].text(8, 8, 'T = 2.4 $[J/k_{B}]$, MC = %d' %MC[i], fontsize=20 ,bbox={'facecolor': 'white', 'pad': 10})

fig.savefig('../figs/configurations_L100_T24.pdf', dpi=300, bbox_inches='tight')


