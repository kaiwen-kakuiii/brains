import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

path = os.getcwd()

df = pd.read_csv(path+'/data/para_names_2d.txt', delim_whitespace = True, skiprows = 4)

para_names = np.array(list(filter(lambda x : x != 'time_series', np.array(df['Par']))))
length = int(len(para_names))
par_min = np.array(df['Min'][0:length])
par_max = np.array(df['Max'][0:length])
post_sam = np.loadtxt(path+'/data/posterior_sample_2d.txt')

fig= plt.figure(figsize = (10,10))
ax = fig.add_subplot(111)

post_sam_e = np.exp(post_sam[:,8])
ax.hist(post_sam_e, density = False, bins = 10, color = 'black')
#ax.set_xlim(par_min[8]-3*np.std(post_sam[:,8]), par_max[8]+3*np.std(post_sam[:,8]))
#ax.text(0.03, 1.02, para_names[8], transform=plt.gca().transAxes, fontsize = 9)
rounded_av = round(np.average(post_sam_e), 4)
rounded_median = round(np.median(post_sam_e),4)
rounded_av = str(rounded_av)
rounded_median = str(rounded_median)
ax.text(0.65, 0.92, rounded_av +' (*10^6 solar masses)', transform=plt.gca().transAxes, fontsize = 10)
ax.text(0.67, 0.95, rounded_median +' (*10^6 solar masses)', transform=plt.gca().transAxes, fontsize = 10)
ax.text(0.53, 0.92, 'Mean (Red): ', transform=plt.gca().transAxes, fontsize = 10)
ax.text(0.53, 0.95, 'Median (Blue): ', transform=plt.gca().transAxes, fontsize = 10)
#plt.axvline(par_min[8], 0, 1, color = 'red'  )
#plt.axvline(par_max[8], 0, 1, color = 'cyan')
plt.axvline(np.average(post_sam_e), 0, 1, color = 'red'  )
plt.axvline(np.median(post_sam_e), 0, 1, color = 'cyan')
ax.set_xlim(min(post_sam_e)-1, max(post_sam_e)+1)
ax.set_title('Central Supermassive Black Hole Mass, NGC7469')
ax.set_xlabel('Millions of Solar Masses')
ax.set_ylabel('Count')
plt.savefig('mass_dist.pdf', format = 'pdf')
plt.show()
