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
#j = 24
#print(f'this is min of {para_names[j]}: {min(post_sam[j])}')

cols = 5
rows = length // cols
if length % cols !=0:
	rows +=1
position = range(1,length+1)

fig = plt.figure(figsize = (12,12))

for i in range(length):
	ax = fig.add_subplot(rows, cols, position[i])
	fig.tight_layout(pad = 2)
	ax.hist(post_sam[:,i], density = False, bins = 100, color = 'black')
	ax.set_xlim(par_min[i]-3*np.std(post_sam[:,i]), par_max[i]+3*np.std(post_sam[:,i]))
	ax.text(0.05, 1.02, para_names[i], transform=plt.gca().transAxes, fontsize = 9)
	plt.axvline(par_min[i], 0, 1, color = 'red'  )
	plt.axvline(par_max[i], 0, 1, color = 'cyan')
plt.savefig('model_param_hist.pdf', format = 'pdf')
plt.show()
