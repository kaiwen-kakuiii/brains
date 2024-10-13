import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

path = os.getcwd()

df = pd.read_csv(path+'/data/para_names_2d.txt', delim_whitespace = True, skiprows = 4)

para_names = np.array(list(filter(lambda x : x != 'time_series', np.array(df['Par']))))

x_axis = str(input(f'What parameter to compare all with?\n'))
length = int(len(para_names))
for i in range(length):
    if x_axis == df['Par'][i]:
        indx = i
        print(f'index found of {x_axis}: {i}')



par_min = np.array(df['Min'][0:length])
par_max = np.array(df['Max'][0:length])
sample = np.loadtxt(path+'/data/sample_2d.txt', usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30))
print(f'this is sample 1 sample: {sample[:,0]}')
cols = 5
rows = length // cols
if length % cols !=0:
	rows +=1
position = range(1,length+1)

fig = plt.figure(figsize = (12,12))

for i in range(length):
	if i != indx:
		ax = fig.add_subplot(rows, cols, position[i])
		fig.tight_layout(pad = 2)
		ax.plot(sample[:,indx], sample[:,i], 'b.', markersize = 1)
		ax.text(0.05, 1.02, para_names[i], transform=plt.gca().transAxes, fontsize = 9)
                
plt.savefig(f'correltions_with_{x_axis}.pdf', format = 'pdf')
plt.show()
