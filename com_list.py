# creator Aidan Ferguson on 3/24/24 for use with BRAINS

import argparse
import pandas as pd

# AJF initialize parsers
parser = argparse.ArgumentParser(description='create continuum files for BRAINS after using PyCali')
parser.add_argument('orpath', metavar='orpath', help='path the object folder (probably "original" directory) ') 
parser.add_argument('season_num', metavar='season_num', help='total number of seasons, i.e. "5"')
arg = parser.parse_args()
for i in range(int(arg.season_num)):
    sn = 'season'+str(int(i)+1)
    df = pd.read_csv(arg.orpath+'flux.lst_season'+str(int(i)+1), delim_whitespace = True)
    out = open(arg.orpath+'combined.txt_'+sn, 'w')
    names = df['name']
    for name in names:
        out.write(str(name)+'\n')
    out.close()
