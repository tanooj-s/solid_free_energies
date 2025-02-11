# python script to read lammps log file and write out to a csv file 
# the csv file can be loaded in a spreadsheet for further analysis

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="write lammps thermo output to a csv file")
parser.add_argument('-i', action="store", dest="input_file")
args = parser.parse_args()

# hack function to confirm tokens are numeric as floats
def is_num(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True

data = []
write_data = False
header_written = False
header_tokens = []


with open(args.input_file, 'r') as in_f:
    lines = in_f.readlines()

for line in lines:
    if 'ERROR' in line:
        break
    tokens = line.split()
    if len(tokens) > 0:
        if 'Total wall time:' in line:
            break
        if (tokens[0] == 'Step') and (header_written == False): # 'Step'
            write_data = True
            header_tokens = tokens
            data.append(tokens)
            header_written = True
        if header_written == True:
            s = [is_num(s) for s in tokens]
            if (sum(s) == len(s)) and (len(tokens) == len(header_tokens)):
                data.append(tokens)

out_file = args.input_file + ".csv"
with open(out_file,'w') as f:
    for line in data:
        f.write(','.join(line)+'\n')

