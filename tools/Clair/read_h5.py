import deepdish as dd
import numpy as np
import sys

file_name = str(sys.argv[1])
print(file_name.split(".")[0])
d = dd.io.load(file_name)
s = 0
for key in d.keys():
    s += sum(d[key])
for key in d.keys():
    print(key, "{:.4f}, {:.2f}%".format(round(sum(d[key]), 4), round(sum(d[key])/s, 4) * 100))
