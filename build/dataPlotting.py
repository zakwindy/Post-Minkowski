import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = 'Files to plot.')
parser.add_argument('file1', type=argparse.FileType('r'), help='first file to be plotted', metavar='FILE')
parser.add_argument('file2', type=argparse.FileType('r'), help='second file to be plotted', metavar='FILE')

args = parser.parse_args()

data1 = pd.read_csv(args.file1)
data2 = pd.read_csv(args.file2)

args.file1.close()
args.file2.close()

fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')

data1.columns = ['x', 'y']
data2.columns = ['x', 'y']
ax.plot(data1.x, data1.y) 
ax.plot(data2.x, data2.y)

plt.show()

