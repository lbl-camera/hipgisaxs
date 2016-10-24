
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt


def abort(msg):
  print msg
  sys.exit(1)
  
def load_matrix(fname):
  ext = fname.split('.')[-1]
  if ext == 'out': return np.loadtxt(fname)
  elif ext == 'edf' or ext == 'EDF':
    e = xio.edf.EDFFile(fname)
    data = e.data
    return np.array(data)
  else: abort('unknown file format '+ext)

def norm(data):
  sum = 0.
  for i in range(0, len(data)):
    for j in range(0, len(data[0])):
      sum += data[i][j] * data[i][j]
  return np.sqrt(sum)

def normalize(data):
  normval = norm(data)
  for i in range(0, len(data)):
    for j in range(0, len(data[0])):
      data[i][j] = data[i][j] / normval
  return data


ap = argparse.ArgumentParser()
ap.add_argument('-i')   ## input matrix 1
ap.add_argument('-j')   ## input matrix 2
ap.add_argument('-o')   ## output plot file
ap.add_argument('-m', default=None)   ## mask file
ap.add_argument('--hor', default=-1, type=int)   ## horizontal line cut
ap.add_argument('--ver', default=-1, type=int)   ## vertical line cut
args = ap.parse_args()

infile1, infile2, outfile, maskfile = None, None, None, None
if args.i: infile1 = args.i
else: abort('error: an input data file 1 is necessary')
if args.i: infile2 = args.j
else: abort('error: an input data file 2 is necessary')
if args.o: outfile = args.o
else: abort('error: an output plot filename is necessary')
if args.m: maskfile = args.m
if args.hor: horizontal = args.hor
if args.ver: vertical = args.ver

if horizontal < 0 and vertical < 0:
  abort('one of \'--hor\' or \'--ver\' is required')
if horizontal >= 0 and vertical >= 0:
  print('warning: both horizontal and vertical values given. selecting horizontal')
  vertical = -1

data1 = load_matrix(infile1)
data2 = load_matrix(infile2)
if maskfile is not None:
  mask = load_matrix(maskfile)
  data1 = data1 * mask
  data2 = data2 * mask

data1 = normalize(data1)
data2 = normalize(data2)

print len(data1[0]), len(data1)
print len(data2[0]), len(data2)
if len(data1[0]) != len(data2[0]) or len(data1) != len(data2):
  abort('error: the input patterns are not of the same size')

if horizontal >= 0:
  print 'horizontal cut', horizontal
  line1 = data1[horizontal,:]
  line2 = data2[horizontal,:]
elif vertical >= 0:
  print 'vertical cut', vertical
  line1 = data1[:,vertical]
  line2 = data2[:,vertical]

# print line1
# print line2

plt.plot(line1, '-', lw=3, color='#A11135', alpha=.7)
plt.plot(line2, '-', lw=3, color='#153AC0', alpha=.7)
plt.gcf().set_size_inches(10, 5)
plt.savefig(outfile, bbox_inches='tight')
