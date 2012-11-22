import os
from string import *
import re
f=os.popen('ls', 'r')
a=f.read()
#=os.popen('ls', 'r')
b=f.flush()
#print a
#print a[0]
#.read().strip()
print 'f', f
print 'a', a
print ' '
print a[0:25]
print 'b', b

g = open('../fireballrun/answer.bas', 'r')
gread =g.read()
gs=str(gread)

h = gs.split()
print h
print h[5]

g = open('../fireballrun/answer.bas', 'r')
in1 = g.readlines()
print 'in1', in1
print 'line1',in1[0]
print 'line 3', in1[2]
print 'line 3', in1[2][1:10]
s3 = str(in1[2])
print 'line 3 word 4', s3.split()[3]

q=os.popen('awk -f awk_Tfinal ../fireballrun/fb', 'r')
temp=float(q.read())
print 'T-final', temp
print 'T-final*2', temp*2

g2 = open('../fireballrun/answer.bas', 'r').readlines()
print 'g2', g2

g3 = str(os.system('more ../fireballrun/answer.bas'))
print 'g3', g3
n=2
#str1 = './runCreateFB.bash ' + str(n)
#p=os.popen(str1, 'r')


