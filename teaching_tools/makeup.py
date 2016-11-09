from numpy import zeros,amax,amin,mean,median,tanh,nonzero, where
from copy import deepcopy


def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

exlines = readfile('exam.csv')
mlines = readfile('makeup.csv')
wtanh = 20.0
ptsAdd = 15.0

maxlen = max([len(exlines),len(mlines)])

data = zeros(maxlen,dtype = [('name', 'S20'),('netid', 'S20'),('exscore', float),('mscore', float),('newscore', float)])

for i,line in enumerate(exlines[1:]):
    info = line.strip().split(',')
    data[i]['name'] = line.split('"')[1]
    data[i]['netid'] = info[3]
    data[i]['exscore'] = info[5] 
           
data2 = deepcopy(data[::-1]) #reversed
    
for i,line in enumerate(mlines[1:]):
    info = line.strip().split(',')
    netid = info[2]
    if netid in data['netid']:
        idata = data['netid'].tolist().index(netid)
        if info[3] != '':
            data[idata]['mscore'] = info[3]

data2 = data[data["mscore"]!=0.0]

# print where(data[:]['exscore'] == 0)[0]#data['exscore'].nonzero()
# print
# print data['exscore'].nonzero()[0]

minscore = amin(data2['exscore'])
maxscore = amax(data2['exscore'])
print 'Old max:', maxscore
print 'Old min:', minscore
print 'Old average', mean(data2['exscore'])
print 'Old median', median(data2['exscore'])

for i,line in enumerate(data2['exscore']):
    oldscore = data2[i]['exscore']
    temp = oldscore*100/maxscore
    data2[i]['newscore'] = temp + ptsAdd *(tanh((100-temp)/wtanh))

print '\n\nNew max:', amax(data2['newscore'])
print 'New min:', amin(data2['newscore'])
print 'New average', mean(data2['newscore'])
print 'New median', median(data2['newscore'])

#Write to file
file = open('newscore.csv','w')
file.write('Last Name,First Name,Netid,score\n')
for i in range(len(data2)):
    file.write('{},{},{}\n'.format(data2[i]['name'],data2[i]['netid'],data2[i]['newscore']))
file.close()
print "Done"