#analysis script
import os, string
#Specify Directory to use
mainDir = '/bluehome/bch/vasprun/graphene.structures/half_graphane/'
#Specify the subdir
subdir = 'initial_relax'
dir = mainDir + subdir
os.chdir(dir)

def nstrip(list):
    import string
    list2 = []
    for string1 in list:   
        string2 = string1.strip("\n")
        list2.append(string2)
    return list2
    


outfile = open('half_graphane_initial.csv','w')

file = open('elements','r')  
elements = nstrip(file.readlines())
file.close()
print elements

file = open('energies','r')
energies = nstrip(file.readlines())
file.close()

file = open('stretch','r')
stretch = nstrip(file.readlines())
file.close()

file = open('distances','r')
distances = nstrip(file.readlines())
file.close()

outfile.write('Element,Calculated Energy,Stretch Energy,Distance\n')
for i in range(len(elements)):
    linei = elements[i]+','+energies[i]+','+stretch[i]+','+distances[i]+'\n'
    outfile.write(linei)
outfile.close()
print "done"
