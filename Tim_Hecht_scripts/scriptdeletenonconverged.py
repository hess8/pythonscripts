import os,subprocess,shutil


file = open('outputanalysis.csv','r')
os.chdir('1element')
file.readline()
lines = file.readlines()
for line in lines:
  if float(line.split(",")[-1]) > 0.1:
    print line.split(",")[0]
    print line.split(",")[-1]
    foldertodelete = line.split(",")[0]
    try:
      shutil.rmtree(foldertodelete)
    except OSError:
      print "Folder not Found."
    