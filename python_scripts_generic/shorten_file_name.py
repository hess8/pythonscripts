import os,sys

### Main script ###)
# *****Frequently changed settings *****

dir = "I:\Bret's docs"
dirList = []
maxLen = 240
i = 0
for path, subdirs, files in os.walk(dir):
    for name in files:
        type = '.'+name.split('.')[-1]
        title = name.replace(type,'')
        newname =  os.path.join(path,title[:maxLen]+type)
        if len(newname)>maxLen:
            print newname
            nt = len(title)
            excess = len(newname)-maxLen
            print nt-excess
            if os.path.exists(newname):
                i+=1
                newname =  os.path.join(path,title[:nt-excess-len(str(i))]+str(i)+type)
            os.rename(os.path.join(path, name),newname)
            print
print 'done'
