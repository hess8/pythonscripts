                    
        #create folders with tag.value as name
        for iSet, varSet in enumerate(varsIn): 
            print varSet
            for tag, values in varSet.iteritems():
                tagtext= tag.lower()[1:] #kills the @ symbol and uses only lowercase for dir name
                for value in values:
                    relpath = tagtext+'.'+ str(value)+'/'
                    os.system('mkdir %s' % relpath)
                    try:
                        os.chdir(path1)
                    except OSError:
                        os.system('mkdir %s' % path1)
                        os.chdir(path2)