    
    oszicar = open("OSZICAR",'r')
    outputFile.write(folder + "\t\t" + oszicar.readlines()[-1].split()[2] + "\n")
    oszicar.close()
    file.close()