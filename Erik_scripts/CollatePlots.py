'''
Created on Jul 30, 2014

@author: eswens13
'''

import sys, os, subprocess
import numpy as num

class CollatePlots:
    '''
    classdocs
    '''


    def __init__(self, structFile):
        ''' Constructor '''
        self.structFile = structFile
        self.structList = []
        
        self.getStructuresFromFile()
        
        self.plotFile = '/allStructPlots/'
        
        self.equivClasses = []
        self.initEquivClassesFromFile('structure_classes')
        
    def getStructuresFromFile(self):
        infile = open(self.structFile, 'r')
        structLines = infile.readlines()
        infile.close()
        
        self.structList = [line.strip() for line in structLines]    
    
    def initializeEquivClasses(self):
        self.equivClasses.append([0])
        self.equivClasses.append([1,2,4,8,16,32,64,128])
        self.equivClasses.append([3,12,48,192])
        self.equivClasses.append([5])
        self.equivClasses.append([6,9,18,33,72,96,132,144])
        self.equivClasses.append([7,11,13,14,19,35,49,50,76,112,140,176,196,200,208,224])
        self.equivClasses.append([10])
        self.equivClasses.append([15])
        self.equivClasses.append([17])
        self.equivClasses.append([20])
        self.equivClasses.append([21,42,69,81,84,138,162,168])
        self.equivClasses.append([22,41,73,97,104,134,146,148])
        self.equivClasses.append([23,43,77,113,142,178,212,232])
        self.equivClasses.append([24,36,66,129])
        self.equivClasses.append([25,26,37,38,70,74,82,88,98,100,133,137,145,152,161,164])
        self.equivClasses.append([27,39,78,114,141,177,216,228])
        self.equivClasses.append([28,44,52,56,67,131,193,194])
        self.equivClasses.append([29,46,53,58,71,83,92,116,139,163,172,184,197,202,209,226])
        self.equivClasses.append([30,45,54,57,75,99,108,120,135,147,156,180,198,201,210,225])
        self.equivClasses.append([31,47,55,59,79,115,143,179,205,206,220,236,241,242,244,248])
        self.equivClasses.append([34])
        self.equivClasses.append([40])
        self.equivClasses.append([51])
        self.equivClasses.append([60])
        self.equivClasses.append([61,62,124,188,199,203,211,227])
        self.equivClasses.append([63,207,243,252])
        self.equivClasses.append([65])
        self.equivClasses.append([68])
        self.equivClasses.append([80])
        self.equivClasses.append([85])
        self.equivClasses.append([86,89,101,106,149,154,166,169])
        self.equivClasses.append([87,93,117,171,174,186,213,234])
        self.equivClasses.append([90])
        self.equivClasses.append([91,94,103,110,118,122,155,157,167,173,181,185,217,218,229,230])
        self.equivClasses.append([95])
        self.equivClasses.append([102])
        self.equivClasses.append([105])
        self.equivClasses.append([107,109,121,151,158,182,214,233])
        self.equivClasses.append([111,123,159,183,222,237,246,249])
        self.equivClasses.append([119])
        self.equivClasses.append([125])
        self.equivClasses.append([126,189,219,231])
        self.equivClasses.append([127,191,223,239,247,251,253,254])
        self.equivClasses.append([130])
        self.equivClasses.append([136])
        self.equivClasses.append([150])
        self.equivClasses.append([153])
        self.equivClasses.append([160])
        self.equivClasses.append([165])
        self.equivClasses.append([170])
        self.equivClasses.append([175])
        self.equivClasses.append([187])
        self.equivClasses.append([190])
        self.equivClasses.append([195])
        self.equivClasses.append([204])
        self.equivClasses.append([215])
        self.equivClasses.append([221])
        self.equivClasses.append([235])
        self.equivClasses.append([238])
        self.equivClasses.append([240])
        self.equivClasses.append([245])
        self.equivClasses.append([250])
        self.equivClasses.append([255])
    
    def initEquivClassesFromFile(self, infile):
        self.equivClasses = []
        
        classFile = open(infile, 'r')
        classLines = classFile.readlines()
        classFile.close()
        
        for line in classLines:
            self.equivClasses.append(line.strip().split())
    
    def getEquivClasses(self):
        return self.equivClasses
     
    def makeEquivPlotsDir(self):
        subprocess.call(['mkdir','equivPlots'])
        subprocess.call(['mkdir','equivPlots/html_files'])
        
        for equivClass in self.equivClasses:
            subprocess.call(['mkdir','equivPlots/' + str(equivClass[0])]) 
    
    def fillEquivPlots(self):
        plotDir = os.getcwd() + '/8CStructPlots/'
        currDir = os.getcwd()
        
        os.chdir('equivPlots')
        for equivClass in self.equivClasses:
            for struct in equivClass:
                subprocess.call(['cp',plotDir + str(struct) + '_plot.png',str(equivClass[0])])  
      
        os.chdir(currDir)
    
    def createEquivHTML(self):
        for equivClass in self.equivClasses:
            stringEquivClass = [str(struct) for struct in equivClass]
            self.collate_plots('equivPlots/html_files/' + str(equivClass[0]) + '_equivplot.html', stringEquivClass)
    
    def collate_plots(self, outfile, structureList):
        '''Creates an HTML page with the plots and labels'''   
        numPerRow = 4
        
        collatefile  = open(outfile,'w')
        
        collatefile.write(' <html>\n <HEAD>\n<TITLE> 8C Structures </TITLE>\n</HEAD>\n')
        collatefile.write(' <BODY>\n <p style="font-size:36px">\n <table border="1">\n <tr>\n') #start row
    
        iImage = 0       
        for struct in structureList:
            iImage += 1  
            #Image and element under it         
            collatefile.write('<td><p><img src="%s%s" ></p><p style="text-align:center;font-size:36px;font-weight:bold">%s</p></td>\n' % 
                              (self.plotFile,struct + '_plot.png', struct + '  (' + '{0:08b}'.format(int(struct)) + ')'))
            if iImage % numPerRow == 0:
                collatefile.write('</tr>\n<tr>\n') #end of row, begin new
        collatefile.write(' </tr></table> \n') #end of row and table
                
        collatefile.write(' </BODY> </html>') #end of file
        collatefile.close()
        
if __name__ == "__main__":
    collate = CollatePlots(sys.argv[1])
    collate.makeEquivPlotsDir()
    collate.fillEquivPlots()
    collate.createEquivHTML()
    
    
    

