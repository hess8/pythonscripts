import sys, os
import matplotlib.pyplot as plt
import numpy as np

class PlotGraphene:
    """ This plots either a POSCAR file or an output file from UNCLE.  When calling this class from
        the command line, call it as follows:
        
            python PlotGraphene.py  /path/to/folder/  filename  [ -p | -u ]  [-z]
        
        The "-p" or "-u" tags tell us whether we are reading a POSCAR from VASP or an output file 
        from UNCLE.  The default tag (if none is supplied) is "-p".  The /path/to/folder/ is the 
        path to the folder containing the file we want to read and the filename is the name of the 
        actual file we want to read.  The -z tag enables differentiation between points above and 
        below the plane. """
        
    def __init__(self, poscarDir, poscarFile, kindTag, zTag):
        if kindTag == "-u":
            self.uncle = True
        else:
            self.uncle = False
        
        self.poscarDir = poscarDir
        os.chdir(self.poscarDir)
        
        self.poscarFile = self.poscarDir + '/' + poscarFile
        
        if zTag == "-z" or zTag == "-Z":
            self.zFunc = True
        else:
            self.zFunc = False
        
        self.direct = False
        
        self.lattVec1 = []
        self.lattVec2 = []
        self.lattVec3 = []
        
        self.atomCounts = []
        self.xshift = 4.26257704
        self.yshift = 2.461
        
        self.origXlist = []
        self.origYlist = []
        self.origZlist = []
        
        self.xlist = []
        self.ylist = []
        self.zlist = []
        
        self.readPOSCAR()
        
        self.Ccirclelist = []
        self.Hcirclelist = []
        self.Mcirclelist = []
        
        self.figure = None
        self.initializeFigure()
    
        self.Hnum = 0
        self.Mnum = 0
        
    def readPOSCAR(self):
        
        poscarFile = open(self.poscarFile, 'r')
        poscarLines = poscarFile.readlines()
        
        # Extract the lattice vectors.
        if self.uncle:
            scrambledVec = []
            scrambledVec = poscarLines[2].strip().split()
            self.lattVec1 = [float(scrambledVec[1]), float(scrambledVec[2]), float(scrambledVec[0])]
            
            scrambledVec = poscarLines[3].strip().split()
            self.lattVec2 = [float(scrambledVec[1]), float(scrambledVec[2]), float(scrambledVec[0])]
        else:
            self.lattVec1 = poscarLines[2].strip().split()
            self.lattVec1 = [float(comp) for comp in self.lattVec1]
            
            self.lattVec2 = poscarLines[3].strip().split()
            self.lattVec2 = [float(comp) for comp in self.lattVec2]
            
            if sum(self.lattVec1[:2]) == 0.0: 
                self.lattVec1 = poscarLines[4].strip().split()
                self.lattVec1 = [float(comp) for comp in self.lattVec1] 
            if sum(self.lattVec2[:2]) == 0.0: 
                self.lattVec2 = poscarLines[4].strip().split()
                self.lattVec2 = [float(comp) for comp in self.lattVec2]      
                
            print'lv1',self.lattVec1  
            print'lv2',self.lattVec2       
        
        # Get the number of each type of atom.
        if self.uncle:
            nonCatomCounts = poscarLines[5].strip().split()
            nonCatomCounts = [int(count) for count in nonCatomCounts]
            
            self.atomCounts = [sum(nonCatomCounts), nonCatomCounts[0], nonCatomCounts[1]]
        else:
            self.atomCounts = poscarLines[5].strip().split()
            self.atomCounts = [int(count) for count in self.atomCounts]
        
        # Determine whether to compute atom positions in direct or cartesian coordinates.
        if poscarLines[6].strip().lower() == "d" or poscarLines[6].strip().lower() == "direct":
            self.direct = True
                
        positionLines = poscarLines[7:]
        
        if self.direct:
            if self.uncle:
                # If it is an UNCLE file, they don't explicitly state the C positions, so we run
                # this twice.  Once for the C positions and once for the H and M positions.
                self.getPositionsFromDirectCoordinates(positionLines)
                self.getPositionsFromDirectCoordinates(positionLines)
            else:
                self.getPositionsFromDirectCoordinates(positionLines)
        else:
            self.getPositionsFromCartesianCoordinates(positionLines)        
            
    def getPositionsFromDirectCoordinates(self, positionLines):
        zcoord = 1.0
        for line in positionLines:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            
            self.origXlist.append(float(comp1[0] + comp2[0]))
            self.xlist.append(float(comp1[0] + comp2[0]))
            self.origYlist.append(float(comp1[1] + comp2[1]))
            self.ylist.append(float(comp1[1] + comp2[1]))
            self.origZlist.append(zcoord)
            self.zlist.append(zcoord)
            
            if zcoord == 1.0:
                zcoord = -1.0
            else:
                zcoord = 1.0
        
    def getPositionsFromCartesianCoordinates(self, positionLines):
        for line in positionLines:
            positions = line.strip().split()
            self.origXlist.append(float(positions[0]))
            self.xlist.append(float(positions[0]))
            self.origYlist.append(float(positions[1]))
            self.ylist.append(float(positions[1]))
            self.origZlist.append(float(positions[2]))
            self.zlist.append(float(positions[2]))
    
    def initializeFigure(self):
        self.figure = plt.gcf()
        self.figure.gca().set_aspect('equal')
        plt.axis([0, 15, 0, 15])
    
    def getAtomCounts(self):
        return self.atomCounts
       
    def periodicByVecs(self, vec1num, vec2num):
        
        self.Ccirclelist = []
        self.Hcirclelist = []
        self.Mcirclelist = []
        
        for i in range(0,sum(self.atomCounts)):
            self.xlist[i] = self.origXlist[i] + (vec1num * self.lattVec1[0]) + (vec2num * self.lattVec2[0])
            self.ylist[i] = self.origYlist[i] + (vec1num * self.lattVec1[1]) + (vec2num * self.lattVec2[1])
            
        numOfC = self.atomCounts[0]
        numOfH = self.atomCounts[1]
        numOfM = self.atomCounts[2]
        
        Cxlist = self.xlist[0:numOfC]
        Cylist = self.ylist[0:numOfC]
        Czlist = self.zlist[0:numOfC]

        Hxlist = self.xlist[numOfC:numOfC + numOfH]
        Hylist = self.ylist[numOfC:numOfC + numOfH]
        Hzlist = self.zlist[numOfC:numOfC + numOfH]

        Mxlist = self.xlist[numOfC + numOfH:sum(self.atomCounts)]
        Mylist = self.ylist[numOfC + numOfH:sum(self.atomCounts)]
        Mzlist = self.zlist[numOfC + numOfH:sum(self.atomCounts)]
    
        for i in range(0,numOfC):
            self.Ccirclelist = self.Ccirclelist + [plt.Circle((Cxlist[i], Cylist[i]), .7, color='#2B65EC')]
    
        for i in range(0,numOfC):
            self.figure.gca().add_artist(self.Ccirclelist[i])
            
        for i in range(0,numOfH):
            if Hzlist[i] < 0 and self.zFunc:
                self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, facecolor='#FFD801', edgecolor='black', lw = 3)]
                self.Hnum += 1
            else:
                self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, color='#FFD801')]
                self.Hnum += 1
        
        for i in range(0,numOfH):
            self.figure.gca().add_artist(self.Hcirclelist[i])

        for i in range(0,numOfM):
            if Mzlist[i] < 0 and self.zFunc:
                self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, facecolor='r', edgecolor='black', lw = 3)]
                self.Mnum += 1
            else:
                self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, color='r')]
                self.Mnum += 1
                
        for i in range(0,numOfM):
            self.figure.gca().add_artist(self.Mcirclelist[i])
    
    def fillByVecs(self, num):
        if num == 0:
            self.periodicByVecs(0, 0)
        else:
            for i in xrange(0, num + 1):
                for j in xrange(0, num + 1):
                    self.periodicByVecs(i, j)
                    self.periodicByVecs(-i, j)
                    self.periodicByVecs(i, -j)
                    self.periodicByVecs(-i, -j)
    
    def addLines(self):
        if self.lattVec1[0] == 0:
            slope1 = 'vertical'
        else:
            slope1 = self.lattVec1[1] / self.lattVec1[0]
        
        if self.lattVec2[0] == 0:
            slope2 = 'vertical'
        else:
            slope2 = self.lattVec2[1] / self.lattVec2[0]
        
        xPoints = np.arange(-20, 20, .1)
        print'slope1, slope2', slope1, slope2
        for i in range(10):
            self.shiftLineByVec1(i, xPoints, slope2)
            self.shiftLineByVec1(-i, xPoints, slope2)
            self.shiftLineByVec2(i, xPoints, slope1)
            self.shiftLineByVec2(-i, xPoints, slope1)
    
    def shiftLineByVec1(self, ntimes, xPoints, slope):
        if slope == 0:
            plt.plot(xPoints, (slope * xPoints) + (-ntimes * self.lattVec1[1]), color='k')
        elif slope == 'vertical':
            plt.plot(((ntimes * self.lattVec1[0]), (ntimes * self.lattVec1[0])), (0, 15), color='k')
        else:
            plt.plot(xPoints, 
                     (slope * xPoints) + (slope * (-ntimes * self.lattVec1[0])) + (ntimes * self.lattVec1[1]), color='k')
    
    def shiftLineByVec2(self, ntimes, xPoints, slope):
        if slope == 0:
            plt.plot(xPoints, (slope * xPoints) + (-ntimes * self.lattVec2[1]), color='k')
        elif slope == 'vertical':
            plt.plot(((ntimes * self.lattVec1[0]), (ntimes * self.lattVec2[0])), (0, 15), color='k')
        else:
            plt.plot(xPoints, 
                     (slope * xPoints) + (slope * (-ntimes * self.lattVec2[0])) + (ntimes * self.lattVec2[1]), color='k')
        
    def saveFigure(self):
        plt.plot()
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.savefig(self.poscarFile + "_plot.png", bbox_inches='tight')
 
#==================================================================================================       
#============================================= MAIN ===============================================
#==================================================================================================

if __name__ == "__main__":
            
    error = False
    if len(sys.argv) > 5 or len(sys.argv) < 3:
        error = True
        print "USAGE: python PlotGraphene.py  /path/to/folder/  filename  [ -p | -u ]  [-z]"
    
    elif len(sys.argv) == 5:
        plotter = PlotGraphene(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    
    elif len(sys.argv) == 4:
        if sys.argv[3] == "-z" or sys.argv[3] == "-Z":
            plotter = PlotGraphene(sys.argv[1], sys.argv[2], "-p", "-z")
        elif sys.argv[3] == "-p" or sys.argv[3] == "-P":
            plotter = PlotGraphene(sys.argv[1], sys.argv[2], "-p", "!")
        elif sys.argv[3] == "-u" or sys.argv[3] == "-U":
            plotter = PlotGraphene(sys.argv[1], sys.argv[2], "-u", "!")
        else:
            error = True
            print "USAGE: python PlotGraphene.py  /path/to/folder/  filename  [ -p | -u ]  [-z]"
    
    elif len(sys.argv) == 3:
        plotter = PlotGraphene(sys.argv[1], sys.argv[2], "-p", "!")
    
    if not error:
        if plotter.getAtomCounts()[0] > 4:
            plotter.fillByVecs(5)
        else:
            plotter.fillByVecs(10)
            
        plotter.addLines()
        plotter.saveFigure()
            
    