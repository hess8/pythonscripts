    def weightPoints(self,eps):
        '''Find the volume of the Voronoi cell around each point, and use it to weight the point.
        Search a sphere of radius a few packing radii for neighbors.  Use the half vectors to these points 
        and the vectors to the walls to define the bounding planes.
        Vectors are first taken from each mesh point as the origin, 
        then displaced to their real positions in the cell for possible display'''
        allMPfacets = []
        self.IBZ.weights = []
        for ip,point in enumerate(self.IBZ.mesh):
            print ip,
            pointCell = cell()
            neighs,neighLbls = self.getNeighbors(point,self.IBZ,eps)
#             print 'neighLbls',neighLbls
            boundVecs = zeros(len(neighs)+ len(self.IBZ.bounds[0]),dtype = [('uvec', '3float'),('mag', 'float')]) 
            for iw, u in enumerate(self.IBZ.bounds[0]):    
                ro = self.IBZ.bounds[1][iw]
                d = ro-dot(point,u)
                boundVecs[iw]['uvec'] = u #unit vector stays the same for the plane
                boundVecs[iw]['mag'] = d
#                 print 'wall',iw,u, vec, norm(vec)
            for j, jpoint in enumerate(neighs):
                vec = (jpoint - point)/2
                mag = norm(vec)
#                 print 'neighs',j,jpoint, vec, norm(vec)
                boundVecs[j+len(self.IBZ.bounds[0])]['uvec'] = vec/mag
                boundVecs[j+len(self.IBZ.bounds[0])]['mag'] = mag
            boundVecs.sort(order = 'mag') 
            pointCell = getVorCell(boundVecs,pointCell,'point',eps)
            self.IBZ.weights.append(pointCell.volume)
            self.IBZ.vorVols.append(pointCell.volume)
             
            #For completeness,could update pointCell.center and pointCell.fpoints.  For brevity, we don't do this. 
 
            allMPfacets.append(pointCell.facets)
         
        print
        self.facetsMeshVCMathFile(self.IBZ,allMPfacets)
        wtot = sum(self.IBZ.vorVols)
        stdev = std(self.IBZ.vorVols)
        meanV = mean(self.IBZ.vorVols)
        volCheck = 0.1
        volErr = wtot - self.IBZ.volume        
        volErrRel = volErr/self.IBZ.volume
         
        print 'Total volume of point Vor cells',wtot,'vs IBZ volume', self.IBZ.volume
        print 'Relative volume error', volErrRel,'Abs volume error', volErr, 'Std dev/mean',stdev/meanV
        if not areEqual(wtot, self.IBZ.volume, volCheck*self.IBZ.volume):
#             print 'Total volume of point Vor cells',wtot,'vs IBZ volume', self.IBZ.volume
            sys.exit('Stop: point Voronoi cells do not sum to the IBZ volume.')
        else:
            print 'Point Voronoi cells volumes sum OK to within factor of {} of IBZ volume OK'.format(volCheck)  
        pf = len(self.IBZ.mesh)*4/3.0*pi*(self.rpacking)**3/self.IBZ.volume
        print 'Packing fraction (can be >1 from points near boundary', pf
        meshDet = open('../meshDetails.csv','a')
        N = len(self.IBZ.mesh)
        meshDet.write('{},{},{:6.3f},{:6.3f},{:6.3f}\n'.format(self.nTarget,N,stdev/meanV,self.meshEnergy/float(N),pf))
        meshDet.flush()
        meshDet.close()
        self.redistrWghts()