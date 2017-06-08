def getHessian(self,p):
    '''
    For a wall force in the form: f = (d/dw)^(-p), 
    D[en[d[x]], x, y] = -p ux uy (d/w)^(-1 - p)/w
   
    For an interparticle force in the form: f = (d/df)^(-p), the pair energy is 
    of the form 
    *Same type* component (x or y or z) on *different particles*:  Here x, y z represent any three components,
    reordered here:
      D[en[d[x1, x2, y1, y2]], x1, x2] =  
        (d/df)^-p/d^3 *(p(x1 - x2)^2 - y1^2 + 2 y1 y2 - y2^2 - z1^2 + 2 z1 z2 - z2^2)
        For example, for particles a, b, c,  d^2E/dx_a dx_b
        This applies to the one pair a-b
    
    Identical component on the *same* particle: D[en[d[x1, x2, y1, y2]], x1, x1] = -(above result)
        For example, for particles a, b, c,  d^2E/(dx_a)^2
        This applies to all pairs that include particle a
    
    *Different type* components (x vs y vs z) on same or different particles (x vs y vs z): 
      D[en[d[x1, x2, y1, y2]], x1, y1] = -(d/df)^-p/d^3 * (1 + p)(x1 - x2)(y1 - y2)
        For example, for particles a, b, c,  d^2E/dx_a dy_b or d^2E/dx_a dy_a 
        If comps are on the *same particle*, this applies to all pairs that include particle a 
        If comps are on *different particles*, this applies to the one pair a-b.  But use y2-1'''
    points = self.points
    N = len(points)
    self.hessian = zeros((3*N,3*N),dtype = real)
    df = self.df
    dw = self.dw 
    for ipoint in range(N):
        for icomp in range(3):
            for jpoint in range(ipoint+1,N): #2nd member of pair
                for jcomp in range(3):
                    print 'pair',ipoint,jpoint, 'comps',icomp,jcom
                    if ipoint == jpoint:
                        #wall forces: only same-particle contributes, but could be different comps
                        for iw, u in enumerate(self.bounds[0]):
                            ro = self.bounds[1][iw]
                            d = ro-dot(points[ipoint],u) #distance from plane to ri
                            if d<0:
                                d = -d/10 #Have crossed boundary. Increase force by shortening effective d. 
                                print 'ri,ro,u, dot(ri,u),d'
                                print ri,ro,u, dot(ri,u), d 
                                sys.exit('Error. Point {} in getHessian is not in the IBZ.'.format(iw))
                            self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += -(p*ux*uy*(d/dw)**(-1 - p))/dw
                        #differentiation coordinates on same particles
                        for mpoint in points:
                            if mpoint != ipoint: #forming all energy pairs that contain ipoint
                                comps = [0,1,2]
                                d = norm(points[ipoint]-points[mpoint])
                                x1 = points[ipoint][icomp]
                                x2 = points[mpoint][icomp]
                                comps.pop(icomp)
                                if icomp == jcomp:
                                    y1 = points[ipoint][comps[0]]
                                    y2 = points[mpoint][comps[0]]                 
                                    z1 = points[ipoint][comps[1]]
                                    z2 = points[mpoint][comps[1]]
                                    self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += -((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                                else:
                                    y1 = points[ipoint][jcomp]
                                    y2 = points[mpoint][jcomp]
                                    comps.pop(jcomp)
                                    z1 = points[ipoint][comps[0]]
                                    z2 = points[mpoint][comps[0]]
                                    self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += -((d/df)**-p)/d**3 * (1 + p)(x1 - x2)(y1 - y2)                 
                    else: #different particles
                        comps = [0,1,2]
                        d = norm(points[ipoint]-points[jpoint])
                        x1 = points[ipoint][icomp]
                        x2 = points[jpoint][icomp]
                        comps.pop(icomp)
                        if icomp == jcomp: #same comp type
                            y1 = points[ipoint][comps[0]]
                            y2 = points[jpoint][comps[0]]                 
                            z1 = points[ipoint][comps[1]]
                            z2 = points[jpoint][comps[1]]
                            self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += ((d/df)**-p)/d**3 *(p*(x1 - x2)**2 - y1**2 + 2*y1*y2 - y2**2 - z1**2 + 2*z1*z2 - z2**2)
                        else: #different comp type
                            y1 = points[ipoint][jcomp]
                            y2 = points[jpoint][jcomp]
                            comps.pop(jcomp)
                            z1 = points[ipoint][comps[0]]
                            z2 = points[jpoint][comps[0]]
                            self.hessian[ipoint*3+icomp,ipoint*3+jcomp] += -((d/df)**-p)/d**3 * (1 + p)(x1 - x2)(y2 - y1)            
                
                                  
                           
                    
                    
                    
    