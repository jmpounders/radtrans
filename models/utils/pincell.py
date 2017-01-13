def annulusPoints(nRad,nAzim,R,Ri,square=False) :
    Atot = np.pi*(R**2-Ri**2)
    Aann = Atot/nRad
    
    # select the radii to conserve annular volume
    rs = np.zeros(nRad)
    rprev = Ri
    for i in range(nRad) :
        rs[i] = np.sqrt(Aann/np.pi + rprev**2)
        rprev = rs[i]
    # select the radii to be equally spaced
    #dr = (R-Ri)/nRad
    #rs = np.linspace(Ri+dr,R,nRad)

    if Ri < 1.0e-9 :
        nRad += 1
        rs = np.insert(rs,0,(rs[0]+Ri)/2.0)    
    
    # decide whether to create a center point
    # and create the coordinate arrays
    thetai = np.pi/2.0/nAzim
    if Ri < 1.0e-9 :
        x = np.zeros(nRad*(nAzim+1)+1)
        y = np.zeros(nRad*(nAzim+1)+1)
        x[0] = 0.0
        y[0] = 0.0
        cntr = 1
    else :
        x = np.zeros(nRad*(nAzim+1))
        y = np.zeros(nRad*(nAzim+1))
        cntr = 0
    
    # mesh the annulus
    for ri in rs :
        for j in range(nAzim+1) :
            x[cntr] = ri*np.cos(j*thetai)
            y[cntr] = ri*np.sin(j*thetai)
            cntr += 1
            
    # fill in the square part
    if square :
        if nRad == 1 :
            div = R - Ri
        else :
            div = rs[nRad-1]-rs[nRad-2]
        for j in range(nAzim+1) :
            if j*thetai < np.pi/4.0 :
                dr = R/np.cos(j*thetai) - R
            else :
                dr = R/np.sin(j*thetai) - R
            if dr < 0.6*div :
                nPnts = 0
                x[(nRad-1)*(nAzim+1)+j] = (R+dr)*np.cos(j*thetai)
                y[(nRad-1)*(nAzim+1)+j] = (R+dr)*np.sin(j*thetai)
            else :
                nPnts = int(np.floor(dr/div))
            for i in range(nPnts) :
                x = np.append(x, (R+(i+1)*dr/nPnts)*np.cos(j*thetai))
                y = np.append(y, (R+(i+1)*dr/nPnts)*np.sin(j*thetai))
        if nAzim%2 == 1 :
            x = np.append(x,R)
            y = np.append(y,R)
    return x,y

def getPinCellMesh(Rf,Rc,qpitch,fuelDisc,cladDisc,modDisc) :
    (x1,y1) = annulusPoints(fuelDisc[0],fuelDisc[1],Rf,0.0)
    (x2,y2) = annulusPoints(cladDisc[0],cladDisc[1],Rc,Rf)
    (x3,y3) = annulusPoints(modDisc[0],modDisc[1],qpitch,Rc, True)
    x = np.concatenate((x1,x2,x3))
    y = np.concatenate((y1,y2,y3))
    
    tri = Delaunay(zip(x,y))
    convex_hull = [[tri.simplices[s,i] for i in range(3) if tri.neighbors[s,i] > -1]
                   for s in range(len(tri.simplices)) if -1 in tri.neighbors[s,:]]
    convex_hull = np.array(convex_hull)

    materials = np.zeros(tri.simplices.shape[0], dtype=np.int32)
    for i,s in enumerate(tri.simplices) :
        xm = (tri.points[s[0],0] + tri.points[s[1],0] + tri.points[s[2],0])/3.0
        ym = (tri.points[s[0],1] + tri.points[s[1],1] + tri.points[s[2],1])/3.0
        r = np.sqrt(xm**2 + ym**2)
        if r < Rf :
            materials[i] = 1
        elif r < Rc :
            materials[i] = 2
        else :
            materials[i] = 3
    
    return simpleMesh(tri.points,
                      tri.simplices,
                      materials,
                      tri.neighbors,
                      convex_hull)
