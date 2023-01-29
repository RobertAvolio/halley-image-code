import sys
import math
import random
'''
Code modified for Halley

     This program calculates the appearance of gaseous jets in a comet as seen
     from Earth assuming radicals are released from a parent species. It is 
     assumed that fluoresence from radicals as the emission source. By setting
     the velocity of radicals and the parent lifetime to very small numbers,
     and the radical lifetime to a sufficiently large number (characteristic
     of lifetimes for dust grains), one should be able to mimic the dust
     jets.

     Nucleus is assumed to be spherical and active regions have finite
     source sizes (circular in shape). In addition, it allows a gaussian 
     dispersion normal to the radial direction.

     The production rate of parent species is proportional to the solar
     zenith angle for the sunward side and the zero production rate is
     assumed for the anti-sunward side. In addition, one can have a 
     component representing uniform emission over the entire nucleus 
     (both day and night sides). By setting this "uniform background" 
     to zero one can turn-off the uniform background emission.

     Image orientation is such that north is up and east is to the left 
     when the (0,0) pixel is at the lower left corner of the image.

     WRITTEN BY NALIN SAMARASINHA 
                    last modified ----- January 2023
'''
def main():
    
    random.seed(-1)
    
    grav=6.673e-8
    solmas=1.989e+33
    au=1.496e+13
    pi=2.0*math.asin(1.0)
    piby2=pi/2.0
    twopi=pi*2.0
    conv=pi/180.0
    rad = []
    aper = []
    axlen = []
    array = []
    avect = []
    
    if(len(sys.argv) == 1):
        print("Please input a file name.")
        return
    # Reading in file:
    
    fileName = sys.argv[1]
    f = open(fileName, "r")

    helio = float(f.readline())
    didist = float(f.readline())
    scale = float(f.readline())
    cenlon = float(f.readline())
    cenlat = float(f.readline())
    raddeg = float(f.readline())
    sigjet = float(f.readline())
    angra = float(f.readline())
    angdec = float(f.readline())
    solra = float(f.readline())
    soldec = float(f.readline())
    earra = float(f.readline())
    eardec = float(f.readline())
    numap = int(f.readline())
    for i in range(numap):
        rad.append(float(f.readline()))
    angthe = float(f.readline())
    ratpsi = float(f.readline())
    ratphi = float(f.readline())
    date = float(f.readline())
    gpsi = float(f.readline())
    gphi = float(f.readline())
    radg = float(f.readline())
    radr = float(f.readline())
    velgr = float(f.readline())
    velr = float(f.readline())
    lifgra = float(f.readline())
    lifrad = float(f.readline())
    timbin = float(f.readline())
    totime = float(f.readline())
    nregion = int(f.readline())
    binnum = float(f.readline())
    bg = float(f.readline())
    '''
    print(helio)
    print (didist)
    print (scale)
    print (cenlon)
    print (cenlat)
    print (raddeg)
    print (sigjet)
    print (angra)
    print (angdec)
    print (solra)
    print (soldec)
    print (earra)
    print (eardec)
    print (numap)
    for i in range(len(rad)):
        print (rad[i])
    print (angthe)
    print (ratpsi)
    print (ratphi)
    print (date)
    print (gpsi)
    print (gphi)
    print (radg)
    print (radr)
    print (velgr)
    print (velr)
    print (lifgra)
    print (lifrad)
    print (timbin)
    print (totime)
    print (nregion)
    print (binnum)
    print (bg)
    '''

    # Calculations:
    pix=725.3*scale*didist

    cenlon*=conv
    cenlat*=conv
    ccenla=math.cos(cenlat)
    scenla=math.sin(cenlat)

    radar=raddeg*conv
    cradar=math.cos(radar)
      
    sigjet*=conv
    constr=1.0-math.exp(-0.125*(pi/sigjet)**2)

    angra*=conv
    angdec*=conv
    cangde=math.cos(angdec)
    sangde=math.sin(angdec)

    solra*=conv
    soldec*=conv
    csolde=math.cos(soldec)
    ssolde=math.sin(soldec)
    xtail=-csolde*math.cos(solra)
    ytail=-csolde*math.sin(solra)
    ztail=-ssolde

    earra*=conv
    eardec*=conv
    cearde=math.cos(eardec)
    searde=math.sin(eardec)
    cearra=math.cos(earra)
    searra=math.sin(earra)
    for i in range(numap):
        rad[i]*=pix
        aper.append(0)
        
    # find the PA of the equatorial cdt system north wrt.
    # ecliptic cdt system
    polelo=90.0*conv
    polela=66.55*conv
    spole=math.sin(polela)
    cpole=math.cos(polela)
    eplon=earra-polelo
    ceplon=math.cos(eplon)
    cospe=searde*spole+cearde*cpole*ceplon
    sinpe=math.sqrt(1.0-cospe**2)
    papole=math.acos((spole*cearde-cpole*searde*ceplon)/sinpe)
    if(eplon < 0.0):
        eplon=eplon+twopi
    if(eplon > pi):
        papole=-papole
    spa=math.sin(papole)
    cpa=math.cos(papole)

    #
    # papole is the negative of the PA of the ecliptic cdt system
    # north wrt. equatorial cdt system
    #

    xnorth=-searde*cearra
    ynorth=-searde*searra
    znorth=cearde
    xeast=searra
    yeast=-cearra
    xyznor=xtail*xnorth+ytail*ynorth+ztail*znorth
    xyzeas=xtail*xeast+ytail*yeast

    angthe=angthe*conv
    cangth=math.cos(angthe)
    sangth=math.sin(angthe)

    ratphi=twopi/(ratphi*86400.0)

    gpsi*=conv
    gphi*=conv

    beta=grav*solmas/au/au
     
    radg*=(1.0e-5)/(helio**2)
    radr*=(1.0e-5)/(helio**2)

    lifgra/=(helio**2)
    lifrad/=(helio**2)

    #
    # calculate the position of a radical as seen from the spacecraft
    #
    
    num=int(totime/timbin)
    timeob=timbin*float(num)
    
    fphi=gphi+ratphi*(date*86400.0-timeob)
    fpsi=gpsi+ratpsi*(date*86400.0-timeob)

    for k in range(65536):
        array.append(0)
    
    for i in range(num):
        for j in range(nregion):
            
            # calculate the cometographic longitude and latitude of the
            # parent/grain to be emitted
            azireg=twopi*random.uniform(0,1)
            cazire=math.cos(azireg)
            disreg=math.acos(1.0-(random.uniform(0,1))*(1.0-cradar))
            cdisre=math.cos(disreg)
            sdisre=math.sin(disreg)
            comlat=math.asin(scenla*cdisre+ccenla*sdisre*cazire)
            scomlo=sdisre*math.sin(azireg)
            ccomlo=cdisre*ccenla-sdisre*scenla*cazire
            comlon=math.atan2(scomlo,ccomlo)
            comlon=comlon+cenlon

            # calculate the body fixed coordinates of the parent/grain origin
            ccomla=math.cos(comlat)
            body1=ccomla*math.cos(comlon)
            body2=ccomla*math.sin(comlon)
            body3=math.sin(comlat)

            time0=timbin*float(i)
            angphi=ratphi*time0+fphi
            angpsi=ratpsi*time0+fpsi
            cangps=math.cos(angpsi)
            sangps=math.sin(angpsi)
            cangph=math.cos(angphi)
            sangph=math.sin(angphi)

            space1=body1*(cangps*cangph-sangps*sangph*cangth)-body2*(sangps*cangph+cangps*sangph*cangth)+body3*sangph*sangth
            space2=body1*(cangps*sangph+sangps*cangph*cangth)-body2*(sangps*sangph-cangps*cangph*cangth)-body2*(sangps*sangph-cangps*cangph*cangth)           
            space3=body1*sangps*sangth+body2*cangps*sangth+body3*cangth

            #
            # calculate the cometocentric latitude and longitude of the
            # point of activity in space fixed momentum frame
            #

            spalat=math.asin(space3)
            cspala=math.cos(spalat)
            sspala=math.sin(spalat)
            spalon=math.atan2(space2,space1)
            cspalo=math.cos(spalon)
            sspalo=math.sin(spalon)

            #
            # calculate ecliptic longitude and latitude of the point of activity
            #
            spadec=math.asin(sangde*sspala-cangde*cspala*cspalo)
            cspade=math.cos(spadec)
            sspade=math.sin(spadec)
            sspara=cspala*sspalo
            cspara=sspala*cangde+cspala*sangde*cspalo
            spara=math.atan2(sspara,cspara)
            spara=spara+angra

            #
            # determine whether the point of activity is in sunlight
            #

            csunan=sspade*ssolde+cspade*csolde*math.cos(solra-spara)
            if(csunan < 0.0):
                csunan=0.0
            csunan=csunan+bg
            n100=i%100
            if(n100 == 0):
                timeto=time0-timeob

            nbin=int(binnum*csunan)

            #
            # Assume production rate is proportional to the insolation to the third power
            #

            disgr=totime-time0

            #
            # calculate for all the parents/grains emitted in the time bin under consideration
            #
            
            for j in range(nbin):
                timgra=-lifgra*math.log(random.uniform(0,1)) # replaced prod() function, only used once and could be condensed
                parentToDaughter=timgra+time0
                #print("Check continue 1, parentToDaughter totime: " + str(parentToDaughter) + " " + str(totime))
                if(parentToDaughter>totime):
                    continue
                
                fracr=random.uniform(0,1)
                timrad=-lifrad*math.log(fracr)
                daughterNotDecayed=parentToDaughter+timrad
                #print("Check continue 2, daughterNotDecayed totime: " + str(daughterNotDecayed) + " " + str(totime))
                if(daughterNotDecayed < totime):
                    continue
                print("Past continues")
                #
                # calculate the specific parent/grain velocity
                #

                velg=velgr
                disg=disgr*velg

                #
                # consider effects due to dispersion in radial outflow
                #
                disp=gasdev()

                #
                # disg=disg0*(1.0+0.05*disp)
                # 5 percent chosen based on Dave's input.
                #

                disg=disg*(1.0+0.05*disp)

                azijet=twopi*random.uniform(0,1)

                cazire=math.cos(azireg)

                disjet=sigjet*sqrt(-2.0*math.log(1.0-constr*random.uniform(0,1)))

                cdisje=math.cos(disjet)

                sdisje=math.sin(disjet)

                decjet=math.asin(sspade*cdisje+cspade*sdisje*cazije)

                cdecje=math.cos(decjet)

                srajet=sdisje*math.sin(azijet)

                crajet=cdisje*cspade-sdisje*sspade*cazije

                rajet=math.atan2(srajet,crajet)

                rajet=rajet+spara

                #
                # calculate the distances travelled
                #
                
                xjet=cdecje*cos(rajet)

                yjet=cdecje*math.sin(rajet)

                zjet=math.sin(decjet)

                gnorth=disg*(xjet*xnorth+yjet*ynorth+zjet*znorth)

                geast=disg*(xjet*xeast+yjet*yeast)

                time3=totime-parentToDaughter

                disr=velr*time3

                disag=radg*timgra*(0.5*timgra+time3)

                disar=0.5*radr*(time3**2)

                disa=disag+disar

                radazi=twopi*random.uniform(0,1)

                cradpr=math.sqrt(1.0-(random.uniform(0,1))**2)

                rnorth=disr*cradpr*math.cos(radazi)

                reast=disr*cradpr*math.sin(radazi)

                anorth=disa*xyznor
                aeast=disa*xyzeas
                dirn=gnorth+rnorth+anorth
                dire=geast+reast+aeast
                totn=dirn*cpa+dire*spa
                tote=dire*cpa-dirn*spa

                totr=math.sqrt(totn*totn+tote*tote)
                
                for n in range(numap):
                    if(totr < rad[n]):
                        aper[n]=aper[n]+1.0

                #
                # write data into an array
                #
                
                numnor=int(128.5+totn/pix)
                if(numnor >= 1 and numnor <= 256):
                    numeas=nint(128.5-tote/pix)
                    if(numeas >= 1 and numeas <= 256):
                        numnor1=numnor-1
                        k=(numnor1*256)+numeas
                        array[k]=array[k]+1.0


            # End of for loop
        # End of for loop
    # End of for loop

    #
    # create an image
    #

    naxis=2

    dtype=6
    axlen.append(256)
    axlen.append(256)

    for kk in range(2,7):
        axlen.append(1)
    # Removed, ier is never defined except through commented out imcrea call. 
#    if(ier != 0):
#        return
    #
    # read the array into the image
    #

    pixelTable = open(fileName + "_pix_table", "w")
    
    for row in range(256):
        if(row != 0):
            pixelTable.write("\n")
        for col in range(256):
            pixelTable.write(str(array[row*256 + col]) + " ")
            

    # Removed, ier is never defined except through commented out imcrea call.
#        if(ier != 0):
#            return

    out = open("output_" + str(date), "w")
    
    
    out.write(str(date))

    for n in range(numap):
        if(aper[n] == 0):
            out.write("999")
        else:
            aper[n]=math.log(aper[n], 10)
            out.write(str(aper[n]))  # moved to for loop for simplicity
        
    f.close()

    return
        

# Sampling Gaussian Distribution
def gasdev():
    
    r = 1
    
    while(r >= 1):
        v1 = random.uniform(-1,1)
        v2 = random.uniform(-1,1)
        r=v1**2+v2**2
        
    
    fac=math.sqrt(-2*math.log(r)/r)
    gset=v1*fac
    gasdev=v2*fac
    
    return gasdev
    
    
if __name__ == '__main__':
    main()
    
