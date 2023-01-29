      program ap_photo_Halley_auto
c
c     Code modified for Halley
c
c     This program calculates the appearance of gaseous jets in a comet as seen
c     from Earth assuming radicals are released from a parent species. It is 
c     assumed that fluoresence from radicals as the emission source. By setting
c     the velocity of radicals and the parent lifetime to very small numbers,
c     and the radical lifetime to a sufficiently large number (characteristic
c     of lifetimes for dust grains), one should be able to mimic the dust
c     jets.
c
c     Nucleus is assumed to be spherical and active regions have finite
c     source sizes (circular in shape). In addition, it allows a gaussian 
c     dispersion normal to the radial direction.
c
c     The production rate of parent species is proportional to the solar
c     zenith angle for the sunward side and the zero production rate is
c     assumed for the anti-sunward side. In addition, one can have a 
c     component representing uniform emission over the entire nucleus 
c     (both day and night sides). By setting this "uniform background" 
c     to zero one can turn-off the uniform background emission.
c
c     Image orientation is such that north is up and east is to the left 
c     when the (0,0) pixel is at the lower left corner of the image.
c
c     WRITTEN BY NALIN SAMARASINHA 
c                    last modified ----- December 2016
c     
      common /first/idu
      character*72 imnew,errmsg
      character*8 fileout
      character*32 fileName
      integer axlen(7)
      integer numnor,numnor1,k1,k2,k3,k11
      integer ier,imn,naxis,dtype
      real lifgra,lifrad
      real array(65536),avect(256)
      real rad(10),aper(10)
      grav=6.673e-8
      solmas=1.989e+33
      au=1.496e+13
      pi=2.0*asin(1.0)
      piby2=pi/2.0
      twopi=pi*2.0
      conv=pi/180.0

      call get_command_argument(1, fileName)
c      print*,'Enter the file name to be read: '
c      read*,fileName
      
      open (1, file = fileName, status='old')
      
      imnew = fileName
      read(1,*) helio
      print*, 'helio: ', helio
      read(1,*) didist
      print*, 'didist: ', didist
      read(1,*) scale
      print*, 'scale', scale
      pix=725.3*scale*didist
      read(1,*) cenlon
      print*, 'cenlon: ', cenlon
      read(1,*) cenlat
      print*, 'cenlat: ', cenlat
      cenlon=cenlon*conv
      cenlat=cenlat*conv
      ccenla=cos(cenlat)
      scenla=sin(cenlat)
      read(1,*) raddeg
      print*, 'raddeg: ', raddeg
      radar=raddeg*conv
      cradar=cos(radar)
      read(1,*) sigjet
      print*, 'sigjet: ', sigjet
      sigjet=sigjet*conv
      constr=1.0-exp(-0.125*(pi/sigjet)**2)
      read(1,*) angra
      print*, 'angra: ', angra
      read(1,*) angdec
      print*, 'angdec: ', angdec
      angra=angra*conv
      angdec=angdec*conv
      cangde=cos(angdec)
      sangde=sin(angdec)
      read(1,*) solra
      print*, 'solra: ', solra
      read(1,*) soldec
      print*, 'soldec: ', soldec
      solra=solra*conv
      soldec=soldec*conv
      csolde=cos(soldec)
      ssolde=sin(soldec)
      xtail=-csolde*cos(solra)
      ytail=-csolde*sin(solra)
      ztail=-ssolde
      read(1,*) earra
      print*, 'earra: ', earra
      read(1,*) eardec
      print*, 'eardec: ', eardec
      earra=earra*conv
      eardec=eardec*conv
      cearde=cos(eardec)
      searde=sin(eardec)
      cearra=cos(earra)
      searra=sin(earra)
      read(1,*) numap
      print*, 'numap: ', numap
      do 100 n=1,numap
         read(1,*) rad(n)
         print*, 'rad ', n, ': ', rad(n)
c     convert aperure radii from pixels to km
         rad(n)=rad(n)*pix
 100  continue
c
      do 150 n=1,numap	
      aper(n)=0.0
 150  continue
c
c     find the PA of the equatorial cdt system north wrt.
c     ecliptic cdt system
      polelo=90.0*conv
      polela=66.55*conv
      spole=sin(polela)
      cpole=cos(polela)
      eplon=earra-polelo
      ceplon=cos(eplon)
      cospe=searde*spole+cearde*cpole*ceplon
      sinpe=sqrt(1.0-cospe**2)
      papole=acos((spole*cearde-cpole*searde*ceplon)/sinpe)
      if(eplon.lt.0.0) eplon=eplon+twopi
      if(eplon.gt.pi) papole=-papole
      spa=sin(papole)
      cpa=cos(papole)
c     papole is the negative of the PA of the ecliptic cdt system
c     north wrt. equatorial cdt system
c
      xnorth=-searde*cearra
      ynorth=-searde*searra
      znorth=cearde
      xeast=searra
      yeast=-cearra
      xyznor=xtail*xnorth+ytail*ynorth+ztail*znorth
      xyzeas=xtail*xeast+ytail*yeast
      read(1,*) angthe
      print*, 'angthe: ', angthe
      angthe=angthe*conv
      cangth=cos(angthe)
      sangth=sin(angthe)
      read(1,*) ratpsi
      print*, 'ratpsi: ', ratpsi
      ratpsi=twopi/(ratpsi*86400.0)
      read(1,*) ratphi
      print*, 'ratphi: ', ratphi
      ratphi=twopi/(ratphi*86400.0)
      read(1,*) date
      print*, 'date: ', date

      timfin=date
      read(1,*) gpsi
      print*, 'gpsi: ', gpsi
      gpsi=gpsi*conv
      read(1,*) gphi
      print*, 'gphi: ', gphi
      gphi=gphi*conv

      beta=grav*solmas/au/au

      read(1,*) radg
      print*, 'radg: ', radg
      radg=radg*(1.0e-5)/(helio**2)
      read(1,*) radr
      print*, 'radr: ', radr
      radr=radr*(1.0e-5)/(helio**2)
      read(1,*) velgr
      print*, 'velgr: ', velgr
c     print*,'give the range (2*FWHM) of parent/grain velocities in km/sec'
c     print*,'(the velocity distribution is assumed to be'
c     print*,' triangular and peaks at the mean velocity)'
c     read*,delvg
c     delvg=0.5*delvg
c      print*,'give radical velocity in km/sec'
      read(1,*) velr
      print*, 'velr: ', velr
      read(1,*) lifgra
      print*, 'lifgra: ', lifgra
      lifgra=lifgra/(helio**2)
      read(1,*) lifrad
      print*, 'lifrad: ', lifrad
      lifrad=lifrad/(helio**2)
      read(1,*) timbin
      print*, 'timbin: ', timbin
c      print*,'give the maximum time jet action to be followed in sec'
      read(1,*) totime
      print*, 'totime: ', totime
      num=nint(totime/timbin)
      read(1,*) nregion
      print*, 'nregion: ', nregion
c      print*,'give the maximum number of parent/grain particles to be'
c      print*,'emitted per sampled location (per time bin)'
c      print*,'[suggestion: > 100]'
      read(1,*) binnum
      print*, 'binnum: ', binnum
c      print*,'give the constant background production rate as a'
c      print*,'fraction of the maximum production rate'
      read(1,*) bg
      print*, 'bg: ', bg
c 
c     calculate the position of a radical as seen from the spacecraft
c
      timeob=timbin*float(num)
      fphi=gphi+ratphi*(timfin*86400.0-timeob)
      fpsi=gpsi+ratpsi*(timfin*86400.0-timeob)
c
      
      do 200 k=1,65536
      array(k)=0.0
  200 continue
      idu=-1
      
      do 800 i=1,num
c
      do 700 ii=1, nregion
c
c     calculate the cometographic longitude and latitude of the
c     parent/grain to be emitted

         azireg=twopi*rand2(idu)
         cazire=cos(azireg)
         disreg=acos(1.0-(rand2(idu))*(1.0-cradar))
         cdisre=cos(disreg)
         sdisre=sin(disreg)
         comlat=asin(scenla*cdisre+ccenla*sdisre*cazire)
         scomlo=sdisre*sin(azireg)
         ccomlo=cdisre*ccenla-sdisre*scenla*cazire
         comlon=atan2(scomlo,ccomlo)
         comlon=comlon+cenlon

c
c     calculate the body fixed coordinates of the parent/grain origin

         ccomla=cos(comlat)
         body1=ccomla*cos(comlon)
         body2=ccomla*sin(comlon)
         body3=sin(comlat)

c
         time0=timbin*float(i)
         angphi=ratphi*time0+fphi
         angpsi=ratpsi*time0+fpsi
         cangps=cos(angpsi)
         sangps=sin(angpsi)
         cangph=cos(angphi)
         sangph=sin(angphi)
         space1=body1*(cangps*cangph-sangps*sangph*cangth)
     +        -body2*(sangps*cangph+cangps*sangph*cangth)
     +        +body3*sangph*sangth

      space2=body1*(cangps*sangph+sangps*cangph*cangth)
     +       -body2*(sangps*sangph-cangps*cangph*cangth)
     +        -body3*cangph*sangth

      space3=body1*sangps*sangth+body2*cangps*sangth+body3*cangth

c
c     calculate the cometocentric latitude and longitude of the
c     point of activity in space fixed momentum frame
c

      spalat=asin(space3)
      cspala=cos(spalat)
      sspala=sin(spalat)
      spalon=atan2(space2,space1)
      cspalo=cos(spalon)
      sspalo=sin(spalon)
c
c     calculate ecliptic longitude and latitude of the point of activity 
c
      spadec=asin(sangde*sspala-cangde*cspala*cspalo)
      cspade=cos(spadec)
      sspade=sin(spadec)
      sspara=cspala*sspalo
      cspara=sspala*cangde+cspala*sangde*cspalo
      spara=atan2(sspara,cspara)
      spara=spara+angra
c
c     determine whether the point of activity is in sunlight
c
      csunan=sspade*ssolde+cspade*csolde*cos(solra-spara)
      if(csunan.lt.0.0) csunan=0.0
      csunan=csunan+bg
      n100=mod(i,100)
      if(n100.ne.0) goto 300
      timeto=time0-timeob
c     write(2,995) timeto,csunan
  300 nbin=nint(binnum*csunan)
c
c     Assume production rate is proportional to the insolation to the
c     third power
c
c     300 nbin=nint(binnum*csunan*csunan*csunan)
      disgr=totime-time0
c
c     calculate the distance travelled based on empirical data for CN
c     (for dust:  disg=2.6 + 0.031181*disgr + 0.000006697*disgr*disgr )
c     disg=2.6 + 0.5*disgr + 2.0*0.000006697*disgr*disgr
c     disg=2.6 + 0.5*disgr
c     if(disgr.lt.8000.0) then
c     disg=2.6 + 0.3*disgr
c     else
c     disg=2.6 + 0.3*8000.0 + 0.7*(disgr-8000.0)
c     endif
c     if(disgr.lt.4444.0) then
c     disg0=2.6 +  0.2*disgr + 0.00005625*disgr*disgr
c     else
c     disg0=2000.0 + 0.7*(disgr-4444.0) 
c     endif
c
c     calculate for all the parents/grains emitted in the time bin
c     under consideration
c
      do 600 j=1,nbin
c
      call prod(lifgra,timgra)

      time1=timgra+time0

      if(time1.gt.totime) goto 600

      fracr=rand2(idu)

      timrad=-lifrad*alog(fracr)

      time2=time1+timrad

      if(time2.lt.totime) goto 600

c
c     calculate the specific parent/grain velocity
c     call vran(delvg,deltav)
c     velg=velgr+deltav-delvg

      velg=velgr

      disg=disgr*velg
c     consider effects due to dispersion in radial outflow
c
      print*,'disp=gasdev(idu)'
      disp=gasdev(idu)
c     disg=disg0*(1.0+0.05*disp)
c     5 percent chosen based on Dave's input.
c

      disg=disg*(1.0+0.05*disp)
c

      azijet=twopi*rand2(idu)

      cazije=cos(azijet)

      disjet=sigjet*sqrt(-2.0*alog(1.0-constr*rand2(idu)))

      cdisje=cos(disjet)

      sdisje=sin(disjet)

      decjet=asin(sspade*cdisje+cspade*sdisje*cazije)

      cdecje=cos(decjet)

      srajet=sdisje*sin(azijet)

      crajet=cdisje*cspade-sdisje*sspade*cazije

      rajet=atan2(srajet,crajet)

      rajet=rajet+spara

c     calculate the distances travelled

      xjet=cdecje*cos(rajet)

      yjet=cdecje*sin(rajet)

      zjet=sin(decjet)

      gnorth=disg*(xjet*xnorth+yjet*ynorth+zjet*znorth)

      geast=disg*(xjet*xeast+yjet*yeast)

      time3=totime-time1

      disr=velr*time3

      disag=radg*timgra*(0.5*timgra+time3)

      disar=0.5*radr*(time3**2)

      disa=disag+disar

      radazi=twopi*rand2(idu)

      cradpr=sqrt(1.0-(rand2(idu))**2)

      rnorth=disr*cradpr*cos(radazi)

      reast=disr*cradpr*sin(radazi)

      anorth=disa*xyznor
      aeast=disa*xyzeas
      dirn=gnorth+rnorth+anorth
      dire=geast+reast+aeast
      totn=dirn*cpa+dire*spa
      tote=dire*cpa-dirn*spa
c
      totr=sqrt(totn*totn+tote*tote)
c
      do 400 n=1,numap
      if(totr.lt.rad(n)) then
      print*,'numap: ', numap
      print*,'aper #: ', n  
      aper(n)=aper(n)+1.0
      endif
 400  continue
c
c     write data into an array
      numnor=nint(128.5+totn/pix)
      if(numnor.lt.1) goto 600
      if(numnor.gt.256) goto 600
      numeas=nint(128.5-tote/pix)
      if(numeas.lt.1) goto 600

      if(numeas.gt.256) goto 600
      numnor1=numnor-1
      k=(numnor1*256)+numeas
      array(k)=array(k)+1.0
 600  continue

 700  continue

 800  continue

c     create an image
      naxis=2

      dtype=6
      axlen(1)=256
      axlen(2)=256
      do 900 kk=3,7
      axlen(kk)=1
  900 continue
c      call imcrea(imnew,axlen,naxis,dtype,ier)
      if(ier.ne.0) goto 999
c     open the image
c      call imopen(imnew,3,imn,ier)
      if(ier.ne.0) goto 999
c     read the array into the image
      do 980 k1=1,256
      do 960 k2=1,256
      k11=k1-1
      k3=(k11*256)+k2
      avect(k2)=array(k3)
  960 continue
c      call impl2r(imn,avect,k1,ier)
      if(ier.ne.0) goto 999
  980 continue
c     close the image
c      call imclos(imn,ier)
      if(ier.ne.0) goto 999
c
      do 990 n=1,numap
      aper(n)=alog10(aper(n))
  990 continue
c

      write(fileout,994)date
      open(11,file='output_'//fileout)
      write(11,996) date, (aper(n), n=1,numap)
  994 format(f6.4)
  995 format(1x,2(1pe12.4))
 996  format(1x,6(1pe13.5))
 999  stop
c     print error messages
c  999 call imemsg(ier,errmsg)
      write(*,'("error:",a80)') errmsg
      end
 
     
c     subroutine to calculate the time taken for a parent
c     to produce a daughter, assuming a linear growth and
c     linear decay (with equal growth and decay times) for
c     the daughter production rate
c     subroutine prod(palife,patime)
c     common /first/idu
c     rannum=rand2(idu)
c     rnum=rand2(idu)
c     snum=sqrt(rnum)
c     if(rannum.gt.0.5) goto 20
c     patime=palife*snum
c     goto 40
c  20 patime=palife*(2.0-snum)
c  40 return
c     end


c     subroutine to calculate the time taken for a parent
c     to produce a daughter, assuming an exponential decay
c     for the daughter production rate
      subroutine prod(palife,patime)
      common /first/idu
      rannum=rand2(idu)
      patime=-palife*alog(rannum)
      return
      end

c     subroutine to calculate the velocity differential
c     with respect to the mean parent/grain velocity
c     subroutine vran(width,veldif)
c     common /first/idu
c     arannu=rand2(idu)
c     arnum=rand2(idu)
c     asnum=sqrt(arnum)
c     if(arannu.gt.0.5) goto 60
c     veldif=width*asnum
c     goto 80
c  60 veldif=width*(2.0-asnum)
c  80 return
c     end


c     Gaussian distribution from NUMERICAL RECIPES
c     by PRESS, FLANNERY, TEUKOLSKY, and VETTERLIN
      function gasdev(idum)
      data iset/0/
      if (iset.eq.0) then
        
    1   v1=2.*rand2(idum)-1.
        v2=2.*rand2(idum)-1.
        r=v1**2+v2**2

        if(r.ge.1.)go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
        
      else
        
        gasdev=gset
        iset=0
        
      endif
      return
      end


c     random number generating function from NUMERICAL RECIPES
c     by PRESS, FLANNERY, TEUKOLSKY, and VETTERLING
c     (ran2 renamed rand2)
      function rand2(idum)
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      dimension ir(97)
      data iff /0/
      if(idum.lt.0.or.iff.eq.0) then
          iff=1
          idum=mod(ic-idum,m)
          do 11 j=1,97
              idum=mod(ia*idum+ic,m)
              ir(j)=idum
   11     continue
          idum=mod(ia*idum+ic,m)
          iy=idum
      endif
      j=1+(97*iy)/m
c      if(j.gt.97.or.j.lt.1) pause
      iy=ir(j)
      rand2=iy*rm
      if(rand2.lt.1.0e-10) rand2=1.0e-10
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
