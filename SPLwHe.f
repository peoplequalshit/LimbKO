! systematic condition
      implicit double precision (a-h,o-z)
      integer reac
      double precision Eg,lEg,dLEg,Ep,lEp,dlEp,y,dummy,norm,gamma,gamma1
      double precision gamma2,Ebreak,mp,powerlaw,powerlaw1,powerlaw2,CHe
      double precision gammaHe,deltagammaHe,sHe,factorHe,frac,sigmafrac
      double precision dNHdR,dNHedR,dNHedR1,dNHedR2,RH,RHe,RHe0,normAll
      double precision par(5)
      character(len=70) fn
      character(len=20) :: arg

      call spec_ini                !initialization

      reac=0                       !pp collision, combined Kamae+QGSJET
      id=0                         !photon spectrum
      sigmafrac=1.77               !fraction of crssection HeN/pN

      do i=1,iargc()
       call getarg(i,arg)
       read(arg,*) par(i)
      enddo
      norm=par(1)
      gamma1=par(2)
      gamma2=par(3)
      Ebreak=par(4)
      normAll=par(5)
! my condition
      mp=0.938
      Nbinsg=50 ! actually *2
      lEming=1
      lEmaxg=3
      Nbinsp=50 ! actually *2
      lEminp=1
      lEmaxp=6
! test multiple file
!      write(fn,fmt='(i0,a)') k, '.dat'
      open(2,file='0.dat')
! Define powerlaw rigidity of Helium
      C=0.0948
      gammaHe=2.780
      deltagammaHe=0.119
      sHe=0.027
      RHe0=245
! start mycode
      dlEg=log(10.0)*(1.0/Nbinsg)
      dlEp=log(10.0)*(1.0/Nbinsp)
      do i=lEming,(Nbinsg*(lEmaxg-lEming))+1
       Eg=10.0**(lEming+((i-1.0)/Nbinsg))
       y=0
       do j=lEminp,(Nbinsp*(lEmaxp-lEminp))+1
        Ep=(10.0)**(lEminp+((j-1.0)/Nbinsp))
        if (Ep>Eg) then
         RH=sqrt(Ep*(Ep+2*mp))
         fff=spec_int(Ep,Eg,id,reac)
         powerlaw1=norm/(RH**gamma1)
         !(Ep*(Ep+(2.0*mp)))**(-0.5*gamma1)
         powerlaw2=(Ep+mp)/sqrt(Ep*(Ep+2.0*mp))
         powerlaw=powerlaw2*powerlaw1
         !!
         RHe=2.0*RH
         dNHedR1=C/((RHe/45.0)**(gammaHe))
         dNHedR2=(1.0+((RHe/RHe0)**(deltagammaHe/sHe)))**sHe
         dNHedR=dNHedR1*dNHedR2
         dNHdR=norm/(RH**gamma1)
         frac=sigmafrac*dNHedR*(1.0/dNHdR)*2.0
         factorHe=1.0+frac
         !!
         sum=(Ep/Eg)*powerlaw*fff*dlEp*factorHe
         y=y+sum
        endif
       enddo
       write(2,*) Eg,y*normAll
      enddo
      close(2)
      end
