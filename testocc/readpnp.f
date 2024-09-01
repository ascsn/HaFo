c---------------------------------------------------------------------
c
      program readpnp
c
      implicit none
c
c     The routine reads data conerning PNP from the SHF package.
c
c     The projection reads
c
c       X  =   Counter / Denominator
c
c       Counter     = int  dphi int dr r**2  ( rkernel *  rotrho**power )
c       Denominator = int  dphi  ( rkernel ) 
c
c     The 'rkernel' is the product of basic phase factors 
c       
c       rkernel = exp(-i*N*phi-beta*N)*(u**2+v**2*exp(beta+i*2*phi))
c
c     The 'rotrho' is the transition density.
c
      integer minshift,maxshift,maxphi,kgrid
      parameter (minshift=-60,maxshift=50)
      parameter (maxphi=71)
      parameter (kgrid=91)
c
      complex*16 rotrho(kgrid,minshift:maxshift,maxphi)
      complex*16 rkernel(minshift:maxshift,maxphi)
      complex*16 addrho,addnorm
      complex*16 accnorm,accrho,accrho3
      real*8     rgrid(kgrid)
      real*8     shift(minshift:maxshift)
      real*8     phigrid(maxphi)
      real*8     ract2,frpirs
      integer    iphi,ir,ishift,ngrid
c
c---------------------------------------------------------------------
c
c     read data from file
c
      open(17,form='unformatted',file='pnp.data')
      read(17) ngrid,frpirs
      read(17) (rgrid(ir),ir=1,ngrid)
      read(17) phigrid
      read(17) shift
      do ishift=minshift,maxshift
      do iphi=0,maxphi-1
        read(17) rkernel(ishift,iphi+1),
     &          (rotrho(ir,ishift,iphi+1),ir=1,ngrid)
      enddo
      enddo
c
c     compute projection as a test
c
      do ishift=minshift,maxshift
        accnorm = DCMPLX(0.0d0,0.0d0)
        accrho  = DCMPLX(0.0d0,0.0d0)
        accrho3 = DCMPLX(0.0d0,0.0d0)
        do iphi=0,maxphi-1
          addnorm  = rkernel(ishift,iphi+1)
          accnorm = accnorm + addnorm
          do ir=2,ngrid
            ract2 = rgrid(ir)**2
            addrho  = rotrho(ir,ishift,iphi+1)
            accrho  = accrho  + addnorm*ract2*addrho**2
            accrho3 = accrho3 + 
     &                addnorm*ract2*addrho**(7.0D0/3.0D0)
          enddo
        enddo
        write(6,'(i5,f8.4,4(1pg18.9))')
     &   ishift,shift(ishift),frpirs*accrho/accnorm,
     &   frpirs*accrho3/accnorm
      enddo
c
      stop
      end
