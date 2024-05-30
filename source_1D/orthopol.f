c-----orthopol-----------------------------------------------------------
c
      subroutine orthopol(rhop,rhon,rstep,ngrid,langmo,maxorder)

      include 'params.inc'
c
c     computes orthogonal polynamials over density 'rhop'
c
      dimension rhop(kg1),rhon(kg1),rhoweight(kg1)
      dimension polyp(kg1,0:maxorder)
c
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
c------------------------------------------------------------------------------
c
c      write(21,*) '# enter ORTHPOL: rstep,ngrid,langmo,maxorder=',
c     &   rstep,ngrid,langmo,maxorder
      if(langmo.le.0) stop ' ORTHOPOL: needs L>0'
      if(maxorder.gt.10) stop  ' ORTHOPOL: too large order'
c
c     prepare weighted density
c
      do i=1,ngrid
        r = (i-1)*rstep
        rhoweight(i) = r**(2*langmo)/rhop(i)
      enddo
c
c     build polynomials
c
      do norder=0,maxorder
        do i=2,ngrid
          r = (i-1)*rstep
          if(norder.eq.0) then
            polyp(i,norder) = 1.0D0
          else
            polyp(i,norder) = r**(2*norder)
          endif
        enddo       
        do nord2=0,norder
          ovab = overlrho(rhoweight,polyp(1,norder),polyp(1,nord2),
     &                    rstep,ngrid)
          do i=2,ngrid
            polyp(i,norder) = polyp(i,norder)-ovab*polyp(i,nord2)
          enddo       
        enddo
        ovaa = overlrho(rhoweight,polyp(1,norder),polyp(1,norder),
     &                  rstep,ngrid)
        ovaa = 1.0D0/sqrt(ovaa)
        do i=2,ngrid
          polyp(i,norder) = polyp(i,norder)*ovaa
        enddo       
      enddo
c
c     print
c
      do i=2,ngrid
        r = (i-1)*rstep
        write(21,'(f6.1,11(1pg12.4))') 
     &    r,rhop(i),(polyp(i,norder)*r**langmo,norder=0,maxorder)
c     &    r,(polyp(i,norder)*sqrt(rhoweight(i)),norder=0,maxorder)
      enddo
c
      return
      end

c-----overlrho-----------------------------------------------------------
c
      function overlrho(rho,fielda,fieldb,rstep,ngrid)

      include 'params.inc'
c
c     'rho' weighted overlap between 'fielda' and 'fieldb'
c
      dimension rho(kg1),fielda(kg1),fieldb(kg1)
c
      data fourpi/12.56637061D0/
c
c------------------------------------------------------------------------
c
      sum = 0.0D0
      do i=2,ngrid
        r = (i-1)*rstep
        sum = (r*r)*rho(i)*fielda(i)*fieldb(i) + sum
      enddo
      overlrho = fourpi*rstep*sum
c
      return
      end

