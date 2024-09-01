c-----mcsear ---------------------------------------------------------
c
      subroutine mcsear(zuord,a,sigmaa,sigma0,nterm1,nterm2,
     *                  nwdr1,nwdr2,ifcurv,chitmp)
c
c     Driver routine for Monte-Carlo search.
c     Calls bevington step for iteration of first 'nterm1'
c     parameters.
c     The subsequent 'nterm2' parameters are minimised by
c     Monte-Carlo sampling.
c     
c
c      the parameters in the calling sequence are:
c       nterm1  = number of varied parameters in Bevington step
c       nterm2  = number of varied parameters in outer loop
c                 these occupy the places 'nterm1'+1...'nterm1'+'nterm2'
c       nwdr1   = maximum number of iterations in inner Bevington step
c       nwdr2   = maximum number of iterations in outer optimisation
c       ifcurv  = number of normal bevington steps which are to be
c                 run before the curved step is switched on. if no
c                 curved step is desired make ifcurve larger than
c                 nwdr.
c       zuord   = names of parameters for printing of protocol 'for003'
c       a       = initial parameters as input, new parameters after
c                 one iteration step given back as output.
c       sigmaa   = field containing the correlated errors on the
c                  model parameters, it is a 1-dimensional field
c                  sigmaa(parameter).
c       sigma0   = field containing the uncorrelated errors.
c
c
c     dimension statement valid for nterms up to 'kparmx'
c
      include 'params.inc'
      include 'comfit.inc'
c
      character*6 zuord(kparmx)
      dimension a(kparmx),sigmaa(kparmx),sigma0(kparmx)
c
      dimension avar(kparmx)
c
cold      include 'numerc.inc'
cold      include 'shells.inc'
cold      include 'forces.inc'
cold      include 'result.inc'
c
      data faclim/4.0/                      ! factor limiting outer step
      data iseed/19067/
cold      data chitmp/1.0/                      ! pseudo-temperature in chi**2
c
c
c--------------------------------------------------------------------
c
c     chi**2 at initial parameter set
c
      write(6,'(a,g12.4)') 
     &   ' MCSEAR entered: chi**2-temperature=',chitmp
      open(4,file='for004')
      write(4,'(a,g12.4)') 
     &   ' MCSEAR entered: chi**2-temperature=',chitmp
c
      call beving(zuord,a,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
      chi00 = chisqr
      write(2,'(//40a2/a,3g12.4)')
     &     ('* ',i=1,40),
     &     ' first BEVING: chisqr,chiq,chimin=',chisqr,chiq,chimin
      write(4,'(a,3g14.6)')
     &     ' first BEVING: chisqr,chiq,chimin=',chisqr,chiq,chimin
c
c     MC loop
c
      do iouter=1,nwdr2
        write(2,'(/a,g13.5)')
     &  ' next step in MC loop starting with chi**2=',chi00
        call flush(4)
c        call spool(4,'for004')
        write(4,'(/a,g13.5)')
     &  ' next step in MC loop starting with chi**2=',chi00
c
c       stochastic variation of parameters
c
        do i=1,kparmx
          avar(i) = a(i)
        enddo
        do npar2=nterm1+1,nterm1+nterm2
          stppar = gasdev(iseed)*dparmo(npar2)
          stpmax  = faclim*abs(dparmo(npar2))
          stppar  = sign(min(stpmax,abs(stppar)),stppar)
          avar(npar2) = a(npar2)+stppar
        enddo
        call flush(4)
c          call spool(4,'for004')
        call beving(zuord,avar,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
        delchi = chisqr - chi00
        write(4,'(/a,i3,a,2g14.6,10(/5x,4g14.6))')
     &    ' after trial step',iouter,',: chi**2,delta chi**2, params=',
     &    chisqr,delchi,(avar(i),i=1,nterm1+nterm2)
        write(2,'(//40a2/a,i3,a,2g14.6,10(/5x,4g14.6))')
     &    ('* ',i=1,40),
     &    ' after trial step',iouter,',: chi**2,delta chi**2, params=',
     &    chisqr,delchi,(avar(i),i=1,nterm1+nterm2)
        if(delchi.lt.0.0D0 .or. exp(-delchi/chitmp).gt.ran1(iseed)) then
          do i=1,kparmx
            a(i) = avar(i)
          enddo
          chi00 = chisqr
          write(4,'(a)') ' new parameters accepted'
        endif
      enddo
c
      return
      end
c-----splitb ---------------------------------------------------------
c
      subroutine splitb(zuord,a,sigmaa,sigma0,nterm1,nterm2,
     *                  nwdr1,nwdr2,ifcurv,jpact)
c
c     Driver routine for split bevington step.
c     Calls bevington step for iteration of first 'nterm1'
c     parameters.
c     The subsequent 'nterm2' parameters are minimised by
c     fully quadratic gradient step always calling the
c     re-optimisation of the 'nterm1' first parameters.
c     
c
c      the parameters in the calling sequence are:
c       nterm1  = number of varied parameters in Bevington step
c       nterm2  = number of varied parameters in outer loop
c                 these occupy the places 'nterm1'+1...'nterm1'+'nterm2'
c       nwdr1   = maximum number of iterations in inner Bevington step
c       nwdr2   = maximum number of iterations in outer optimisation
c       ifcurv  = number of normal bevington steps which are to be
c                 run before the curved step is switched on. if no
c                 curved step is desired make ifcurve larger than
c                 nwdr.
c       zuord   = names of parameters for printing of protocol 'for003'
c       a       = initial parameters as input, new parameters after
c                 one iteration step given back as output.
c       sigmaa   = field containing the correlated errors on the
c                  model parameters, it is a 1-dimensional field
c                  sigmaa(parameter).
c       sigma0   = field containing the uncorrelated errors.
c
c
c     dimension statement valid for nterms up to 'kparmx'
c
      include 'params.inc'
      include 'comfit.inc'
c
      character*6 zuord(kparmx)
      dimension a(kparmx),sigmaa(kparmx),sigma0(kparmx)
c
      dimension chip(kparmx),chim(kparmx),avar(kparmx),astep(kparmx),
     &          aparp(kparmx,kparmx),aparm(kparmx,kparmx),
     &          jpact(kparmx)
c
cold      include 'numerc.inc'
cold      include 'shells.inc'
cold      include 'forces.inc'
cold      include 'result.inc'
c
      data faclim/2.0/                      ! factor limiting outer step
      data ntry/4/
      data ifauto/0/                        ! rescale DPARMO automatially
      data astep/KPARMX*0.D0/
c
c
c--------------------------------------------------------------------
c
c     chi**2 at initial parameter set
c
      open(4,file='for004')
      call beving(zuord,a,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
      chi00 = chisqr
      write(2,'(//40a2/a,3g12.4)')
     &     ('* ',i=1,40),
     &     ' first BEVING: chisqr,chimin=',chisqr,chimin
      write(4,'(a,3g14.6)')
     &     ' first BEVING: chisqr,chimin=',chisqr,chimin
c
c     outer loop
c
      do iouter=1,nwdr2
      write(2,'(/a,g13.5)')
     &  ' start outer loop with chi**2=',chi00
      call flush(4)
c      call spool(4,'for004')
      write(4,'(/a,i4,g14.6)')
     &  ' start outer loop with iouter,chi**2=',iouter,chi00
c
c       variation of outer parameters
c       maintain sum of spin-orbit
c
        sumlsparams = a(jpact(5))+a(jpact(15))              ! new
        do npar2=nterm1+1,nterm1+nterm2
c
          do i=1,kparmx
            avar(i) = a(i)
          enddo
          avar(npar2) = a(npar2)+dparmo(npar2)
          avar(jpact(5)) = sumlsparams - avar(jpact(15))    ! new
          call flush(4)
c          call spool(4,'for004')
          call beving(zuord,avar,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
          chip(npar2) = chisqr
          do i=1,kparmx
            aparp(i,npar2) = avar(i)
          enddo
          write(4,'(/a,i3,a,g14.6,10(/5x,4g14.6))')
     &     ' after vary + par#',npar2,',: chi**2, params=',
     &     chip(npar2),(avar(i),i=1,nterm1+nterm2)
          write(2,'(//40a2/a,i3,a,22g13.5)')
     &     ('* ',i=1,40),
     &     ' vary + par#',npar2,',: chi**2, params=',
     &     chip(npar2),(avar(i),i=1,nterm1+nterm2)
c
          do i=1,kparmx
            avar(i) = a(i)
          enddo
          avar(npar2) = a(npar2)-dparmo(npar2)
          call flush(4)
c          call spool(4,'for004')
          call beving(zuord,avar,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
          chim(npar2) = chisqr
          do i=1,kparmx
            aparm(i,npar2) = avar(i)
          enddo
          write(4,'(/a,i3,a,g14.6,10(/5x,4g14.6))')
     &     ' after vary - par#',npar2,',: chi**2, params=',
     &     chim(npar2),(avar(i),i=1,nterm1+nterm2)
          write(2,'(//40a2/a,i3,a,22g13.5)')
     &     ('* ',i=1,40),
     &     ' vary - par#',npar2,',: chi**2, params=',
     &     chim(npar2),(avar(i),i=1,nterm1+nterm2)
        enddo
c
c       compute gradient to next set of outer parameters
c
        call flush(4)
c        call spool(4,'for004')
        write(4,'(//40a2//a/a)')
     &     ('* ',i=1,40),
     &     ' Protocol of large step:',
     &     '  npar2  chi+      chi-      dchi   d2chi    step'
        do npar2=nterm1+1,nterm1+nterm2
          dchi    = (chip(npar2)-chim(npar2))/(2.0*dparmo(npar2))
          d2chi   = (chip(npar2)-2.0*chi00+chim(npar2))/dparmo(npar2)**2b
          stppar  = -dchi/abs(d2chi)
          stpmax  = faclim*abs(dparmo(npar2))
          stppar  = sign(min(stpmax,abs(stppar)),stppar)
          astep(npar2) = stppar
c???          avar(npar2) = a(npar2)+stppar
          write(4,'(1x,i4,5g13.5)')
     &       npar2,chip(npar2),chim(npar2),dchi,d2chi,stppar
          do i=1,kparmx
            aparp(i,npar2) = (aparp(i,npar2)-aparm(i,npar2))/
     &                       (2.0*dparmo(npar2))
          enddo
          if(ifauto.eq.1 .and. d2chi.gt.0.0D0) then
            dparmo(npar2) = freep/sqrt(d2chi)
            write(4,'(a,i2,a,1pg13.5)')
     &        ' DPARMO(',npar2,') reset to:',dparmo(npar2)
          endif
        enddo
c
c       trial shots with possible reduction of steps
c
        do itry=0,ntry
          do i=1,kparmx
            avar(i) = a(i)
          enddo
          do npar2=nterm1+1,nterm1+nterm2
            stppar = astep(npar2)*0.5D0**itry
            avar(npar2) = a(npar2)+stppar
            do i=1,nterm1
              avar(i) = avar(i) + stppar*aparp(i,npar2)
            enddo
          enddo        
          call beving(zuord,avar,sigmaa,sigma0,nterm1,nwdr1,ifcurv)
          call flush(4)
c          call spool(4,'for004')
          write(4,'(/a,i3,a,g14.6,10(/5x,4g14.6))')
     &     ' trial shot #',itry,',: chi**2, params=',
     &     chisqr,(avar(i),i=1,nterm1+nterm2)
          write(2,'(//40a2/a,i2,a,i3,a,22g13.5)')
     &     ('* ',i=1,40),
     &     ' trial',itry,' par#',npar2,',: chi**2, params=',
     &     chisqr,(avar(i),i=1,nterm1+nterm2)
          if(chisqr.le.chi00) then
            chi00 = chisqr
            do i=1,kparmx
              a(i) = avar(i)
            enddo
            goto 99
          endif
        enddo
        write(4,*) ' SPLITB: no improvement in outer variation'
        stop ' SPLITB: no improvement in outer variation'
 99     continue
c        
      enddo
c
      return
      end
c-----gasdev------------ part of tdks package -------------------------------
c
      function gasdev(idum)
c
c     produces a Gaussian distribution with unit deviate,
c     i.e. a distribution proportional to exp(-x**2/2),
c     see 'press et al, numerical recipes', sect. 7.2.
c
c     'idum' is the seed for the random generator
c     the result is returned as function value 'gasdev'
c
      include 'params.inc'
c
      data iset/0/
      save gset
c
c-----------------------------------------------------------------------
c
ctest      write(6,'(a,i10)') ' GASDEV: iseed=',idum
      if(iset.eq.0) then
 1      continue
        v1     = 2.0*ran1(idum)-1.0
        v2     = 2.0*ran1(idum)-1.0
        r      = v1*v1+v2*v2
        if(r.ge.1.0D0 .or. r.eq.0.0D0) goto 1
        fac    = sqrt(-2.0*log(r)/r)
        gset   = v1*fac
        gasdev = v2*fac
        iset   = 1
      else
        gasdev = gset
        iset   = 0
      endif
c
      return
      end
c-----spool ------------------------------------------------------------
c
      subroutine spool(iunit,fname)
c
c     closes and reopens unit 'iunit', spool to end of file.
c     routine serves to override buffering.
c
      character*(*) fname
c
c-----------------------------------------------------------------------
c
      write(6,'(a,i4,1x,a)')
     &   ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SPOOL ',
     &   iunit,fname
      close(iunit)
      open(iunit,file=fname)
      do i1=1,999
      do i2=1,999
        read(iunit,'(1x)',err=99,end=99)
      enddo
      enddo
      stop ' SPOOL: file count exhausted'
c
 99   return
      end


