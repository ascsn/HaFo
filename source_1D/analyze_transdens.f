      frpirs = fourpi*rstep
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926D0/rleng
      delql  = 0.5D0*delq
      do iq=1,ngrid
        q1     = iq*delql-delql
        if(q1.eq.zero) q1 = delql/10D0
        r      = zero
        sumd   = zero
        sumr   = zero
        sump   = zero
        sumn   = zero
        do i=2,ngrid
          r      = rstep+r
          vol    = r*r*bessel(r*q1,lrpa)
          sump   = vol*prho(i)+sump
          sumn   = vol*prho(i+ngrid)+sumn
          sumd   = vol*dqcur(i)+sumd
          sumr   = vol*rqcur(i,1)+sumr
        enddo
        sump   = frpirs*sump
        sumn   = frpirs*sumn
        sumd   = frpirs*sumd
        sumr   = frpirs*sumr
        forlon = sumd/q1
        fortra = sumr/q1
        write(6,'(1x,i3,f7.4,5g13.4)')
     &       iq-1,q1,sump,sumd,sumr,forlon,fortra
        write(14,'(1x,i3,f7.4,5(1pg20.12))')
     &       iq-1,q1,sump,sumn,sump+sumn,facpr*sump+facne*sumn
        write(293,'(1x,f7.4,5(1pg13.4))')
     &       q1,sump,sumn,sump+sumn,facpr*sump+facne*sumn
      enddo
      stop
      end
c-----bessel------------------------------------------------------
c
      function bessel(x,l)
c
      include 'params.inc'
c
c     spherical bessel function of angular momentum 'l' at
c     argument 'x'.
c
c     uses series expansion for x<l=1 and asymptotic expression
c     for x>l-1.
c
c     precision is tested only up to 5 digits.
c
      include 'constn.inc'
c
      data xargmx/1.0D10/
c
c------------------------------------------------------------------------------
c
cc       if(l.gt.7) stop ' l>7 in "bessel"'
c
c      case of zero argument
c
       if(x.eq.zero) then
        if(l.eq.0) then
         b      = one
        else
         b      = zero
        endif
c
c      very large argument  --> set result to zero
c
       elseif(abs(x).gt.xargmx) then
         b      = zero
c
c      small arguments  --> use series expansion
c
       elseif(abs(x).le.l-one) then
        add   = one
        do i=1,l
          add   = x/(2*i+1)*add
        enddo
        arg2  = -0.25*x*x
        b     = add
        do i=1,l+1
          add   = arg2*add/(i*(l+i+0.5D0))
          b     = add + b
        enddo
c
c      large arguments  --> use expansion in sin and cos
c
       else
        s      = sin(x)/x
        c      = cos(x)/x
        if(l.eq.0) then
         b      = s
        elseif(l.eq.1) then
         b      = s/x-c
        elseif(l.ge.2) then
         xinv   = one/x
         bjm    = s
         b      = s/x-c
         do lrec=2,l
           bjp    = ((lrec+lrec)*xinv-xinv)*b-bjm
           bjm    = b
           b      = bjp
         enddo
         b      = b
        endif
       endif
c
       bessel = b
       return
       end
