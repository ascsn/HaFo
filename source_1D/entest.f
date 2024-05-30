      if(testen) then
        sumt0 = 0.0
        sumt1 = 0.0
        sumls = 0.0
        sumkin = 0.0
        do i=2,ngrid
          j      = i+ngrid
          r      = (i-1)*rstep
          vol    = r*r
          rp     = rho(i)
          rm     = rho(j)
          rs     = rp+rm
          tp     = tau(i)
          tm     = tau(j)
          ts     = tp+tm
          ssp    = sjs(i)
          ssm    = sjs(j)
          sss    = ssp+ssm
          sp     = sj(i)
          sm     = sj(j)
          ss     = sp+sm
          d1p    = d1r(i)
          d1m    = d1r(j)
          d1s    = d1p+d1m
          d2p    = d2r(i)
          d2m    = d2r(j)
          d2s    = d2p+d2m
          sumt0  = vol*(t0a*rs*rs-t0b*(rp*rp+rn*rn))     + sumt0  
          sumt1  = vol*(t1p2*ts*rs-t1m2*(tp*rp+tn*rn)
     &                  -t13m2*d2s*d2s+t13p2*(d2p*d2p+d2n*d2n))
     &            + sumt1
          sumls  = vol*( -b4*sss*sss- b4p*(ssp*ssp+ssn*ssn)) + sumls
          sumkin = vol*(amfc(1)*tp+amfc(2)*tn)    + sumkin
        enddo
        sumt0  = frpirs*sumt0*0.5
        sumt1  = frpirs*sumt1
        sumls  = frpirs*sumls
        sumkin = frpirs*sumkin
        write(6,'(a,4g13.5)')
     &    ' epot_0,epot_1,epot_ls,ekin=',sumt0,sumt1,sumls,sumkin
      endif
