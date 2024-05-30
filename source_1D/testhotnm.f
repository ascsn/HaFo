c-----hotnm --------------------------------------------------------
c
      program testhotnm
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
c
c
c
      data f1d3  /0.333333333/
      data itmax/50/
      data delfp/0.01/
      data epsfp/1.0E-5/
      data ipgrid/4000/
      data delp/0.001/
      data pi/3.1415926D0/
c
c-----------------------------------------------------------------------
c
      amfca0 = 20.75D0
      sumk   = 0.0D0
      sum0   = 0.0D0
      do itemp=0,20,5
      do ie=0,20
        temper = max(itemp*1.0D0,1.0D-6)
        eferm = ie*1.0D0
        do ip=1,ipgrid
          pact2  = ip*delp
          pact2  = pact2*pact2
          vol    = pact2/(1.0+exp((amfca0*pact2-eferm)/temper))
          sum0   = vol       + sum0
          sumk   = vol*pact2 + sumk
        enddo
        sum0   = delp*sum0/(pi*pi)
        sumk   = delp*sumk/(pi*pi)
        taunmp = sumk
        write(6,'(1x,5(1pg12.4))')
     &   temper,sum0,sumk,sumk/sum0,sumk/sum0*amfca0
      enddo
      enddo
c
      return
      end
