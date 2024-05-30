      c0 =-477.170
      c3 =-2783.68
      c0p =-718.617
      c3p =-1836.30
      power = 0.25
c
         rhonm = 0.159
         powfac = 3.0/rhonm**power
         b0     = 0.5*(c0+c3)
         b3     = powfac*(c0-c3)
         b0p    = 0.5*(c0p+c3p)
         b3p    = powfac*(c0p-c3p)
c
         rhonm = 0.160
        powfac = rhonm**power/6.0
        c0     = b0+powfac*b3
        c3     = b0-powfac*b3
        c0p    = b0p+powfac*b3p
        c3p    = b0p-powfac*b3p
c
        write(6,'(1pg18.8,a,1pg14.4,a)')
     &   c0,',',' 0.03',', 1, c0',
     &   c0p,',',' 3.0',', 6, c0p',
     &   c3,',',' 1.6',', 4, c3',
     &   c3p,',',' 60.0',', 9, c3p'
c
       stop
       end

