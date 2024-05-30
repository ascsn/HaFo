From reinhard@ect.unitn.it  Fri Mar 27 19:25:06 1998
Return-Path: <reinhard@ect.unitn.it>
Received: from faui45.informatik.uni-erlangen.de (root@faui45.informatik.uni-erlangen.de [131.188.2.45])
	by faupt212.physik.uni-erlangen.de (8.8.8/8.8.8) with ESMTP id TAA04601
	for <reinhard@faupt212.physik.uni-erlangen.de>; Fri, 27 Mar 1998 19:25:06 +0100
Received: from ECTstar.ect.unitn.it (ectstar.ect.unitn.it [193.205.206.2])
	by faui45.informatik.uni-erlangen.de (8.8.8/8.1.16-FAU) with SMTP id TAA06605
	for <reinhard@faupt212.physik.uni-erlangen.de>; Fri, 27 Mar 1998 19:25:03 +0100 (MET)
Received: (from reinhard@localhost) by ECTstar.ect.unitn.it (8.6.12/8.6.12) id TAA11748 for reinhard@theorie2.physik.uni-erlangen.de; Fri, 27 Mar 1998 19:24:21 +0100
From: Paul-Gerhard Reinhard <reinhard@ect.unitn.it>
Message-Id: <199803271824.TAA11748@ECTstar.ect.unitn.it>
Subject: eosn.f
To: reinhard@faupt212.physik.uni-erlangen.de (pgr)
Date: Fri, 27 Mar 1998 19:24:21 +0100 (MET)
X-Mailer: ELM [version 2.4 PL22]
MIME-Version: 1.0
Content-Type: text/plain; charset=ISO8859-1
Content-Transfer-Encoding: 8bit
Status: O

      dimension rhos(15),es(15),eks(15),eps(15)
      dimension slop(14),rhosl(14)
c
c     direct quantities
c
      read(5,'(1x)')
      do i=1,15
        read(5,*) rho,e
        rhos(i) = log(rho)
        es(i)   = log(e)
        ekin = 0.6*20.75*(3.0*3.14596**2)**0.66667*rho**0.66667
        eks(i)  = log(ekin)
        eps(i)  = log(abs(ekin-e))
        write(6,'(1x,8(1pg12.4))') 
     &    rho,e,ekin,e-ekin,log(rho),log(e),log(ekin),log(abs(e-ekin))
      enddo
c
c     slopes
c
      rold = rhos(1)
      eold = es(1)
      ekold = eks(1)
      epold = eps(1)
      do i=2,15
        fac = 1.0/(rhos(i)-rold)
        rhoact = exp(0.5*(rhos(i)+rold))
        rhosl(i-1) = log(rhoact)
        slope = (eps(i)-epold)*fac
        slop(i-1)  = log(abs(slope-1.0))
        write(6,'(1x,4(1pg12.4))')
     &  rhoact,(es(i)-eold)*fac,(eks(i)-ekold)*fac,slope
        eold  = es(i)
        ekold = eks(i)
        epold = eps(i)
        rold  = rhos(i)
      enddo
      
c
c     slope of slope
c
      rold = rhosl(1)
      sold = slop(1)
      do i=2,14
        fac = 1.0/(rhosl(i)-rold)
        rhoact = exp(0.5*(rhosl(i)+rold))
        write(6,'(1x,2(1pg12.4))') rhoact,(slop(i)-sold)*fac
        rold = rhosl(i)
        sold = slop(i)
      enddo
c
c     calibrated potential energy
c
      write(6,'(/a)') ' calibrated slopes'
      do i=1,15
        rho  = exp(rhos(i))
        epc  = exp(eps(i))/rho
        write(6,'(1x,4(1pg12.4))')
     &  rho,epc,log(rho),log(abs(epc))
      enddo
c
c     calibrated slope of potential energy
c
      rold = exp(rhos(1))
      epcold = exp(eps(1))/rold
      write(6,'(/a)') ' calibrated slopes'
      do i=2,15
        rho  = exp(rhos(i))
        epc  = exp(eps(i))/rho
        slope = (epc-epcold)/(rho-rold)
        rhoact = 0.5*(rho+rold)
        rhosl(i-1) = log(rhoact)
        slop(i-1)  = log(abs(slope))
        write(6,'(1x,4(1pg12.4))')
     &  rhoact,slope,log(rhoact),log(abs(slope))
        epcold = epc
        rold  = rho
      enddo
      
c
c     slope of slope
c
      rold = rhosl(1)
      sold = slop(1)
      do i=2,14
        fac = 1.0/(rhos(i)-rold)
        rhoact = exp(0.5*(rhosl(i)+rold))
        write(6,'(1x,2(1pg12.4))') rhoact,(slop(i)-sold)*fac
        rold = rhosl(i)
        sold = slop(i)
      enddo
      stop
      end


