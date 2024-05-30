      program testeq
c
      implicit double precision (a-h,o-z)
c
      data d2/1.0/
      data iteqdf/20/
      data epeqf/1.0D-8/
      data epeqdx/1.0D-8/
      data alpha/0.4/
c
c
      do itest=-10,10
      e = 2.0D0*itest
      eq=sqrt(e*e+d2)
      v2=0.5*(1.-e/eq)
c
c      old iteration
c
c
c         initial guess for ep
c
      x  = e
      eq=sqrt(e*e+d2)
cold      v2=0.5*(1.-e/eq)
cold      f=v2*alpha
      f=e+alpha*v2-x
      derf=-0.5*alpha*d2/eq**3-1.
c
cold      if (abs(f).lt.epeqf) go to 7
c      
cold      x=e
c
c         start iteration loop
c
      do 6 i=1,iteqdf
      itact = i
c
      dx=f/derf
      x=x-dx
c
      eq=sqrt(d2+x*x)
      v2=0.5*(1.-x/eq)      
      f=e+alpha*v2-x
      derf=-alpha*0.5*d2/eq**3-1.
c
c        test on satisfactory accuracy
c
      if(abs(dx).ge.epeqdx) go to 6
c
      if(abs(f).le.epeqf) go to 7
c
    6 continue
c
      write(6,100)
 100  format(1x,'maximal number of iterations exceeded. stop in eqdef') 
c
c     stop                                                    ! mod
c
    7 continue
      write(6,'(a,3g12.4)')
     &  ' old iteration: ITER,DX,F=',itact,dx,f
c
c     new iteration strategy #1
c
      eq=sqrt(e*e+d2)
      v2=0.5*(1.-e/eq)
        v2old = v2
        do it=1,iteqdf
        itact = it
        x  = e+alpha*v2
        eq=sqrt(d2+x*x)
        v2=0.5-0.5*x/eq
        if(abs(v2-v2old).lt.epeqf) goto 999
        v2old = v2
        enddo
c        write(6,'(a,i5)')
c     &   ' no convergence in EQDEF at K=',k
  999   continue
      write(6,'(a,3g12.4)')
     &  ' new iteration: ITER,DX=',itact,v2-v2old
c
c
c     new iteration strategy #2
c
      eq=sqrt(e*e+d2)
      v2=0.5*(1.-e/eq)
        v2old = v2
        do it=1,iteqdf
        itact = it
        x  = e+alpha*v2
        eq=sqrt(d2+x*x)
        v2=0.5-0.5*x/eq
        derv2 = alpha*d2*0.5/(eq*eq*eq)
        if(abs(v2-v2old).lt.epeqf) goto 998
        delv  = (v2old-v2)/derv2
        v2old = v2
        v2    = v2+delv
        enddo
c        write(6,'(a,i5)')
c     &   ' no convergence in EQDEF at K=',k
  998   continue
      write(6,'(a,3g12.4)')
     &  ' 3. iteration: ITER,DX=',itact,v2-v2old
c
      enddo
c
      stop
      end
