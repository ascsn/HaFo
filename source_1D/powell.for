c-----powell---------------------------------------------------------
c
      subroutine powell(a,nterms,nwdr)
c
c     chi**2 minimization by powell's comjugate directions method
c
c     dimension statement valid for nterms up to 'kparmx'
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
c
      external ffverg
c
      dimension a(kparmx)
c
c     workspaces
c
      dimension direct(kparmx,kparmx)
      dimension parmom(kparmx)
      equivalence (direct,fderij),(chimom,parmom)
c
!      data  eps/0.001/
c
c----------------------------------------------------------------------
c
      eps=endcon
      WRITE(2,*) 'POWELL with termination criterion eps=',eps
c
c     initialize starting vector and matrix of directions
c
      do i=1,kparmx
       parmom(i)   = a(i)
      enddo
      do i=1,kparmx
      do j=1,kparmx
       if(i.eq.j) then
         direct(i,j) = dparmo(i)
       else
         direct(i,j) = 0.0
       endif
      enddo
      enddo
c
c     powell's iteration
c
      call powels(parmom,direct,nterms,kparmx,eps,iter,nwdr,chiobt)
      do i=1,kparmx
        a(i) = parmom(i)
      enddo
c
      return
c
      end
c-----powels-----------------------------------------------------------
c
c     package from numercical recipes, slightly modified
c
      SUBROUTINE powels(p,xi,n,np,ftol,iter,itmax,fret)
      INTEGER iter,n,np,NMAX,ITMAX
c      DOUBLE PRECISION fret,ftol,p(np),xi(np,np),func
c      EXTERNAL func
      DOUBLE PRECISION fret,ftol,p(np),xi(np,np),ffverg
      EXTERNAL ffverg
      PARAMETER (NMAX=40)
c     PARAMETER (NMAX=20,ITMAX=10)
CU    USES func,linmin
      INTEGER i,ibig,j
      DOUBLE PRECISION del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      double precision pcom(NMAX),xicom(NMAX)
      COMMON /f1com/ pcom,xicom,chiminsav,ncom
!      data chiminsav/0D0/
c
      do i=1,np
        pcom(i) = p(i)
        pt(i)   = p(i)
        ptt(i)  = p(i)
      enddo
c      fret=func(p)
      fret=ffverg(p)
      chiminsav = fret
        write(2,'(a,i5,1pg12.4,5(/1x,5g14.6))')
     &   ' 0. iter,chi**2=',iter,fret,
     &   (p(i),i=1,n)
        write(6,'(a,i5,1pg12.4)')
     &   ' 0. iter,chi**2=',iter,fret
      CALL flush(2)
      do 11 j=1,np
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.d0
      do 13 i=1,n
        do 12 j=1,np
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        CAll setwrite(-2)
!        WRITE(2,'(a,i5)') 'CALL n=',n
        call linmin(p,xit,n,fret)
        write(2,'(a,2i5,1pg14.6,5(/1x,5g14.6))')
     &   ' I. iter,ipar,chi**2=',iter,i,fret,
     &   (p(i2),i2=1,n)
        CALL flush(2)
        write(6,'(a,i5,1pg14.6)')
     &   ' I. iter,chi**2=',iter,fret
        IF(fret<chiminsav) THEN
          chiminsav = fret
          chmin=fret
          CALL write_params(p,chiminsav,iter)
        END IF 
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.d0*abs(fp-fret).le.ftol*(abs(fp)+abs(fret))) return
      if(iter.eq.ITMAX) return 
      do 14 j=1,np
        ptt(j)=2.d0*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
c      fptt=func(ptt)
!      CAll setwrite(-1)
      fptt=ffverg(ptt)
      IF(fptt<chiminsav) THEN
        chiminsav = fptt
        chimin = fptt
        CALL write_params(ptt,chiminsav,iter)
      END IF 
      CAll setwrite(-2)
      if(fptt.ge.fp)goto 1
      t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.d0)goto 1
      call linmin(p,xit,n,fret)
        write(2,'(a,i5,1pg14.6,5(/1x,5g14.6))')
     &   ' II. iter,chi**2=',iter,fret,
     &   (p(i),i=1,n)
        CALL flush(2)
        write(6,'(a,i5,1pg14.6)')
     &   ' II. iter,chi**2=',iter,fret
        IF(fret<chiminsav) THEN
          chiminsav = fret
          chimin = fret
          CALL write_params(p,chiminsav,iter)
        END IF 
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END
      SUBROUTINE linmin(p,xi,n,fret)
      INCLUDE 'comfit.inc'
      INTEGER n,NMAX
      DOUBLE PRECISION fret,p(*),xi(*),TOL
      PARAMETER (NMAX=40,TOL=1.d-3)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
     *,brent
      COMMON /f1com/ pcom,xicom,chiminsav,ncom
      EXTERNAL f1dim
!      WRITE(2,'(a,2i5)') 'CALLED n=',n,kparmx
      ncom=n
      do 11 j=1,kparmx
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
      SUBROUTINE setwrite(level)
      include 'numpar.inc'
      iwr = level      
      END
      FUNCTION f1dim(x)
      INCLUDE 'comfit.inc'
      INTEGER NMAX
c      DOUBLE PRECISION f1dim,func,x
      DOUBLE PRECISION f1dim,ffverg,x
      PARAMETER (NMAX=40)
!      parameter (kparmx=17)
CU    USES func
      INTEGER j,ncom
      DOUBLE PRECISION pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,chiminsav,ncom
      do j=1,kparmx
        xt(j) = pcom(j)
      enddo
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
c      write(6,'(a,20(1pg14.6))')
c     &  ' F1DIM: params=',(xt(j),j=1,ncom)
c      write(2,'(a,20(1pg14.6))')
c     &  ' F1DIM: params=',(xt(j),j=1,ncom)
c      f1dim=func(xt)
      chiact = ffverg(xt)
      IF(chiact<chiminsav) THEN
        chiminsav = chiact
        CALL write_params(xt,chiminsav,iter)
      END IF 
c      write(2,'(a,1pg12.4)') ' chi**2=',chiact
c      write(6,'(a,1pg12.4)') ' chi**2=',chiact
      f1dim=chiact
      return
      END
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      DOUBLE PRECISION brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=4,CGOLD=.3819660d0,ZEPS=1.0d-10)
      INTEGER iter
      DOUBLE PRECISION a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x
     *,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)
     *) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        write(2,'(a,i5,g13.5)') ' BRENT: iter,chi=',iter,fu
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
c      pause 'brent exceed maximum iterations'
      write(6,'(a)') 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      write(2,'(a,i5)') ' BRENT: nr. iterations=',iter
      write(6,'(a,i5)') ' BRENT: nr. iterations=',iter
      return
      END
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      DOUBLE PRECISION ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
      DOUBLE PRECISION dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
      WRITE(2,'(a,3(1pg13.5))') 'MNBRAK: fa,fb,fc=',fa,fb,fc
      icount=0
1     if(fb.ge.fc)then
        icount=1+icount
        IF(icount>20) THEN
          WRITE(2,'(a)') 'MNBRAK: repetition count exhausted --> return'
          STOP "repetition count exhausted"
       END IF
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        WRITE(2,'(a,2(1pg13.5))') 'MNBRAK: fu=',fu
        goto 1
      endif
      return
      END
      
