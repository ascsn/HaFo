c-----centroids-----------------------------------------------------
c
      subroutine centroid_window(omega,maxout)
c
      include 'params.inc'
c
c     determins appropriate windows for computing energy centroids
c
      include 'constn.inc'
      include 'wavfun.inc'
      include 'linres.inc'
      dimension omega(krpmax)
      REAL(8),PARAMETER :: epsconv=1D-3
      INTEGER,PARAMETER :: ite1max=50
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
!      WRITE(14,*) 'centroid_window: lrpa,nprot=',lrpa,nprot
      IF(nprot==8) THEN
        IF(lrpa==0) THEN
          energmin = 11D0
          energmax = 40D0
        ELSE IF(lrpa==1) THEN
          energmin = 15D0
          energmax = 30D0
        ElSE IF(lrpa==2) THEN
          energmin = 11D0
          energmax = 40D0
        ELSE
          return 
        END IF
      ELSE IF(nprot==82) THEN
        IF(lrpa==0) THEN
          energmin = 10D0
          energmax = 35D0
        ELSE IF(lrpa==1) THEN
          energmin = 8D0
          energmax = 25D0
        ElSE IF(lrpa==2) THEN
          energmin = 8.6D0 !8D0
          energmax = 16D0
          e1scold = 12D0
          DO iter=1,ite1max
           CALL centroids(omega,maxout,energmin,energmax,e1pr,e1sc,e1vc)
           IF(ABS(e1scold-e1sc)<epsconv) EXIT
           WRITE(6,'(a,i3,1pg13.5)') ' iter,e1=',iter,e1sc
           e1scold = e1sc
           energmax = 2*e1sc-energmin
          END DO
        ELSE
          return 
        END IF
      ELSE
        return 
      END IF
      CALL centroids(omega,maxout,energmin,energmax,e1pr,e1sc,e1vc)
      write(14,'(2(a,5x,3(9x,g13.5)/))')
     &   ' e1-ct:  ',e1pr,e1sc,e1vc
c      
      RETURN
      END

c-----centroids-----------------------------------------------------
c
      subroutine centroids(omega,maxout,energmin,energmax,
     &                     e1pr,e1sc,e1vc)
c
      include 'params.inc'
c
c     computes energy centroids in window [energmin:energmax]
c
      include 'constn.inc'
      dimension omega(krpmax)
      common /bels/belsavpr(krpmax),belsavne(krpmax),
     &             belsavsc(krpmax),belsavvc(krpmax)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
!      WRITE(14,*) 'centroids: energmin,energmax=',energmin,energmax
      ewsrpr = 0D0
      ewsrsc = 0D0
      ewsrvc = 0D0
      sumblp = 0D0
      sumbls = 0D0
      sumblv = 0D0
      do imode=1,maxout
        energy = omega(imode)
        IF(energmin<energy .AND. energy<energmax) THEN
          bl2pr  = belsavpr(imode)
          bl2sc  = belsavsc(imode)
          bl2vc  = belsavvc(imode)
          ewpr   = bl2pr*energy
          ewsc   = bl2sc*energy
          ewvc   = bl2vc*energy
          ewsrpr = ewpr+ewsrpr
          ewsrsc = ewsc+ewsrsc
          ewsrvc = ewvc+ewsrvc
          sumblp = bl2pr+sumblp
          sumbls = bl2sc+sumbls
          sumblv = bl2vc+sumblv
        END IF
      enddo
      e1pr = ewsrpr/sumblp
      e1sc = ewsrsc/sumbls
      e1vc = ewsrvc/sumblv
c
      return
      end
