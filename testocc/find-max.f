! finds maximum of strength (second column) in a file
program find_max
!
  implicit none
!
  integer i
  real energ,str,enrmax,strmax
!
!--------------------------------------------------------------------
!
  enrmax = 0.0
  strmax = 0.0
  do i=1,9
    read(5,'(x)')
  enddo
  do i=1,9999
    read(5,*,err=9,end=9) energ,str
    if(str.gt.strmax) then
      enrmax = energ
      strmax = str
    endif
  enddo
 9 continue
  write(6,'(f7.3)') enrmax
!
  stop
end 
