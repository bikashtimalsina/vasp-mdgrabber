module handleoutcar
implicit none
! This is the possible minimal structure
type structure
integer(8) nAtom
real(8) energy, volume
real(8), dimension(3) :: posAtom, forceAtom
real(8), dimension(6) :: virial
real(8),dimension(3,3) :: box
end type structure
! The number of configuration or snapshots if present in OUTCAR
type, extends(structure) :: snapshots
integer(4) nconfig
end type snapshots
contains
! method to match a word from the text file
subroutine findword(word,line,found,i)
logical found
character(*) line, word
integer(4) i,l,k
l = len_trim(line)
k=len_trim(word)
found = .False.
do i=1,l-k+1
  if(line(i:i+k-1) .eq. word(1:k)) then
      found = .True.
      exit
  endif
enddo
if (line .eq. word) then
  found=.True.
endif
end subroutine findword
! routine to find the number of atom on the OUTCAR
end module handleoutcar

program test
use handleoutcar
character line*199, word*59, newformat*8
type(snapshots) snap
integer(8) newindex,ntot
integer(4) i,j,outread,maxl,reason,k,p,num1,num2,ncount,outfile,r
real(8),allocatable :: energy(:), positionFirst(:,:),forceFirst(:,:)
logical found
maxl=99999999
outread=432
open(outread,file="./OUTCAR")
j=0
outfile=123
!-------------------------------------------------------------------------------
! This is the loop to find number of atom on OUTCAR file
ncount=0
do i=1,maxl
word='free  energy   TOTEN'
call findword(word,line,found,p)
if (found) then
ncount=ncount+1
endif
word='NIONS'
read(outread,'(A)',iostat=reason) line
call findword(word,line,found,p)
if (found) then
j=i
read(line(p+8:),'(I10)') snap%nAtom
write(*,*) "Number of ions in the structure: ", snap%nAtom
!exit
endif
if (reason < 0) then
exit
endif
enddo
close(outread)
snap%nconfig=ncount
allocate(energy(snap%nconfig))
ntot=ncount*snap%nAtom
allocate(positionFirst(snap%nAtom,3),forceFirst(snap%nAtom,3))
! ------------------------------------------------------------------------------
! This is the loop to find number of snapshots from an OUTCAR FILE
k=1
r=1
open(outread,file="./OUTCAR")
open(outfile,file="configuration.txt")
do i=1,maxl
word='free  energy   TOTEN'
read(outread,'(A)',iostat=reason) line
call findword(word,line,found,p)
if (found) then
j=i
num1=len_trim(line(len_trim(word)+11:len_trim(line)-3))
num1=num1-1
num2=len_trim(line(len_trim(word)+11:len_trim(line)-3))-9
num2=num2-1+num1
write(newformat,'(A,I2.2,A,I2.2,A)') '(F',num2,'.',num1,')'
read(line(len_trim(word)+11:len_trim(line)-3),newformat) snap%energy
energy(k)=snap%energy
k=k+1
endif
! Below here is to obtain the position and force on each atom from the OUTCAR file
word="POSITION"
call findword(word,line,found,p)
if (found) then
j=i
read(outread,'(A)',iostat=reason) line
if (r .gt. 1) then
write(outfile,*) snap%energy
endif
do newindex=1,snap%nAtom
read(outread,*,iostat=reason) snap%posAtom(1), snap%posAtom(2), snap%posAtom(3), &
snap%forceAtom(1), snap%forceAtom(2), snap%forceAtom(3)
if (r .eq. 1) then
positionFirst(newindex,1)=snap%posAtom(1)
positionFirst(newindex,2)=snap%posAtom(2)
positionFirst(newindex,3)=snap%posAtom(3)
forceFirst(newindex,1)=snap%forceAtom(1)
forceFirst(newindex,2)=snap%forceAtom(2)
forceFirst(newindex,3)=snap%forceAtom(3)
endif
if(r .gt. 1) then
write(outfile,*) snap%posAtom(1), snap%posAtom(2), snap%posAtom(3), &
snap%forceAtom(1), snap%forceAtom(2), snap%forceAtom(3)
endif
enddo
if (r .eq. 1) then
do newindex=1,snap%nAtom
write(outfile,*) positionFirst(newindex,1),positionFirst(newindex,2), &
positionFirst(newindex,3),forceFirst(newindex,1),forceFirst(newindex,2), &
forceFirst(newindex,3)
enddo
endif
r=r+1
endif
if (reason < 0) then
write(*,*) "END OF FILE"
exit
endif
enddo
close(outread)
close(outfile)
write(*,*) "The number of configuration is: ", snap%nconfig
! ------------------------------------------------------------------------------
deallocate(energy,positionFirst,forceFirst)
end program test
!free  energy   TOTEN  =     -1165.77741107 eV  (F17.13)
