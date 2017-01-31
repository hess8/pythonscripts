MODULE vector_matrix_utilities
use num_types
use numerical_utilities
implicit none
private
public matrix_inverse, determinant, cross_product, volume, norm, reduce_C_in_ABC, &
       orthogonality_defect,&
       test_two, test_one, test_three, test_four, minkowski_reduce_basis,  &
       writefloat3x3, writeint3x3, writefloat3x3N, writeint3x3N, number_nonzero, svmesh

INTERFACE determinant
   module procedure determinant_real, determinant_integer
END INTERFACE

INTERFACE norm
   module procedure norm_real_vector, norms_real_vector_list
END INTERFACE norm
CONTAINS

!***************************************************************************************************
! This function checks the minkowski conditions for a 3D lattice basis. 
function minkowski_conditions_check(basis,eps)
logical  :: minkowski_conditions_check
real(dp), dimension(3,3), intent(IN) :: basis
real(dp), intent(IN) :: eps
real(dp), dimension(3) :: b1, b2, b3
b1 = basis(:,1)
b2 = basis(:,2)
b3 = basis(:,3)

minkowski_conditions_check = .true.
if (norm(b1) > norm(b2)+eps) then
   minkowski_conditions_check = .false.
   write(*,*) "Minkowski_condition 1 failed: b1 > b2"
endif
if (norm(b2) > norm(b3)+eps)        then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 2 failed: b2 > b3"
endif
if (norm(b2) > norm(b1+b2)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 3 failed: b2 > b1+b2"
endif
if (norm(b2) > norm(b1-b2)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 4 failed: b2 > b1-b2"
endif
if (norm(b3) > norm(b1+b3)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 5 failed: b3 > b1+b3"
endif
if (norm(b3) > norm(b3-b1)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 6 failed: b3 > b3-b1"
endif
if (norm(b3) > norm(b2+b3)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 7 failed: b3 > b2+b3"
endif
if (norm(b3) > norm(b3-b2)+eps)     then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 8 failed: b3 > b3-b2"
endif
if (norm(b3) > norm(b1+b2+b3)+eps)  then
   minkowski_conditions_check = .true.
!   write(*,*) 'b1,b2,b3',b1,b2,b3
!   write(*,*) 'norm(b3)',norm(b3)
!   write(*,*) 'norm(b1+b2+b3)',norm(b1+b2+b3)
   write(*,*) "Minkowski_condition 9 failed: b3 > b1+b2+b3"
endif
if (norm(b3) > norm(b1-b2+b3)+eps)  then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 10 failed: b3 > b1-b2+b3"
endif
if (norm(b3) > norm(b1+b2-b3)+eps)  then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 11 failed: b3 > b1+b2-b3"
endif
if (norm(b3) > norm(b1-b2-b3)+eps)  then
   minkowski_conditions_check = .true.
   write(*,*) "Minkowski_condition 12 failed: b3 > b1-b2-b3"
endif
end function minkowski_conditions_check

!***************************************************************************************************
!  This routine takes three vectors, A,B,C, defining a lattice, and reduces the last one so that
!  it is as close as possible to the origin while remaining in an affine plane, which is parallel
!  to the A-B plane but which passes through the end of the C vector. See Lecture notes in computer
!  science, ISSN 0302-974, ANTS - VI : algorithmic number theory, 2004, vol. 3076, pp. 338-357
!  ISBN 3-540-22156-5
subroutine reduce_C_in_ABC(A,B,C,eps)
real(dp), dimension(3), intent(inout) :: A, B, C
real(dp), intent(in) :: eps
real(dp), dimension(3) :: T  ! projection of C into the A-B plane
real(dp), dimension(3,3) :: ABC, ABCinv, oldABC ! Matrices of ABC basis vectors and inverse
real(dp), dimension(3)   :: cpdAB ! unit vector perpendicular to the A-B plane
real(dp) :: dist(4) ! the distances from T to enclosing lattice points of A,B (4 corners of the ppiped)
integer LC(3) ! lattice coordinates of C, in the affine plane, using the A,B basis vectors
integer idx(1) ! index of the smallest distance from T to a lattice point in A,B
logical err

ABC = reshape((/A,B,C/),(/3,3/))
oldABC = ABC

! Use Gaussian reduction to reduce the A,B 2D basis so that it is itself Minkowski reduced. If this
! is done, then the closest lattice point (in A,B plane) to the projection of C (into the A,B plane) is
! guaranteed to be one of the corners of the unit cell enclosing the projection of C
call gaussian_reduce_two_vectors(A,B,eps)

! First thing to do is find the (real, not lattice) point in the affine plane A,B + C that is
! nearest the origin. Call this T.
cpdAB = cross_product(A,B)/norm(cross_product(A,B))
T = C - cpdAB*dot_product(C,cpdAB)
if(.not. equal(dot_product(T,cross_product(A,B)),0._dp,eps)) then
   print *,dot_product(T,cross_product(A,B))
   stop "Projection of C into A,B plane failed"
endif

! Now find the four points of the A,B lattice, in the affine plane, that enclose the point T
ABC = reshape((/A,B,C/),(/3,3/))
call matrix_inverse(ABC,ABCinv,err)
if(err)stop "A,B,C vectors in reduce_C_in_ABC are co-planar"
LC = floor(matmul(ABCinv,T) + eps)
! Compute the distance from T to each of the four points and pick the one that is the closest.
dist(1) = norm(T-matmul(ABC,LC)) 
dist(2) = norm(T-matmul(ABC,(/LC(1)+1,LC(2),LC(3)/)))
dist(3) = norm(T-matmul(ABC,(/LC(1),LC(2)+1,LC(3)/)))
dist(4) = norm(T-matmul(ABC,(/LC(1)+1,LC(2)+1,LC(3)/)))

idx = minloc(dist)
select case(idx(1))
case(1)
   C = C - matmul(ABC,LC)
case(2)
   C = C - matmul(ABC,(/LC(1)+1,LC(2),LC(3)/))
case(3)
   C = C - matmul(ABC,(/LC(1),LC(2)+1,LC(3)/))
case(4)
   C = C - matmul(ABC,(/LC(1)+1,LC(2)+1,LC(3)/))
case default
   print *, "Case failed in reduce_C_in_ABC"
   write(*,'("Lattice coordinates in the A,B plane: ",2(i2,1x))') LC
   stop
end select
ABC = reshape((/A,B,C/),(/3,3/))
call matrix_inverse(ABC,ABCinv,err)
if(any(abs(matmul(ABCinv,oldABC)-nint(matmul(ABCinv,oldABC)))>eps)) stop "Lattice was not preserved &
& in reduce_A_in_ABC"

if(err)stop "A,B,C vectors in reduce_C_in_ABC are co-planar (after Minkowski)"

endsubroutine reduce_C_in_ABC

!***************************************************************************************************
! This routine takes two vectors (in three-space) and reduces them to form a shortest set (Minkowski
! reduced). The idea is to subtract multiples of U from V so that the new V is as close to the
! origin as any lattice point along the line that passes through U in the direction of V. The
! process is repeated until the new vector isn't shorter than the other. It's pretty obvious if you
! do an example by hand. Also see 3.1 of Lecture notes in computer science, ISSN 0302-974, ANTS - VI
! : algorithmic number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
!
! Fixes made Apr 2012 GLWH (not sure if they made a practical difference though)
!
subroutine gaussian_reduce_two_vectors(U,V,eps)
real(dp), dimension(3) :: U, V, R
real(dp), intent(in) :: eps

real(dp) temp(3)
integer it

it = 0
if (norm(U) > norm(V) - eps) then 
   ! Make sure that the {U,V} are listed in ascending order; ||U||<||V||
   temp = U; U = V; V = temp ! Keep V as the longest vector
end if
do; it = it + 1
   if (it > 10) stop "gaussian_reduce_two_vectors failed to converge in 10 iterations"
   R = V-nint(dot_product(U,V)/dot_product(U,U))*U !Shorten V as much as possible
   V = U ! Swap U and V (so U remains the shortest)
   U = R
   if (norm(U) >= norm(V) - eps) exit
end do
! Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
temp = U; U = V; V = temp

endsubroutine gaussian_reduce_two_vectors

!***************************************************************************************************
SUBROUTINE minkowski_reduce_basis(IN,OUT,eps)
real(dp), dimension(3,3), intent(in) :: IN !bch The three vectors are the columns of IN.
real(dp), dimension(3,3), intent(out) :: OUT
real(dp), intent(in) :: eps
real(dp)             :: norms(3), temp(3,3)
integer              :: i, idx(1), it, limit

limit = 10

if (equal(determinant(IN),0._dp,eps)) stop "Input basis for 'minkowski_reduce_basis' was not linearly independent"
OUT = IN

! Keep applying the greedy algorithm until the vectors come out already sorted
do it = 1, limit
   ! Sort the three vectors into ascending order
   temp = OUT
   norms = norm(temp)
   do i = 3,1,-1
      idx = maxloc(norms)
      temp(:,i) = OUT(:,idx(1))
      norms(idx(1)) = 0
   enddo
   OUT = temp ! Copy the sorted vectors back to OUT
   call reduce_C_in_ABC(OUT(:,1),OUT(:,2),OUT(:,3),eps)
   if (norm(OUT(:,3))>=norm(OUT(:,2))-eps) exit
   
end do
if (it==limit+1) stop "Too many iterations in 'minkowski_reduce_basis'"
if (.not. minkowski_conditions_check(OUT,eps)) stop "ERROR in minkowski_reduce_basis: Minkowski conditions not met."

! we want to make sure that the det is positive.
! NOTE: This *destroys* the mathematical picture of a "greedy qreduced basis" (Minkowski), but
!       from a physical point of view we don't care ;-)
!       Either way, the basis is as orthogonal as possible.
if (determinant(OUT)<0) then 
  temp(:,1) = OUT(:,2)
  OUT(:,2) = OUT(:,3)
  OUT(:,3) = temp(:,1)
endif
!call writefloat3x3(OUT,'matrix OUT, mink routine') !bch

END SUBROUTINE minkowski_reduce_basis

!***************************************************************************************************
! This function calculates the "orthogonality defect" of the given basis of a 3D lattice.
function orthogonality_defect(basis)
real(dp), dimension(3,3), intent(in) :: basis
real(dp) :: orthogonality_defect
real(dp) :: od
integer j
od = 1._dp
do j=1,3
   od = od*norm(basis(:,j))
enddo
od = od/abs(determinant(basis))
orthogonality_defect = od
endfunction

!!***************************************************************************************************
!! This routine takes a set of basis vectors (that form a lattice) and reduces them so that they form
!! the shortest possible basis. See Lecture Notes in Computer Science, ISSN 0302-974, ANTS - VI :
!! algorithmic number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
!SUBROUTINE reduce_to_shortest_basis(IN,OUT,eps)
!real(dp), dimension(3,3), intent(in) :: IN
!real(dp), dimension(3,3), intent(out) :: OUT
!real(dp), intent(in) :: eps
!
!real(dp), dimension(3) :: A,B,C ! the basis vectors
!real(dp) :: check(3,3), temp(3), old(3,3)
!logical err
!real(dp) :: od, odnew
!integer it!, j
!it = 0
!A = IN(:,1); B = IN(:,2); C = IN(:,3)
!odnew = orthogonality_defect(IN)
!!open(15,file="debug_reduction.out")
!!write(15,'("Before reduction, the orthogonality defect of the basis was",f7.3)') odnew 
!!write(*,'(f7.3)') odnew 
!!write(*,'(i5,1x,3("out: ",3(f11.5,1x)))') it,A,B,C
!
!
!do ; it = it + 1
!!   write(*, '("iter: ",i3)') it
!   od = odnew
!   old = reshape((/A,B,C/),(/3,3/))
!   call reduce_A_in_ABC(A,B,C,eps)
!   call reduce_A_in_ABC(B,C,A,eps)
!   call reduce_A_in_ABC(C,A,B,eps)
!!write(*,'(i5,1x,3("out: ",3(f11.5,1x)))') it,A,B,C
!   OUT = reshape((/A,B,C/),(/3,3/))
!   !print *,"outer loop"
!!do j = 1,3
!!   write(*,'(3(f11.5,1x))') OUT(j,:)
!!enddo
!!print *
!   odnew = orthogonality_defect(OUT)
!!   write(*,'("it: ",i3,1x,"OD: ",2(f7.3,1x))') it, odnew, od
!   if (equal(od,odnew,eps)) exit
!   if(odnew > od) then ! this check is just a temporary fix. Sometimes (very rarely!) the reduction
!      ! makes the new basis *less* orthogonal. Not sure why this is...
!      print *,"Orthogonality defect in reduce_to_shortest_basis *increased*. Not deadly but fix it someday"
!      write(*,'("final OD was ",f7.3)') od
!      A = old(:,1); B = old(:,2); C = old(:,3)
!      exit
!   endif
!   if (equal(old,OUT,eps)) exit
!   if (it>10) stop "More than 10 iterations in outer loop of reduce_to_shortest_basis"
!
!!!! debug
!!   call matrix_inverse(OUT,check,err)
!!   if(.not. equal(matmul(check,IN)-nint(matmul(check,IN)),0._dp,eps)) stop "ERROR: changed the basis"
!!!! end debug
!
!enddo
!call matrix_inverse(OUT,check,err)
!if(err) stop "OUT matrix is singular in reduce_to_shortest_basis"
!! Check that the conversion from old to new lattice vectors is still an integer matrix
!if(.not. equal(matmul(check,IN)-nint(matmul(check,IN)),0._dp,eps)) stop "ERROR: Reduced lattice vectors &
!& in reduce_to_shortest_basis changed the original lattice"
!if (determinant(OUT) < 0._dp) then ! Left-handed coordinate system so switch to right-handed
!   temp = OUT(:,3); OUT(:,3) = OUT(:,2); OUT(:,2) = temp ! exchange vectors 2 and 3
!endif
!!write(15,'(" After reduction, the orthogonality defect of the basis is ",3f7.3)') orthogonality_defect(OUT)
!!close(15)
!!write(*,'(f7.3)') orthogonality_defect(OUT)
!ENDSUBROUTINE reduce_to_shortest_basis



!****************************************************************************************
! Given the 3x3 matrix a, finds its inverse b
subroutine matrix_inverse(a, b, err_, eps_)
real(dp), intent(in):: a(3,3)
real(dp),intent(out):: b(3,3)
real(dp) c,avec(9)
real(dp)           :: eps

logical,  optional :: err_
real(dp), optional :: eps_


if(present(err_)) err_ = .false.
if(present(eps_)) then; eps=eps_; else; eps=10d-14; endif

   c=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
      a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
      a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
   avec=(/-a(2,3)*a(3,2)+a(2,2)*a(3,3),a(1,3)*a(3,2)-a(1,2)*a(3,3), &
          -a(1,3)*a(2,2)+a(1,2)*a(2,3),a(2,3)*a(3,1)-a(2,1)*a(3,3), &
          -a(1,3)*a(3,1)+a(1,1)*a(3,3),a(1,3)*a(2,1)-a(1,1)*a(2,3), &
          -a(2,2)*a(3,1)+a(2,1)*a(3,2),a(1,2)*a(3,1)-a(1,1)*a(3,2), &
          -a(1,2)*a(2,1)+a(1,1)*a(2,2)/)
   if(abs(c) < eps) then; if(present(err_)) err_ = .true.
   else; b=1/c*transpose(reshape(avec,(/3,3/))); endif ! transpose due to column major storage
end subroutine matrix_inverse

!****************************************************************************************
! Given vectors a and b, c = a x b
function cross_product(a,b)
real(dp) a(3), b(3), cross_product(3)

cross_product = (/a(2)*b(3) - a(3)*b(2), &
                  a(3)*b(1) - a(1)*b(3), &
                  a(1)*b(2) - a(2)*b(1)/)
end function cross_product

!****************************************************************************************
! This routine takes the norm of a vector
pure function norm_real_vector(vector)
real(dp) norm_real_vector
real(dp), intent(in):: vector(:)

norm_real_vector = sqrt(dot_product(vector,vector))
end function norm_real_vector

! This routine takes the norms of vectors in a list
pure function norms_real_vector_list(vector_list)
real(dp), intent(in)               :: vector_list(:,:)
real(dp) :: norms_real_vector_list(size(vector_list,2))
integer i

do i = 1, size(vector_list,2)
   norms_real_vector_list(i) = sqrt(dot_product(vector_list(:,i),vector_list(:,i)))
end do
end function norms_real_vector_list

!****************************************************************************************
! This routine finds the determinant of 2x2 or 3x3 matrices
function determinant_real(a)
real(dp) determinant_real,a(:,:)

integer n
n=size(a,dim=1)
if (n==2) then
   determinant_real=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
else if (n==3) then
   determinant_real=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
                    a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
                    a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
else
   stop "DETERMINANT ERROR: Matrix dimension exceeds 3."
end if
end function determinant_real

function determinant_integer(a)
integer determinant_integer,a(:,:)
integer n
n=size(a,dim=1)
if (n==2) then
   determinant_integer=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
else if (n==3) then
   determinant_integer=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
                       a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
                       a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
else
   stop "DETERMINANT ERROR: Matrix dimension exceeds 3."
end if
end function determinant_integer


!****************************************************************************************
! This function takes three vectors and returns a "signed" volume of the parallelpiped
! that they form
function volume(a1, a2, a3)
real(dp) volume
real(dp) a1(3), a2(3), a3(3)

volume = dot_product(a1, cross_product(a2,a3))
end function volume


!***************************************************************************************************
! This function takes a number, prints it, returns
function test_one(a)
integer test_one, a
print *,"Hello from test_one "
write(*,'(2(i3,1x))') a
test_one = a
print *,"Goodbye from test_one"
end function test_one


!***************************************************************************************************
! This function calculates the product of two numbers
function test_two(a,b)
integer test_two, a, b
print *,"Hello from test_two_multiply"
write(*,'(2(i3,1x))') a,b
test_two = a*b
print *,"Goodbye from test_two_multiply"
end function test_two

!***************************************************************************************************                                                           
! This function calculates the product of two numbers                                                                                                          
function test_three(a,b)
real(dp) test_three, a, b
print *,"Hello from test_three"
write(*,'(2(f7.4,1x))') a,b
test_three = a*b
print *,"Goodbye from test_three"
end function test_three

!***************************************************************************************************
! This function calculates the product of two arrays                                                                                                         
function test_four(dim,a)
integer dim
real(dp) a(dim)
real(dp) test_four
print *,"Hello from test_four"
print *, "dim",dim
print *,"first element",a(1)
print *,"size",size(a)
print *,shape(a)
write(*,'(20(f7.4,1x))') a
test_four = sum(a)
print *,"Goodbye from test_four"
end function test_four

!***************************************************************************************************                                                           
! This function calculates the product of two arrays                                                                                                           
function test_five(r,c,a)
integer r,c
real(dp) a(r,c)
real(dp) test_five
print *,"Hello from test_five"
write(*,'(2(i3,1x))') r,c
test_five = sum(a)
print *,"Goodbye from test_five"
end function test_five


!BCH  Write out integer 3x3xN matrix*******************
    Subroutine writeint3x3N(mat,N,label)
    use num_types
    IMPLICIT REAL (A-H,O-Z)
    integer irow, jmat
    integer, intent(in):: N
    integer, dimension(3,3,N):: mat
    character(len=*) :: label
    WRITE(*,*) ' '
    WRITE(*,*) trim(label)
    WRITE(*,*) ' '
    do jmat = 1,N
    	if (jmat < 10) WRITE(*,'(A,i1,A)') 'mat',jmat,'='
    	if (jmat >= 10) WRITE(*,'(A,i2,A)') 'mat',jmat,'='
      	WRITE(*,400) ' {{',mat(1,1,jmat),',',mat(1,2,jmat),',',mat(1,3,jmat),'},'
      	WRITE(*,400) '  {',mat(2,1,jmat),',',mat(2,2,jmat),',',mat(2,3,jmat),'},'
      	WRITE(*,400) '  {',mat(3,1,jmat),',',mat(3,2,jmat),',',mat(3,3,jmat),'}}'
      	WRITE(*,*) ' '
	end do
    400 FORMAT(A,I8,A,I8,A,I8,A)
	end subroutine writeint3x3N
! BCH end

!BCH  Write out float 3x3xN matrix*******************
    Subroutine writefloat3x3N(mat,N,label)
    use num_types
    implicit real (A-H,o-Z)
    integer irow, jmat
    integer, intent(in):: N
	real(dp), intent(in) :: mat(3,3,N)
    character(len=*) :: label
    WRITE(*,*) ' '
    WRITE(*,*) trim(label)
    WRITE(*,*) ' '
    do jmat = 1,N
    	if (jmat < 10) WRITE(*,'(A,i1,A)') 'mat',jmat,'='
    	if (jmat >= 10) WRITE(*,'(A,i2,A)') 'mat',jmat,'='
      	WRITE(*,400) ' {{',mat(1,1,jmat),',',mat(1,2,jmat),',',mat(1,3,jmat),'},'
      	WRITE(*,400) '  {',mat(2,1,jmat),',',mat(2,2,jmat),',',mat(2,3,jmat),'},'
      	WRITE(*,400) '  {',mat(3,1,jmat),',',mat(3,2,jmat),',',mat(3,3,jmat),'}}'
      	WRITE(*,*) ' '
	end do
    400 FORMAT(A,F15.10,A,F15.10,A,F15.10,A)
	end subroutine writefloat3x3N
! BCH end

!BCH  Write out integer 3x3 matrix*******************
    Subroutine writeint3x3(mat,label)
    use num_types
    IMPLICIT REAL (A-H,O-Z)
    integer irow
	integer, intent(in) :: mat(3,3)
    character(len=*) :: label
    WRITE(*,*) ' '
    WRITE(*,*) trim(label)
    WRITE(*,*) ' '
    WRITE(*,400) ' {{',mat(1,1),',',mat(1,2),',',mat(1,3),'},'
    WRITE(*,400) '  {',mat(2,1),',',mat(2,2),',',mat(2,3),'},'
    WRITE(*,400) '  {',mat(3,1),',',mat(3,2),',',mat(3,3),'}}'
    WRITE(*,*) ' '
    400 FORMAT(A,I8,A,I8,A,I8,A)
	end subroutine writeint3x3
!!! BCH end

!BCH  Write out float 3x3 matrix*******************
    Subroutine writefloat3x3(mat,label)
    use num_types
    IMPLICIT REAL (A-H,O-Z)
    integer irow
	real(dp), intent(in) :: mat(3,3)
    character(len=*) :: label
    WRITE(*,*) ' '
    WRITE(*,*) trim(label)
    WRITE(*,*) ' '
    WRITE(*,400) ' {{',mat(1,1),',',mat(1,2),',',mat(1,3),'},'
    WRITE(*,400) '  {',mat(2,1),',',mat(2,2),',',mat(2,3),'},'
    WRITE(*,400) '  {',mat(3,1),',',mat(3,2),',',mat(3,3),'}}'
	WRITE(*,*) 'norm1', norm(mat(:,1))
	WRITE(*,*) 'norm2', norm(mat(:,2))
	WRITE(*,*) 'norm3', norm(mat(:,3))
	write(*,*)
    400 FORMAT(A,F15.10,A,F15.10,A,F15.10,A)
	end subroutine writefloat3x3
! BCH end

!BCH
function number_nonzero(lattpg_op)!find number of nonzero 3x3 matrices in a 3x3x48 matrix BCH
	real(dp):: lattpg_op(3,3,48)
    integer :: number_nonzero
    integer :: irow,icol,iop
    logical :: zeros
!	real(dp), intent(in) :: lattpg_op(3,3,nop)
	do iop = 48,1,-1
	  	zeros = .true.
	  	do irow = 1,3
	  		do icol = 1,3
	  			if (abs(lattpg_op(icol,irow,iop)) .gt. 1e-10) then
	  				zeros = .false.
	  			end if
	  		end do
	  	end do
	  	if (zeros .eq. .false.) then
	  		number_nonzero = iop
	  		exit
	  	end if
	end do
end function number_nonzero

function svmesh(lat,N)
  ! Provides (the most?) compact mesh for Monkhorst Pack method (dividing lattice vectors by
  ! integers to create a mesh), by minimizing the surface/volume ratio of the mesh cell
  ! (actually only minimizes the surface as the volume of the cell is fixed.  Should work
  ! best when the lattice vectors are chosen as perpendicular as possible (Minkowski reduced).  Bret Hess
  ! Input:    N: number of mesh points desired.  lat: the lattice vectors as numpy array. The vectors are in columns
  ! Output:  m_mesh the number of divisions along each lattice vectors for the mesh

  ! The routine first calculates three real numbers (realmesh) that define the minimum surface of the cell.
  ! The real numbers aren't useful because we desire a periodic mesh commensurate with the lattice given.
  ! They are then converted in to near integers, while maintaining any multiple relations
  ! between the the real numbers
    implicit real (A-H,O-Z)
    implicit integer (i,j,m,n)
    real, dimension(3) :: r_mesh, r_pqr
    real(dp), dimension(3,3) ::lat
    integer, dimension(3) :: i_mesh, i_pqr, svmesh, pqr

    u = norm(cross_product(lat(:,1),lat(:,2)))
    v = norm(cross_product(lat(:,2),lat(:,3)))
    w = norm(cross_product(lat(:,3),lat(:,1)))
    r_mesh(1) = N**(1/3.0) * w**(1/3.0) * u**(1/3.0) / v**(2/3.0)
    r_mesh(2) = N**(1/3.0) * u**(1/3.0) * v**(1/3.0) / w**(2/3.0)
    r_mesh(3) = N**(1/3.0) * v**(1/3.0) * w**(1/3.0) / u**(2/3.0)
    write(*,*) 'n1,n2,n3 from min s/v',r_mesh
	write(*,*)
    r_pqr(1) = r_mesh(2)/r_mesh(1)
    r_pqr(2) = r_mesh(3)/r_mesh(2)
    r_pqr(3) = r_mesh(1)/r_mesh(3)
    write(*,*) 'real mesh ratios', r_pqr
	write(*,*)
!   Define triangle as
!               m1
!            P      R
!
!        m2     Q     m3
!    The multiple relations (P,Q,R) 'point' to the larger number.  If they are CCW, record +.
!    If CW, -
	pqr = 0
!    write(*,*) 'initial pqr', pqr
    if (intcheck(r_pqr(1))) then
    	pqr(1) = -nint(r_pqr(1))
    end if
    if (intcheck(1/r_pqr(1))) then
        pqr(1) = nint(1/r_pqr(1))
    end if
    if (intcheck(r_pqr(2))) then
    	pqr(2) = nint(r_pqr(2))
    end if
    if (intcheck(1/r_pqr(2))) then
        pqr(2) = -nint(1/r_pqr(2))
    end if
    if (intcheck(r_pqr(3))) then
    	pqr(3) = nint(r_pqr(3))
    end if
    if (intcheck(1/r_pqr(3))) then
        pqr(3) = -nint(1/r_pqr(3))
    end if
    write(*,*) 'integer mesh ratios', pqr
	nrels = 0 !number of integer relations
	do i = 1,3
	    if (pqr(i) /= 0) then
	    nrels = nrels + 1
	    end if
	end do
    if (nrels > 3) then
    	write(*,*),'Error in SVMESH'
    end if
    write(*,*),'nrels',nrels
    !form final mesh m's
    if (nrels == 0) then
        do i = 1,3
            i_mesh(i) = nint(r_mesh(i))
        	write(*,*) 'i_mesh','r_mesh',i,i_mesh(i),r_mesh(i)
        end do
    end if
    if (nrels == 1) then
    	write(*,*), 'pqr, maxloc(abs(pqr))', r_pqr
        index1 = maxloc(abs(pqr),1) !want index of the only nonzero element
    	write(*,*) 'maxloc(abs(pqr),1), index1', maxloc(abs(pqr),1), index1
        mratio1 = pqr(index1)
        if (mratio1 > 0) then !CCW, so the higher index is greater m
            i_mesh(index1) = nint(r_mesh(index1)) !round smaller one
            i_mesh(icy(index1,1)) = i_mesh(index1) * mratio1 !form larger by relation
        else !CW
            i_mesh(icy(index1,1)) = nint(r_mesh(icy(index1,1))) !round smaller one
            i_mesh(index1) = i_mesh(icy(index1,1)) * abs(mratio1)      !form larger by relation
           write(*,*) 'test', nint(r_mesh(icy(index1,1))), i_mesh(icy(index1,1)),abs(mratio1)
           write(*,*) 'index1,icy(index1,1))',index1,icy(index1,1)
        end if
        i_mesh(icy(index1,-1)) = nint(r_mesh(icy(index1,-1)))
    end if
    if ((nrels == 2) .or. (nrels == 3)) then
        !find the triangle side with no multiple relation.  Find the largest
        rmax = maxval(abs(r_mesh))
        imax = maxloc(abs(r_mesh),1) !index of largest real mesh number.  Create m from this first
        mratio1 = abs(pqr(icy(imax,-1)))
        mratio2 = abs(pqr(imax))
        i_mesh(imax) = mratio1 * mratio2 * nint(rmax/mratio1/mratio2)
        i_mesh(icy(imax,-1)) = i_mesh(imax)/mratio1
        i_mesh(icy(imax,1)) = i_mesh(imax)/mratio2
    end if
    write(*,*) 'i_mesh', i_mesh
    svmesh = i_mesh
end function svmesh

function intcheck(x) !Checks if a real number is very close to an integer
	logical intcheck
	real x,delta
    delta = 1.0e-6
    write(*,*)'intcheck', x, abs(nint(x)-x)
    if (abs(nint(x)-x) < delta) then
        write(*,*),'true'
        intcheck =.true.
    else
        intcheck =.false.
    end if
end function intcheck

function icy(i,ichange) !for cycling indices 1,2,3.  change is either 1 or -1
	integer icy
	integer i,ichange
    i = i + ichange
    if (i == 4) i=1
    if (i == 0) i=3
    icy = i
end function icy

END MODULE
