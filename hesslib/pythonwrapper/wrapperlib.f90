!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_numerical_utilities_equal_rank1 (equal_rank1f2&
     &pywrap, a, b, tolerance, f2py_a_d0, f2py_b_d0)
      use numerical_utilities, only : equal_rank1
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_b_d0
      real(kind=dp) a(f2py_a_d0)
      real(kind=dp) b(f2py_b_d0)
      logical equal_rank1f2pywrap
      equal_rank1f2pywrap = .not.(.not.equal_rank1(a, b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_rank1
      subroutine f2pywrap_numerical_utilities_equal_rank2 (equal_rank2f2&
     &pywrap, a, b, tolerance, f2py_a_d0, f2py_a_d1, f2py_b_d0, f2py_b_d&
     &1)
      use numerical_utilities, only : equal_rank2
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_b_d0
      integer f2py_b_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      real(kind=dp) b(f2py_b_d0,f2py_b_d1)
      logical equal_rank2f2pywrap
      equal_rank2f2pywrap = .not.(.not.equal_rank2(a, b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_rank2
      subroutine f2pywrap_numerical_utilities_equal_rank3 (equal_rank3f2&
     &pywrap, a, b, tolerance, f2py_a_d0, f2py_a_d1, f2py_a_d2, f2py_b_d&
     &0, f2py_b_d1, f2py_b_d2)
      use numerical_utilities, only : equal_rank3
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_a_d2
      integer f2py_b_d0
      integer f2py_b_d1
      integer f2py_b_d2
      real(kind=dp) a(f2py_a_d0,f2py_a_d1,f2py_a_d2)
      real(kind=dp) b(f2py_b_d0,f2py_b_d1,f2py_b_d2)
      logical equal_rank3f2pywrap
      equal_rank3f2pywrap = .not.(.not.equal_rank3(a, b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_rank3
      subroutine f2pywrap_numerical_utilities_equal_scalar (equal_scalar&
     &f2pywrap, a, b, tolerance)
      use numerical_utilities, only : equal_scalar
      real(kind=dp) a
      real(kind=dp) b
      real(kind=dp) tolerance
      logical equal_scalarf2pywrap
      equal_scalarf2pywrap = .not.(.not.equal_scalar(a, b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_scalar
      subroutine f2pywrap_numerical_utilities_equal_scalar_real_int (equ&
     &al_scalar_real_intf2pywrap, a, b, tolerance)
      use numerical_utilities, only : equal_scalar_real_int
      real(kind=dp) a
      integer b
      real(kind=dp) tolerance
      logical equal_scalar_real_intf2pywrap
      equal_scalar_real_intf2pywrap = .not.(.not.equal_scalar_real_int(a&
     &, b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_scalar_real_int
      subroutine f2pywrap_numerical_utilities_equal_scalar_int_int (equa&
     &l_scalar_int_intf2pywrap, a, b, tolerance)
      use numerical_utilities, only : equal_scalar_int_int
      integer a
      integer b
      real(kind=dp) tolerance
      logical equal_scalar_int_intf2pywrap
      equal_scalar_int_intf2pywrap = .not.(.not.equal_scalar_int_int(a, &
     &b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_scalar_int_int
      subroutine f2pywrap_numerical_utilities_equal_rank1_rank0 (equal_r&
     &ank1_rank0f2pywrap, a, b, tolerance, f2py_a_d0)
      use numerical_utilities, only : equal_rank1_rank0
      real(kind=dp) b
      real(kind=dp) tolerance
      integer f2py_a_d0
      real(kind=dp) a(f2py_a_d0)
      logical equal_rank1_rank0f2pywrap
      equal_rank1_rank0f2pywrap = .not.(.not.equal_rank1_rank0(a, b, tol&
     &erance))
      end subroutine f2pywrap_numerical_utilities_equal_rank1_rank0
      subroutine f2pywrap_numerical_utilities_equal_rank2_rank0 (equal_r&
     &ank2_rank0f2pywrap, a, b, tolerance, f2py_a_d0, f2py_a_d1)
      use numerical_utilities, only : equal_rank2_rank0
      real(kind=dp) b
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      logical equal_rank2_rank0f2pywrap
      equal_rank2_rank0f2pywrap = .not.(.not.equal_rank2_rank0(a, b, tol&
     &erance))
      end subroutine f2pywrap_numerical_utilities_equal_rank2_rank0
      subroutine f2pywrap_numerical_utilities_equal_rank2_real_int (equa&
     &l_rank2_real_intf2pywrap, a, b, tolerance, f2py_a_d0, f2py_a_d1, f&
     &2py_b_d0, f2py_b_d1)
      use numerical_utilities, only : equal_rank2_real_int
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_b_d0
      integer f2py_b_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      integer b(f2py_b_d0,f2py_b_d1)
      logical equal_rank2_real_intf2pywrap
      equal_rank2_real_intf2pywrap = .not.(.not.equal_rank2_real_int(a, &
     &b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_rank2_real_int
      subroutine f2pywrap_numerical_utilities_equal_rank1_real_int (equa&
     &l_rank1_real_intf2pywrap, a, b, tolerance, f2py_a_d0, f2py_b_d0)
      use numerical_utilities, only : equal_rank1_real_int
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_b_d0
      real(kind=dp) a(f2py_a_d0)
      integer b(f2py_b_d0)
      logical equal_rank1_real_intf2pywrap
      equal_rank1_real_intf2pywrap = .not.(.not.equal_rank1_real_int(a, &
     &b, tolerance))
      end subroutine f2pywrap_numerical_utilities_equal_rank1_real_int
      
      subroutine f2pyinitnumerical_utilities(f2pysetupfunc)
      interface 
      subroutine f2pywrap_numerical_utilities_equal_rank1 (equal_rank1f2&
     &pywrap, equal_rank1, a, b, tolerance, f2py_a_d0, f2py_b_d0)
      logical equal_rank1
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_b_d0
      real(kind=dp) a(f2py_a_d0)
      real(kind=dp) b(f2py_b_d0)
      logical equal_rank1f2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank1 
      subroutine f2pywrap_numerical_utilities_equal_rank2 (equal_rank2f2&
     &pywrap, equal_rank2, a, b, tolerance, f2py_a_d0, f2py_a_d1, f2py_b&
     &_d0, f2py_b_d1)
      logical equal_rank2
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_b_d0
      integer f2py_b_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      real(kind=dp) b(f2py_b_d0,f2py_b_d1)
      logical equal_rank2f2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank2 
      subroutine f2pywrap_numerical_utilities_equal_rank3 (equal_rank3f2&
     &pywrap, equal_rank3, a, b, tolerance, f2py_a_d0, f2py_a_d1, f2py_a&
     &_d2, f2py_b_d0, f2py_b_d1, f2py_b_d2)
      logical equal_rank3
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_a_d2
      integer f2py_b_d0
      integer f2py_b_d1
      integer f2py_b_d2
      real(kind=dp) a(f2py_a_d0,f2py_a_d1,f2py_a_d2)
      real(kind=dp) b(f2py_b_d0,f2py_b_d1,f2py_b_d2)
      logical equal_rank3f2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank3 
      subroutine f2pywrap_numerical_utilities_equal_scalar (equal_scalar&
     &f2pywrap, equal_scalar, a, b, tolerance)
      logical equal_scalar
      real(kind=dp) a
      real(kind=dp) b
      real(kind=dp) tolerance
      logical equal_scalarf2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_scalar 
      subroutine f2pywrap_numerical_utilities_equal_scalar_real_int (equ&
     &al_scalar_real_intf2pywrap, equal_scalar_real_int, a, b, tolerance&
     &)
      logical equal_scalar_real_int
      real(kind=dp) a
      integer b
      real(kind=dp) tolerance
      logical equal_scalar_real_intf2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_scalar_real_int 
      subroutine f2pywrap_numerical_utilities_equal_scalar_int_int (equa&
     &l_scalar_int_intf2pywrap, equal_scalar_int_int, a, b, tolerance)
      logical equal_scalar_int_int
      integer a
      integer b
      real(kind=dp) tolerance
      logical equal_scalar_int_intf2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_scalar_int_int 
      subroutine f2pywrap_numerical_utilities_equal_rank1_rank0 (equal_r&
     &ank1_rank0f2pywrap, equal_rank1_rank0, a, b, tolerance, f2py_a_d0)
      logical equal_rank1_rank0
      real(kind=dp) b
      real(kind=dp) tolerance
      integer f2py_a_d0
      real(kind=dp) a(f2py_a_d0)
      logical equal_rank1_rank0f2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank1_rank0 
      subroutine f2pywrap_numerical_utilities_equal_rank2_rank0 (equal_r&
     &ank2_rank0f2pywrap, equal_rank2_rank0, a, b, tolerance, f2py_a_d0,&
     & f2py_a_d1)
      logical equal_rank2_rank0
      real(kind=dp) b
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      logical equal_rank2_rank0f2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank2_rank0 
      subroutine f2pywrap_numerical_utilities_equal_rank2_real_int (equa&
     &l_rank2_real_intf2pywrap, equal_rank2_real_int, a, b, tolerance, f&
     &2py_a_d0, f2py_a_d1, f2py_b_d0, f2py_b_d1)
      logical equal_rank2_real_int
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_a_d1
      integer f2py_b_d0
      integer f2py_b_d1
      real(kind=dp) a(f2py_a_d0,f2py_a_d1)
      integer b(f2py_b_d0,f2py_b_d1)
      logical equal_rank2_real_intf2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank2_real_int 
      subroutine f2pywrap_numerical_utilities_equal_rank1_real_int (equa&
     &l_rank1_real_intf2pywrap, equal_rank1_real_int, a, b, tolerance, f&
     &2py_a_d0, f2py_b_d0)
      logical equal_rank1_real_int
      real(kind=dp) tolerance
      integer f2py_a_d0
      integer f2py_b_d0
      real(kind=dp) a(f2py_a_d0)
      integer b(f2py_b_d0)
      logical equal_rank1_real_intf2pywrap
      end subroutine f2pywrap_numerical_utilities_equal_rank1_real_int
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_numerical_utilities_equal_rank1,f2pywr&
     &ap_numerical_utilities_equal_rank2,f2pywrap_numerical_utilities_eq&
     &ual_rank3,f2pywrap_numerical_utilities_equal_scalar,f2pywrap_numer&
     &ical_utilities_equal_scalar_real_int,f2pywrap_numerical_utilities_&
     &equal_scalar_int_int,f2pywrap_numerical_utilities_equal_rank1_rank&
     &0,f2pywrap_numerical_utilities_equal_rank2_rank0,f2pywrap_numerica&
     &l_utilities_equal_rank2_real_int,f2pywrap_numerical_utilities_equa&
     &l_rank1_real_int)
      end subroutine f2pyinitnumerical_utilities

      subroutine f2pywrap_vector_matrix_utilities_minkowski_conditions_c&
     &heck (minkowski_conditions_checkf2pywrap, basis, eps)
      use vector_matrix_utilities, only : minkowski_conditions_check
      real(kind=dp) eps
      real(kind=dp) basis(3,3)
      logical minkowski_conditions_checkf2pywrap
      minkowski_conditions_checkf2pywrap = .not.(.not.minkowski_conditio&
     &ns_check(basis, eps))
      end subroutine f2pywrap_vector_matrix_utilities_minkowski_conditio&
     &ns_check
      
      subroutine f2pyinitvector_matrix_utilities(f2pysetupfunc)
      use vector_matrix_utilities, only : reduce_C_in_ABC
      interface 
      subroutine f2pywrap_vector_matrix_utilities_minkowski_conditions_c&
     &heck (minkowski_conditions_checkf2pywrap, minkowski_conditions_che&
     &ck, basis, eps)
      logical minkowski_conditions_check
      real(kind=dp) eps
      real(kind=dp) basis(3,3)
      logical minkowski_conditions_checkf2pywrap
      end subroutine f2pywrap_vector_matrix_utilities_minkowski_conditio&
     &ns_check
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_vector_matrix_utilities_minkowski_cond&
     &itions_check,reduce_C_in_ABC)
      end subroutine f2pyinitvector_matrix_utilities


