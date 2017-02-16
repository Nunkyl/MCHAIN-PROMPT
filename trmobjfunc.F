! PROMPT Toolbox for MATLAB
!
! By Gaik Tamazian, 2016.
! gaik (dot) tamazian (at) gmail (dot) com

#include "fintrf.h"

      ! Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use prompt, only: objFunc, TrModel, Cell_1d_ar, Cell_2d_ar
      implicit none

      mwPointer     :: plhs(*), prhs(*)
      integer       :: nlhs, nrhs

      ! Function declarations
      mwPointer        :: mxCreateDoubleMatrix
      mwPointer        :: mxGetPr, mxGetField
      mwPointer        :: mxGetM, mxGetN
      mwPointer        :: mxSetCell, mxGetCell
      mwPointer        :: mxCreateCellMatrix
      mwpointer        :: mxDuplicateArray
      integer          :: mxIsStruct, mxIsNumeric
      mwSize, external :: mxGetNumberOfElements

      ! Pointers to input data
      type array_of_pointers
          mwPointer :: p
      end type array_of_pointers
      
      type cell_array
          real(kind=8), allocatable, dimension(:,:) :: t
      end type cell_array

      mwPointer               :: model_ptr, angle_indices_ptr, angle_values_ptr,
     + rot_mat_ptr, bond_ind_ptr
      type(array_of_pointers) :: masses_ptr, r_ptr, alpha_ptr, psi_ptr, 
     +  str_crd_ptr, bond_r_ptr, bond_alpha_ptr, bond_psi_ptr
      
      ! Pointers to output data
      mwPointer      :: fval_ptr

      ! Arguments for computational routines
      type(TrModel)                            :: m
      type(array_of_pointers)                  :: atom_num
      mwPointer                                :: conf_num, chain_num !p_num, t_num
      real(kind=8), allocatable, dimension(:)  :: angle_indices
      !real(kind=8), allocatable, dimension(:) :: t_indices
      real(kind=8), allocatable, dimension(:)  :: angle_values
      real(kind=8)                             :: func_val
      !real(kind=8), allocatable, dimension(:) :: grad_vec
      type(cell_array)                         :: temp_coords(:)
      
      ! Check for the proper number of function arguments
      if (nrhs .ne. 3) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:nInput',
     +                         'Three input arguments required.')
      endif

      ! Check input types
      if (mxIsStruct(prhs(1)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonStructure',
     +                         'The first argument must be a structure
     +                          representing a transformation model.')
      endif

      if (mxIsNumeric(prhs(2)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonNumeric',
     +                         'The second argument must be a numeric 
     +                          array.')
      endif

      if (mxIsNumeric(prhs(3)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonNumeric',
     +                         'The third argument must be a numeric 
     +                          array.')
      endif

      ! Get pointers for input data
      model_ptr         = prhs(1)
      angle_indices_ptr = mxGetPr(prhs(2))
      angle_values_ptr  = mxGetPr(prhs(3))
      !angle_values_ptr = mxGetPr(prhs(4))
      
      ! Determine the model parameters: the number of atoms, the
      ! number of configurations and the number of chains 
      conf_num  = mxGetN(mxGetCell(mxGetField(model_ptr, 1, 'r'))
      chain_num = mxGetDimensions(mxGetField(model_ptr, 1, 'm'), 1))
      integer(kind=4) :: i
      do i = 1, chain_num
        atom_num(i)%p  = mxGetM(mxGetField(model_ptr, i, 'r')) + 1
        !atom_num(i)%p  = mxGetM(mxGetCell(mxGetField(model_ptr, i, 'r'), i)) + 1
      end do
      !p_num = mxGetM(prhs(2)) * mxGetN(prhs(2)) ! The amount of planar angles
      !t_num = mxGetM(prhs(3)) * mxGetN(prhs(3)) ! The amount of torsion angles
      a_num = mxGetM(prhs(2)) * mxGetN(prhs(2))  ! The amount of all (torsion) angles
      
      do i = 1, chain_num
        masses_ptr(i)%p  = mxGetField(model_ptr, i, 'm')
        r_ptr(i)%p       = mxGetField(model_ptr, i, 'r')
        alpha_ptr(i)%p   = mxGetField(model_ptr, i, 'alpha')
        psi_ptr(i)%p     = mxGetField(model_ptr, i, 'psi')
        str_crd_ptr(i)%p = mxGetField(model_ptr, i, 'StartCoords')
      end do
      
      do i = 1, chain_num - 1
        bond_r_ptr(i)%p     = mxGetField(model_ptr, i, 'bond_r')
        bond_alpha_ptr(i)%p = mxGetField(model_ptr, i, 'bond_alpha')
        bond_psi_ptr(i)%p   = mxGetField(model_ptr, i, 'bond_psi') 
      end do
      
      rot_mat_ptr    = mxGetField(model_ptr, 1, 'U')
      bond_ind_ptr   = mxGetField(model_ptr, 1, 'bondInd')
      

      ! Create result arrays
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0.0)
      !plhs(2) = mxCreateDoubleMatrix((p_num + t_num) * 
     !+  (conf_num - 2), 1, 0.0)

      ! Get pointers to output data
      fval_ptr = mxGetPr(plhs(1))
      !grad_ptr = mxGetPr(plhs(2))

      ! Allocate Fortran arrays
      do i = 1, chain_num
        allocate(m(i)%atom_num(chain_num))
        allocate(m(i)%atom_masses(atom_num(i)%p))
        allocate(m(i)%r(atom_num(i)%p - 1, conf_num))
        allocate(m(i)%alpha(atom_num(i)%p - 2, conf_num))
        allocate(m(i)%psi(atom_num(i)%p - 3, conf_num))
        allocate(m(i)%start_coords(3, atom_num(i)%p))
        allocate(temp_coords(i)%t(atom_num(i)%p, 3))
      end do
      
      do i = 1, chain_num - 1
        allocate(m(i)%bond_r(1, conf_num))
        allocate(m(i)%bond_alpha(2, conf_num))
        allocate(m(i)%bond_psi(3, conf_num)) 
      end do
      
      allocate(m(1)%rot_mat(3, 3, conf_num - 1))
      allocate(m(1)%bond_ind(chain_num - 1, 4))
      
      
      !allocate(p_indices(p_num))
      !allocate(angle_values((p_num + t_num) * (conf_num - 2)))
      !allocate(t_indices(t_num))
      allocate(angle_indices(a_num))
      allocate(angle_values(a_num))
      !allocate(grad_vec((p_num + t_num) * (conf_num - 2)))

      m(1)%conf_num = int(conf_num, 4)
      !m(1)%atom_num = int(atom_num, 4)
      m(1)%chain_num = int(chain_num, 4)

      ! Copy the input data to Fortran arrays
      do i = 1, chain_num
        call mxCopyPtrToReal8(masses_ptr(i)%p, m(i)%atom_masses, atom_num(i)%p)
        call mxCopyPtrToReal8(r_ptr(i)%p, m(i)%r, (atom_num(i)%p - 1) * conf_num)
        call mxCopyPtrToReal8(alpha_ptr(i)%p, m(i)%alpha, (atom_num(i)%p - 2) * conf_num)
        call mxCopyPtrToReal8(psi_ptr(i)%p, m(i)%psi, (atom_num(i)%p - 3) * conf_num)
        call mxCopyPtrToReal8(str_crd_ptr(i)%p, temp_coords(i)%t, atom_num(i)%p * 3)
        m(i)%start_coords = transpose(temp_coords(i)%t)
        call mxCopyPtrToReal8(rot_mat_ptr, m(1)%rot_mat, conf_num * 3 * 3)
      end do 
      
      !call mxCopyPtrToReal8(p_index_ptr, p_indices, p_num)
      !call mxCopyPtrToReal8(t_index_ptr, t_indices, t_num)
      call mxCopyPtrToReal8(angle_indices_ptr, angle_indices, a_num)
      call mxCopyPtrToReal8(angle_values_ptr, angle_values, a_num * (conf_num - 2))
      
      do i = 1, chain_num
        call mxCopyPtrToReal8(bond_r_ptr(i)%p, m(i)%bond_r, conf_num)
        call mxCopyPtrToReal8(bond_alpha_ptr(i)%p, m(i)%bond_alpha, 2 * conf_num)
        call mxCopyPtrToReal8(bond_psi_ptr(i)%p, m(i)%bond_psi, 3 * conf_num)
      end do 

      ! Launch the computational routines
      call objFunc(m, int(a_num, 4), int(angle_indices, 4),
     !+  int(t_num, 4), int(t_indices, 4),
     +  angle_values, 2) ! nlhs > 1, func_val, grad_vec

      ! Copy the output data back to MATLAB data structures
      call mxCopyReal8ToPtr(func_val, fval_ptr, 1)
      !call mxCopyReal8ToPtr(grad_vec, grad_ptr, (p_num + t_num) *
     !+  (conf_num - 2))

      ! Deallocate Fortran arrays
      do i = 1, chain_num
        deallocate(m(i)%atom_num)
        deallocate(m(i)%atom_masses)
        deallocate(m(i)%r)
        deallocate(m(i)%alpha)
        deallocate(m(i)%psi)
        deallocate(m(i)%start_coords)
        deallocate(temp_coords(i)%t)
      end do
      
      do i = 1, chain_num - 1
        deallocate(m(i)%bond_r)
        deallocate(m(i)%bond_alpha)
        deallocate(m(i)%bond_psi) 
      end do
      
      deallocate(m(1)%rot_mat)
      deallocate(m(1)%bond_ind)
      deallocate(angle_indices)
      deallocate(angle_values)

      return
      end subroutine mexFunction