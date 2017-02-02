#include "redef.h"
!>Computes all eigenvalues and left or right eigenvectors of a general complex matrix
!!
!!This routine is taken (in slightly adapted form) from Mario's Parallel Library (MPL) 
!!written by Mario Sempf (THANKS!!). For a full version which contains solvers for various 
!!eigenvalue problems, please contact msempf AT ipp.mpg.de
subroutine c_eigensolve(wantleft,wantright,wantq,size,loc_matrix, &
     desc_matrix,eigval,loc_eigvec_left,desc_eigvec_left,loc_eigvec_right, &
     desc_eigvec_right,loc_q,desc_q)

!!!!!
  !  Computes eigenvalues and, optionally, all left and/or right eigenvectors
  !  of a complex, dense, non-Hermitian matrix using ScaLAPACK and LAPACK.
  !
  !  The main steps are i) transformation to upper Hessenberg form, ii)
  !  transformation to upper triangular (i.e., complex Schur) form, revealing
  !  the eigenvalues on the diagonal, and iii) computing the eigenvectors of  
  !  the triangular matrix and backtransforming them, using the accumulated
  !  similarity transforms which lead to the Schur form, to obtain the
  !  eigenvectors of the original matrix.
  !
  !  A small part of the code is not parallelized. This part is executed
  !  only if eigenvectors are needed. This should not be too bad if not
  !  much more than, say, 8 or 16 processors are used.
  !
  !  The routine is also well suited for real, nonsymmetric matrices
  !  (converted to complex ones).
  !
!!!!!

  implicit none

  logical wantleft,wantright,wantq
  integer size
  complex loc_matrix(*)
  integer, dimension(9) :: desc_matrix
  complex eigval(size)
  complex loc_eigvec_left(*)
  integer, dimension(9) :: desc_eigvec_left
  complex loc_eigvec_right(*)
  integer, dimension(9) :: desc_eigvec_right
  complex loc_q(*)
  integer, dimension(9) :: desc_q
  logical :: wantq2

!!!!!!
  ! PARAMETERS:
  !
  ! (input)  wantleft         if wantleft=.true., compute left eigenvectors
  !
  ! (input)  wantright        if wantright=.true., compute right eigenvectors
  !
  ! (input)  wantq            if wantq=.true., compute the unitary
  !                           transformation matrix converting the input
  !                           matrix A to Schur form: B = Q^-1 * A * Q,
  !                           where B is in Schur form
  !
  ! (input)  size             global size of the matrix A of which the
  !                           eigenproblem is to be solved
  !
  ! (input)  loc_matrix    the matrix A of which the eigenproblem is to be
  !                        solved. It should be distributed block-cyclically,
  !                        as described by the array descriptor desc_matrix.
  !                        loc_matrix is destroyed on output.
  !                        !!! IMPORTANT NOTE!!! The code has been tested only
  !                        for the case that, for all processors, loc_matrix
  !                        is allocated with the same dimensions
  !                        (loc_rows,loc_cols), where
  !               loc_rows = iceil(iceil(size,blockrows),pg_rows)*blockrows,
  !               loc_cols = iceil(iceil(size,blockcols),pg_cols)*blockcols,
  !                        the global dimension of the distributed matrix is
  !                        (size,size), blockrows and blockcols are the 
  !                        blocking constants of the distribution, pg_rows and
  !                        pg_cols are the number of rows and columns in the
  !                        process grid, and iceil() is a ScaLAPACK function.
  !                                           
  ! (input)  desc_matrix        array descriptor belonging to loc_matrix
  !
  ! (output) eigval             global, non-distributed vector containing the
  !                             eigenvalues (not ordered), obtained by every
  !                             processor
  !                        
  !
  ! (output) loc_eigvec_left    Distributed array which should be allocated in
  !                             the same way as loc_matrix. If wantleft=.true.,
  !                             the columns of the global version of this
  !                             array contain the left eigenvectors
  !                             on output, ordered in the same way as the
  !                             eigenvalues. The normalization of the
  !                             eigenvectors is rather random: the i-th
  !                             component of the i-th vector is equal to 1.
  !                             If wantleft=.false., loc_eigvec_left is not
  !                             referenced.
  !
  ! (input)  desc_eigvec_left   array descriptor belonging to loc_eigvec_left
  !
  ! (output) loc_eigvec_right   Distributed array which should be allocated in
  !                             the same way as loc_matrix. If
  !                             wantright=.true., the columns of the global
  !                             version of this array contain the right
  !                             eigenvectors on output, ordered in the same
  !                             way as the eigenvalues. The normalization of
  !                             the eigenvectors is rather random: the i-th
  !                             component of the i-th vector is equal to 1.
  !                             If wantright=.false., loc_eigvec_right is not
  !                             referenced.
  !
  ! (input)  desc_eigvec_right  array descriptor belonging to loc_eigvec_right
  !
  ! (output) loc_q              if wantq=.true., loc_q contains the
  !                             transformation matrix Q on output, cf. the
  !                             commentary for the parameter wantq. loc_q is
  !                             block-cyclically distributed and has to be
  !                             allocated in the same way as loc_matrix. Even
  !                             if wantq=.false., loc_q has to be alloctated.
  !
  ! (input)  desc_q             array descriptor belonging to loc_q
  !                             
!!!!!

!!!!!
  ! LOCAL VARIABLES
!!!!!

  ! BLACS stuff
  integer iam,nprocs,context,pg_rows,pg_cols,myrow,mycol,info
  integer blockrows,blockcols,loc_rows,loc_cols

  ! auxiliary arrays
  complex,allocatable :: loc_h(:,:)
  complex,allocatable :: work(:),loc_tau(:)
  complex,allocatable :: loc_eigvec_aux(:,:)
  integer, dimension(9) :: desc_eigvec_aux

  ! auxiliary local array descriptors
  integer, dimension(9) :: desc_h,desc_tau

  ! more auxiliary stuff
  integer i,j,k,l,m,lwork,ilo,ihi
  logical myblock

  ! ScaLAPACK functions
  integer iceil
  PERFON('c_eigen')
  ! reconstruct block cyclic distrubution information from array descriptor
  context=desc_matrix(2)    ! the BLACS context
  blockrows=desc_matrix(5)  ! number of rows in the blocks used for
  ! distribution of matrices
  blockcols=desc_matrix(6)  ! number of columns in the blocks used for
  ! distribution of matrices

  ! pg_rows = rows in process grid, pg_cols = colums in process grid,
  ! (myrow,mycol) = processor's coordinates in process grid
  call BLACS_GRIDINFO(context,pg_rows,pg_cols,myrow,mycol)

  ! number of rows in the local part of a distributed array:
  loc_cols=iceil(iceil(size,blockcols),pg_cols)*blockcols

  ! number of columns in the local part of a distributed array:
  loc_rows=iceil(iceil(size,blockrows),pg_rows)*blockrows

  ! iam = processor's number, nprocs = total number of processors
  call BLACS_PINFO(iam,nprocs)

  ! build descriptors for auxiliary arrays
  call DESCINIT(desc_h,size,size,blockrows,blockcols,0,0,context, &
       loc_rows,info)
  call DESCINIT(desc_q,size,size,blockrows,blockcols,0,0,context, &
       loc_rows,info)
  call DESCINIT(desc_tau,1,size-1,blockrows,blockcols,0,0,context, &
       1,info)
  ! allocate auxiliary arrays
  allocate(loc_h(loc_rows,loc_cols))

  !
  ! loc_matrix is transformed to upper Hessenberg form by PCGEHRD. Information
  ! about the corresponding similarity transformation is stored in the
  ! lower part of loc_matrix, and in loc_tau.
  ! 

  ilo=1      ! first matrix row/column to operate on
  ihi=size   ! last matrix row/column to operate on

  PERFON('hessenberg')
  ! workspace query
  allocate(work(1),loc_tau(loc_cols))
  lwork = -1
  call PCGEHRD(size,ilo,ihi,loc_matrix,1,1,desc_matrix,loc_tau,work,lwork, &
       info)
  lwork=2*int(work(1))
  deallocate(work)
  allocate(work(lwork))
  ! the actual PCGEHRD computation:
  call PCGEHRD(size,ilo,ihi,loc_matrix,1,1,desc_matrix,loc_tau,work,lwork, &
       info)

  ! copy loc_matrix into loc_h, will be used later for Schur form computation
  call PCGEMR2D(size,size,loc_matrix,1,1,desc_matrix,loc_h,1,1,desc_h,context)

  do j=1,size             ! make loc_h strictly upper Hessenberg,
     do i=j+2,size        ! otherwise the Schur form computation (PCLAHQR)
        ! will be inaccurate
        call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
             myrow,mycol,i,j,myblock,l,m)
        if(myblock) loc_h(l,m)=(0.0,0.0)  
     enddo
  enddo
  PERFOFF
  !
  ! Next, the information encrypted in the lower part of loc_matrix and in
  ! loc_tau will be extracted to form the transformation matrix q, which
  ! transforms the original matrix to upper Hessenberg form. This is done
  ! by the serial LAPACK subroutine ZUNGHR. Unfortunately, there is no
  ! parallel equivalent (PZUNGHR), so let processor 0 do the ZUNGHR job
  ! alone.
  !
  !q has to be computed if eigenvectors are requested
  wantq2=wantq
  if (wantright .or. wantleft) wantq2=.true.

  PERFON('kopier')
  if(wantq2) then
     !loc_q has been initialized with the unit matrix
     lwork = -1  
     call PCUNMHR('L','N',size,size,ilo,ihi,loc_matrix,1,1,desc_matrix, &
          loc_tau,loc_q,1,1,desc_q,work,lwork,info)
     lwork=int(work(1))
     lwork=max(lwork,size)
     deallocate(work)
     allocate(work(lwork))
     call PCUNMHR('L','N',size,size,ilo,ihi,loc_matrix,1,1,desc_matrix, &
          loc_tau,loc_q,1,1,desc_q,work,lwork,info)
  endif
  PERFOFF
  PERFON('schur1')
  !
  ! Compute complex Schur form by PCLAHQR.
  ! The transformation matrix loc_q will be multiplied with the similarity
  ! transform needed to obtain the Schur form from the Hessenberg form, i.e.,
  ! it will contain the accumulated similarity transform which generates
  ! the Schur form from the original matrix.
  ! loc_h will be replaced by the Schur form.
  !


  ! workspace query

  lwork = -1
  call PCLAHQR(.true.,.true.,size,ilo,ihi,loc_h,desc_h,eigval,1,size, &
       loc_q,desc_q,work,lwork,0,0,info)
  lwork=int(work(1))
  deallocate(work)
  allocate(work(lwork))
PERFOFF

PERFON('schur2')
  ! The actual PCLAHQR computation. Attention: the first two parameters
  ! should be always set to true, otherwise the results cannot be trusted!!!
  call PCLAHQR(.true.,.true.,size,ilo,ihi,loc_h,desc_h,eigval,1,size, &
       loc_q,desc_q,work,lwork,0,0,info)
  PERFOFF

  !
  ! last step: eigenvector computation
  !
  PERFON('eigenvec')
  if((wantleft) .or. wantright) then
     allocate (loc_eigvec_aux(loc_rows,loc_cols))
     if(wantleft) call DESCINIT(desc_eigvec_left,size,size,blockrows, &
          blockcols,0,0,context,loc_rows,info)
     if(wantright) call DESCINIT(desc_eigvec_right,size,size,blockrows, &
          blockcols,0,0,context,loc_rows,info)
     call DESCINIT(desc_eigvec_aux,size,size,blockrows,blockcols,0,0, &
          context,loc_rows,info)

     ! loc_h is in Schur form, i.e., upper triangular, with the eigenvalues
     ! on the diagonal. Eigenvectors of loc_h are computed first by solving
     ! (loc_h - eigval(i))^H * x = 0 (left eigenvectors), or.
     ! (loc_h - eigval(i))   * x = 0 (right eigenvectors).
     ! A particular solution of this underdetermined system is sought by
     ! adding the equation "(i-th component of eigenvector i) = 1" to the
     ! i-th equation of the system.
     ! Afterwards, the eigenvectors of the original matrix are obtained
     ! by multiplication with the transformation matrix loc_q.

     if(wantleft) then   ! left eigenvectors

        do i=1,size      ! loop over eigenvalues/eigenvectors

           ! generate matrix loc_h - eigval(i)
           do j=1,size
              call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
                   myrow,mycol,j,j,myblock,k,l)
              if(myblock) then
                 if(j .eq. i) then      ! equation #i?
                    loc_h(k,l)=(1.0,0.0)! add "1*component i" on left hand side
                 else
                    loc_h(k,l)=eigval(j)-eigval(i)
                 endif
              endif
           enddo

           ! generate a zero right-hand side vector eigvec_aux(:,i)
           call PCSCAL(size,(0.0,0.0),loc_eigvec_aux,1,i,desc_eigvec_aux,1)

           ! set right-hand side to unity for equation #i
           call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
                myrow,mycol,i,i,myblock,k,l)
           if(myblock) loc_eigvec_aux(k,l)=(1.0,0.0)

           ! solve (loc_h - eigval(i))^H * x = 0 with the modified equation #i;
           ! result in right-hand side vector
           call PCTRSV('U','C','N',size,loc_h,1,1,desc_h,loc_eigvec_aux,1,i, &
                desc_eigvec_aux,1)

        enddo   ! loop over eigenvalues/eigenvectors

        ! backtransform: multiply the resulting vectors with loc_q
        call PCGEMM('N','N',size,size,size,(1.0,0.0),loc_q,1,1,desc_q, &
             loc_eigvec_aux,1,1,desc_eigvec_aux,(0.0,0.0),loc_eigvec_left, &
             1,1,desc_eigvec_left)

     endif   ! left eigenvectors done

     if(wantright) then  ! right eigenvectors

        do i=1,size   ! loop over eigenvalues/eigenvectors

           ! generate loc_h - eigval(i)
           do j=1,size
              call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
                   myrow,mycol,j,j,myblock,k,l)
              if(myblock) then
                 if(j .eq. i) then      ! equation #i?
                    loc_h(k,l)=(1.0,0.0)! add "1*component i" on left hand side
                 else
                    loc_h(k,l)=eigval(j)-eigval(i)
                 endif
              endif
           enddo

           ! generate a zero right-hand side vector eigvec_aux(:,i)
           call PCSCAL(size,(0.0,0.0),loc_eigvec_aux,1,i,desc_eigvec_aux,1)

           ! set right-hand side to unity for equation #i
           call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
                myrow,mycol,i,i,myblock,k,l)
           if(myblock) loc_eigvec_aux(k,l)=(1.0,0.0)

           ! solve (loc_h - eigval(i)) * x = 0, with the modified equation #i;
           ! result in right-hand side vector
           call PCTRSV('U','N','N',size,loc_h,1,1,desc_h,loc_eigvec_aux,1,i, &
                desc_eigvec_aux,1)

        enddo   ! loop over eigenvalues/eigenvectors

        ! backtransform: multiply the resulting vectors with loc_q
        call PCGEMM('N','N',size,size,size,(1.0,0.0),loc_q,1,1,desc_q, &
             loc_eigvec_aux,1,1,desc_eigvec_aux,(0.0,0.0),loc_eigvec_right, &
             1,1,desc_eigvec_right)

     endif   ! right eigenvectors done

     deallocate(loc_eigvec_aux)
#if 0
     ! restore Schur form of loc_h
     do j=1,size
        call loc_index(blockrows,blockcols,pg_rows,pg_cols, &
             myrow,mycol,j,j,myblock,k,l)
        if(myblock) loc_h(k,l)=eigval(j)
     enddo
#endif
  endif ! if(wantleft .or. wantright)) 

  ! Now, the left/right eigenvectors are contained in
  ! loc_eigvec_left/loc_eigvec_right, if wanted. Done.

  ! copy loc_h into loc_matrix, so that the latter is returned in Schur
  ! form
!  call PCGEMR2D(size,size,loc_h,1,1,desc_h,loc_matrix,1,1,desc_matrix,context) 
  PERFOFF
  deallocate(loc_h,work,loc_tau)
  PERFOFF
  return

end subroutine c_eigensolve


!>Auxiliary routine to compute local indices from global from global ones
!!
!!This routine is taken (in slightly adapted form) from Mario's Parallel Library (MPL) 
!!written by Mario Sempf (THANKS!!). For a full version which contains solvers for various 
!!eigenvalue problems, please contact msempf AT ipp.mpg.de
subroutine loc_index(blockrows,blockcols,pg_rows,pg_cols, &
     pgrow,pgcol,glob_i,glob_j,myblock,loc_i,loc_j)

  ! This routine computes the row and column index in the local part of
  ! a distributed array, given the global row and column index, the process
  ! grid coordinates of a process under consideration, and the blocking
  ! constants. It is also checked if that particular process is actually
  ! 'responsible' for storing the element having the given global indices, and
  ! the output of the routine is meaningful only if this is so.

  ! NOTE: The routine assumes that the local index pair (1,1) for processor
  ! (0,0) corresponds to the global index pair (1,1).

  implicit none

  integer blockrows,blockcols,pg_rows,pg_cols
  integer pgrow,pgcol,glob_i,glob_j
  logical myblock
  integer loc_i,loc_j

!!!!!!
  ! PARAMETERS:
  !
  ! (input)  blockrows     number of rows in the blocks used for block-cyclic
  !                        distribution of matrices
  !
  ! (input)  blockcols     number of columns in the blocks used for block-
  !                        cyclic distribution of matrices
  !                        
  ! (input)  pg_rows       number of rows in the process grid
  !
  ! (input)  pg_cols       number of columns in the process grid
  !
  ! (input)  pgrow         row coordinate of the process under consideration
  !
  ! (input)  pgrow         column coordinate of the process under consideration
  !
  ! (input)  glob_i        global row index
  !
  ! (input)  glob_j        global column index
  !
  ! (output) myblock       is set to .true. if an array element with the global
  !                        index pair (glob_i,glob_j) would be 'accessible' to
  !                        the process with the grid coordinate (pgrow,pgcol)
  !                        if the array was block-cyclically distributed
  !                        with the blocking constants blockrows and blockcols;
  !                        otherwise, it is set to .false. .
  !
  ! (output) loc_i         local row index with which the process
  !                        (pgrow,pgcol) would access the distributed array
  !                        element corresponding to the global coordinates
  !                        (glob_i,glob_j); the output is meaningful only if
  !                        myblock=.true. on output; check this before
  !                        interpreting the result
  !
  ! (output) loc_j         local column index with which the process
  !                        (pgrow,pgcol) would access the distributed array
  !                        element corresponding to the global coordinates
  !                        (glob_i,glob_j); the output is meaningful only if
  !                        myblock=.true. on output; check this before
  !                        interpreting the result
  !

  loc_i=mod(glob_i-1,blockrows)+1+((glob_i-1)/(pg_rows*blockrows))*blockrows
  loc_j=mod(glob_j-1,blockcols)+1+((glob_j-1)/(pg_cols*blockcols))*blockcols
  myblock=((mod(glob_i-1,pg_rows*blockrows)/blockrows .eq. pgrow) .and. &
       (mod(glob_j-1,pg_cols*blockcols)/blockcols .eq. pgcol))

  return

end subroutine loc_index

