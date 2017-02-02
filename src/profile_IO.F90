#include "intrinsic_sizes.h"
#include "redef.h"
!>Read profiles from ASCII files
module profile_io
  use par_mod
  use lagrange_interpolation
  use file_io, only: get_unit_nr, get_nr_of_lines,&
       &get_nr_of_substr
  use par_geom, only: geomdir, geomfile, magn_geometry, x_def
  implicit none

  public:: calc_rhostar, mref, Tref, nref, omegatorref
  public:: iterdb_file, iterdb_time, norm_index
  public:: get_entrysize_gene, read_profiles_gene
  public:: get_entrysize_iterdb, read_profiles_iterdb
  public:: get_entrysize_iterdb_d3d, read_profiles_iterdb_d3d
#ifdef WITHFUTILS
  public:: get_entrysize_chease, read_profiles_chease
#endif
  !moved to par_mod.F90
  !public:: m_proton, m_deuteron, m_electron
  !REAL, PARAMETER :: m_proton   = 1.672621637E-27 !kg
  !REAL, PARAMETER :: m_deuteron = 3.34358320E-27 !kg
  !REAL, PARAMETER :: m_electron = 9.10938215E-31 !kg
  integer:: norm_index=-1 !<index of species for temp/dens norm.
                          !(-1: search for electrons)
  character(len=FILENAME_MAX) :: iterdb_file=''
  real:: iterdb_time=-1.0
  real:: mref, Tref=0., nref=0., omegatorref=0.0
  private

contains


  !>Routine for automatic calculation of rhostar from profile and geometry input. 
  !!If the simulation considers kinetic electrons, rhostar will be rho_s/a, otherwise
  !!it will be set to rho_i/a (using the first ion species).
  subroutine calc_rhostar(Bref,Lref,minor_r,rhostar)

    real, intent(in):: Bref, Lref, minor_r
    real, intent(inout):: rhostar
    integer:: n, nt_norm_index, nx0_file
    real, dimension(:),allocatable:: temp_prof, dens_prof, xaxis_file, omt_prof,&
         !for chease
         te_prof, ne_prof, ti_prof, ni_prof
    real :: old_Tref
    character(len=50):: entry
    type(spectype) :: sp

    if (norm_index.eq.-1) then
       nt_norm_index = 0 !normalize to first species if 
                         !electrons are not found
       
       do n=0, n_spec-1
          if (spec(n)%charge.lt.0) nt_norm_index = n
       enddo
    else
       nt_norm_index = norm_index
    endif

    if (mref.lt.0.0) then
       if (norm_index.ge.0) then
          mref = spec(norm_index)%mass
          do n=0, n_spec-1
             spec(n)%mass = spec(n)%mass/mref
          enddo
       else
          mref = m_deuteron/m_proton
       endif
       if (mype==0)  Write(*,'(A,E13.6)') 'mref/mproton set to  ',mref
    endif

    if (spec(nt_norm_index)%prof_type.gt.0) then !analytical profile
       sp = spec(nt_norm_index)
       if (Tref.lt.0.0) stop 'Tref has to be given at center position'
       old_Tref = Tref
       select case(spec(nt_norm_index)%prof_type)
       case(1)
          Tref = Tref*exp(sp%kappa_T*minor_r/(sinh(sp%LT_center/sp%LT_width)**2)*(&
               x0-sp%LT_center-sp%LT_width*cosh(sp%LT_center/sp%LT_width)**2*&
               tanh((x0-sp%LT_center)/sp%LT_width)))
       case(2)
          Tref = Tref*exp(-sp%kappa_T*minor_r*sp%LT_width*tanh((x0-sp%LT_center)/sp%LT_width))
       case(3)
          Tref = Tref*(cosh((x0-sp%LT_center+sp%delta_x_T)/sp%LT_width)/&
               &cosh((x0-sp%LT_center-sp%delta_x_T)/sp%LT_width))**&
               &(-0.5*sp%kappa_T*minor_r*sp%LT_width)
       case(4) !Falchetto et al., PPCF 50, 124015
          Tref = Tref*exp(-sp%kappa_T*minor_r*(x0-sp%LT_center-sp%LT_width*&
               &(tanh((x0-sp%LT_center-0.5*sp%delta_x_T)/sp%LT_width)+&
               &tanh((x0-sp%LT_center+0.5*sp%delta_x_T)/sp%LT_width))))
       case(5) !exponential profiles
          Tref = Tref*exp(-sp%kappa_T*minor_r*(x0-sp%LT_center))
       case default
          stop 'rhostar computation not implemented for this profile type!'
       end select
    elseif (Tref.lt.0.0) then
       select case(spec(nt_norm_index)%prof_type)
       case(-1)
          
          if (trim(spec(nt_norm_index)%prof_file).eq.'') &
               &spec(nt_norm_index)%prof_file='./profiles_'//trim(spec(nt_norm_index)%name)
          
          call get_entrysize_gene(spec(nt_norm_index)%prof_file,nx0_file)
          allocate(xaxis_file(1:nx0_file), temp_prof(1:nx0_file), dens_prof(1:nx0_file))
          call read_profiles_gene(spec(nt_norm_index)%prof_file, nx0_file, temp_prof, dens_prof, xaxis_file)
          
       case(-2)
          
          if (spec(nt_norm_index)%charge.lt.0) then
             entry='TE'
          else
             entry='TI'
          endif
          
          !Here, we first read the grid size for the temperature, so we can preallocate
          !arrays with the correct sizes in order to get the file data directly (no interpolation).
          !read_profiles_iterdb knows that it should do this if xaxis_file==0.
          call get_entrysize_iterdb(iterdb_file, iterdb_time, entry, nx0_file)
          allocate(xaxis_file(1:nx0_file), temp_prof(1:nx0_file), dens_prof(1:nx0_file), omt_prof(1:nx0_file)) 
          xaxis_file=0.0
          call read_profiles_iterdb(iterdb_file, iterdb_time, entry, nx0_file, temp_prof, omt_prof, xaxis_file)
          temp_prof=temp_prof*1.e-3 !convert to GENE convention (keV)
          
       case(-3)
          
          if (spec(nt_norm_index)%charge.lt.0) then
             entry='electron temperature, keV'
          else
             !typo as in the D3D IterDB files
             entry='ion temperatue, keV'
          endif
          
          call get_entrysize_iterdb_d3d(iterdb_file, nx0_file)
          allocate(xaxis_file(1:nx0_file), temp_prof(1:nx0_file), dens_prof(1:nx0_file))       
          call read_profiles_iterdb_d3d(iterdb_file, entry, nx0_file, temp_prof, xaxis_file)
          
       case(-4)
          
#ifdef WITHFUTILS
          call get_entrysize_chease(nx0_file)
          allocate(xaxis_file(1:nx0_file),te_prof(1:nx0_file),ne_prof(1:nx0_file),ti_prof(1:nx0_file),&
               ni_prof(1:nx0_file),temp_prof(1:nx0_file),dens_prof(1:nx0_file))       
          call read_profiles_chease(nx0_file,xaxis_file,te_prof,ti_prof,ne_prof,ni_prof)
          if (spec(nt_norm_index)%charge.lt.0) then
             temp_prof=te_prof*1.e-3 !convert to GENE convention (keV)
          else
             temp_prof=ti_prof*1.e-3 !convert to GENE convention (keV)
          endif
          xaxis_file=xaxis_file/maxval(xaxis_file)
#else 
          print*, "You have to compile with FUITLS switch set to yes to use CHEASE profiles"  
          stop
#endif
          
       case(0)
          stop 'rhostar computation only possible with given Tref!'

       case default
          stop 'rhostar computation not implemented for this profile type!'
          
       end select
       
       !interpolate to reference position
       call lag3interp(temp_prof,xaxis_file,nx0_file,Tref,x0)
       if (mype==0) write(*,'(A,E13.6)') 'Tref/keV set to       ',Tref
       deallocate(xaxis_file,temp_prof,dens_prof)
    endif
    
    if (.not.any((/Tref,mref,Bref,Lref,minor_r/).le.0.0)) then
       if (rhostar.lt.0.0) then
          rhostar= 3.2255E-3 * sqrt((mref)*Tref)/Bref/(minor_r*Lref)
          if (mype.eq.0) write (*,"(A,ES12.5,A,F8.3)") &
               &'Using profile and geometry input to compute rhostar. Result: ',&
               &rhostar,' = 1 /',1./rhostar
       endif
    else
       stop 'Cannot compute rhostar due to missing reference values!'
    endif

    if(spec(nt_norm_index)%prof_type.gt.0) Tref = old_Tref !need to revert Tref
 
  end subroutine calc_rhostar



  !>Get number of radial gridpoints in GENE profile file
  subroutine get_entrysize_gene(filename,nx0_file)

    Character(Len=*),intent(in)  :: filename
    integer, intent(out):: nx0_file
    INTEGER :: fileunit, iostat, count

    call get_unit_nr(fileunit)
    count=10
    iostat=1
    do while (iostat.ne.0.and.count.gt.0)
       OPEN(fileunit,file=trim(filename),iostat=iostat)
       count=count-1
    enddo

    IF (iostat.ne.0) then
       WRITE(*,"(3A)") "ERROR reading profile file ",trim(filename)
       STOP
    ENDIF

    nx0_file = get_nr_of_lines(fileunit)
    !skip header lines
    nx0_file = nx0_file - 2
    if (nx0_file.lt.1) then
       WRITE(*,"(3A)") "ERROR: found less than one radial grid point in ",trim(filename) 
       STOP
    endif
    close(fileunit)

  end subroutine get_entrysize_gene



  !>Reads ASCII files containing temperature and density profiles
  !!The format should be the following:
  !!1.) Two header lines beginning with the # symbol
  !!    where the first line should contain as many string entries
  !!    (without space characters) as data columns are present.
  !!    The second line is left for comments (e.g., time stamp)
  !!2.) The data should be arranged as follows:
  !!    1.column: x-axis using the same definition as GENE
  !!    2.column: arbitrary x-axis (e.g. for plotting) which is 
  !!              ignored
  !!    3.column: temperature in keV
  !!    4.column: density in 1E19 m^{-3}
  !!    + as many dummy columns as you wish.
  !!3.) The file may contain several data blocks (e.g., time steps)
  !!    being separated by two blank lines + one comment line
  subroutine read_profiles_gene(filename,nx0_file,temp_file,dens_file,xaxis_file)

    Character(Len=*),intent(in)  :: filename
    integer, intent(in):: nx0_file
    INTEGER :: ncols, fileunit, iostat, count
    REAL, dimension(:,:),allocatable :: data_arr
    REAL, dimension(1:nx0_file) :: temp_file, dens_file,&
         xaxis_file
    LOGICAL :: read_entry=.false.
    Character(Len=128) :: header
    INTEGER :: tstep, tstep_in = 1 !tstep_in defines data block to be read

    call get_unit_nr(fileunit)
    count=10
    iostat=1
    do while (iostat.ne.0.and.count.gt.0)
       OPEN(fileunit,file=trim(filename),iostat=iostat)
       count=count-1
    enddo

    IF (iostat.ne.0) then
       WRITE(*,"(3A)") "ERROR reading profile file ",trim(filename)
       STOP
    ENDIF

    !read and ignore header lines which just contain general
    !information
    read(fileunit,"(A)",iostat=iostat) header
    !determine number of columns
    ncols = get_nr_of_substr(header)
    !ignore '#' symbol:
    ncols = ncols - 1

    if (ncols.lt.4) stop 'Less than 4 columns in profile file'

    allocate(data_arr(1:ncols,1:nx0_file))
    tstep=0
    do while ((iostat.eq.0).and.(tstep.lt.tstep_in))
       read(fileunit,"(A)",iostat=iostat) header
       read (fileunit,*,iostat=iostat) data_arr
       tstep=tstep+1
       if (tstep.lt.tstep_in) then
          !read the two blank lines
          read(fileunit,"(A)",iostat=iostat) header
          read(fileunit,"(A)",iostat=iostat) header
       else
          read_entry = .true.
       endif
    end do

    close(fileunit)
    if (.not.read_entry) then
       WRITE(*,"(3A)") "Could not read profile from ",trim(filename)
       STOP
    ENDIF

    xaxis_file = data_arr(1,:)
    temp_file  = data_arr(3,:)
    dens_file  = data_arr(4,:)

    deallocate(data_arr)

  end subroutine read_profiles_gene


  !>Read number of radial gridpoints from IterDB file
  subroutine get_entrysize_iterdb(infile, target_time, var_name, nrho_db)

    character (*), intent (in) :: infile
    real, intent (in) :: target_time
    character (*), intent (in) :: var_name
    integer :: fileunit, it, ierr, ntime
    integer, intent(inout) :: nrho_db
    character (10) :: var
    character (1000) :: line
    logical:: found

    found = .false.
    
    call get_unit_nr(fileunit)
    open (unit=fileunit, file=trim(infile), status="old", action="read",&
         & iostat=ierr)
    if (ierr.ne.0) then
       write(*,'(3A)') "ERROR: file ",trim(infile)," not found!"
       STOP
    endif

    line = ''
    ierr = 0.0
    do while (.not. found)
       !skip potential unofficial header and comment lines
       do while ((INDEX(line,';-SHOT ').LE.0).and.(ierr.EQ.0)) 
          read (fileunit,'(A)',IOSTAT=ierr) line
       enddo

       ! if end of file, print error message and stop if variable is mandatory
       if (ierr < 0) then
          write (*,*) "Variable ", trim(var_name), " not found in ", TRIM(infile)
       end if

       do it=1,4
          read (fileunit,'(A)') line
       enddo
       read (fileunit,'(2A)') var, line
       read (fileunit,'(A)') line
       read (fileunit,*) nrho_db
       read (fileunit,*) ntime
       if (trim(adjustl(var)).EQ.trim(var_name)) then
          found = .true.
       endif
    enddo

    close(fileunit)

  end subroutine get_entrysize_iterdb


  !>Read variable defined by var_name from ITER data base file
  !! \param infile string containing the file name
  !! \param target_time data will be read from time stamp being closest to target_time
  !! \param var_name variable to be read, see http://tokamak-profiledb.ccfe.ac.uk/DOCS/pr08vardoc.htm
  !! \param rho GENE's rho_tor grid for profile representation
  !! \param var_out variable on GENE grid
  !! \param grad_out negative logarithmic gradient of variable on GENE grid
  !! \param declare variable as optional (no stop if not found)
  subroutine read_profiles_iterdb(infile, target_time, var_name, nrho, &
       &var_out, grad_out, rho, opt_var)

    character (*), intent (in) :: infile
    real, intent (in) :: target_time
    character (*), intent (in) :: var_name
    real, dimension (:), intent (out) :: var_out, grad_out
    real, dimension (:), intent (inout) :: rho
    logical, intent(in), optional:: opt_var
    integer, intent(inout):: nrho
    logical :: found
    integer :: fileunit, it, ierr, r, l, nrho_db
    integer:: ntime !,shot
    integer, dimension (1) :: arraymin
    character (10) :: var
    character (14) :: fmt !11
    character (1000) :: line
    real, dimension (:), allocatable :: rho_db, time, vartmp_db, var1d, tmp_arr
    real, dimension (:,:), allocatable :: var_db

    found = .false.
    
    call get_unit_nr(fileunit)
    open (unit=fileunit, file=trim(infile), status="old", action="read",&
         & iostat=ierr)
    if (ierr.ne.0) then
       write(*,'(3A)') "ERROR: file ",trim(infile)," not found!"
       STOP
    endif

    line = ''
    ierr = 0.0
    do while (.not. found)
       !skip potential unofficial header and comment lines
       do while ((INDEX(line,';-SHOT ').LE.0).and.(ierr.EQ.0)) 
          read (fileunit,'(A)',IOSTAT=ierr) line
       enddo

       ! if end of file, print error message and stop if variable is mandatory
       if (ierr < 0) then
          write (*,*) "Variable ", trim(var_name), " not found in ", TRIM(infile)
          if (present(opt_var)) then
             if (.not.opt_var) then
                stop
             else
                var_out = 0.0
                grad_out = 0.0
                exit
             endif
          else
             stop
          endif
       end if

       do it=1,4
          read (fileunit,'(A)') line
       enddo
       read (fileunit,'(2A)') var, line
       read (fileunit,'(A)') line
       read (fileunit,*) nrho_db
       read (fileunit,*) ntime

       allocate (time(ntime))
       allocate (rho_db(nrho_db))
       allocate (var1d(nrho_db))
       allocate (var_db(nrho_db,ntime))
       allocate (vartmp_db(nrho_db*ntime))

       !the iter data base uses the normalized rho_tor
       !see http://tokamak-profiledb.ccfe.ac.uk/DOCS/PR08MAN/pdbman.html
       read (fileunit,*,iostat=ierr) rho_db
       if (ierr.lt.0) then
          write(*,*) 'failed reading rho grid for ',TRIM(var), 'in ',infile
          stop
       endif
       read (fileunit,*,iostat=ierr) time
       if (ierr.lt.0) then
          write(*,*) 'failed reading time for ',TRIM(var),' in ', infile
          stop
       endif

       do l = 1, (nrho_db*ntime)/6
          read (fileunit,'(A1,6(ES13.6))',iostat=ierr) line, vartmp_db(1+(l-1)*6:l*6)
          if (ierr.lt.0) then
             write(*,*) 'failed reading ',TRIM(var),' in line ',l,' of the block in ', infile
             stop
          endif
       enddo
       r = MODULO(nrho_db*ntime,6)
       if (r.NE.0) then
          write(fmt,"(a,i1,a)") '(A1,',r,'(ES13.6))'
          read (fileunit,fmt,iostat=ierr) line, vartmp_db(1+(l-1)*6:(l-1)*6+r)
       endif

       if (ierr.lt.0) then
          write(*,*) 'failed reading data for ',var, 'in ',infile
          stop
       endif

       call cast2darr1d(vartmp_db,var_db,nrho_db*ntime) !type cast
       deallocate(vartmp_db)

       if (trim(adjustl(var)).EQ.trim(var_name)) then
          arraymin = minloc(abs(time-target_time))
          var1d=var_db(:,arraymin(1))
          found = .true.
       endif

       read (fileunit,'(A)') line
       read (fileunit,'(A)') line
       read (fileunit,'(A)') line

       if (found) then
          if (all(rho==0.)) then 
             !in this case we don't know the GENE grid yet, and we just
             !deliver the profiles as present in the file
             nrho=nrho_db
             var_out=var1d
             grad_out=0.0
             rho=rho_db
          else
             call lag3interp(var1d,rho_db,nrho_db,var_out,rho,nrho)
             if (all(var1d.ge.0.0)) then
                allocate(tmp_arr(nrho_db))
                tmp_arr = -log(var1d)
                call lag3deriv(tmp_arr,rho_db,nrho_db,&
                     grad_out,rho,nrho)
             else
                allocate(tmp_arr(nrho))
                call lag3deriv(var1d,rho_db,nrho_db,&
                     tmp_arr,rho,nrho)
                grad_out = -tmp_arr/var_out
             endif
             deallocate(tmp_arr)
          endif
       endif
       deallocate(rho_db,var1d)
       deallocate (time, var_db)

    end do

    close (fileunit)



  end subroutine read_profiles_iterdb

  !>Read number of radial gridpoints from D-IIID IterDB file
  subroutine get_entrysize_iterdb_d3d(infile, nrho_db)

    character (*), intent (in) :: infile
    integer :: fileunit, ierr
    integer, intent(out) :: nrho_db
    character(len=80) :: line=''

    call get_unit_nr(fileunit)
    open (unit=fileunit, file=trim(infile), status="old", action="read",&
         & iostat=ierr)
    if (ierr.ne.0) then
       write(*,'(3A)') "ERROR: file ",trim(infile)," not found!"
       stop
    endif
    do while (index(line,'nj').eq.0)
       read(fileunit, "(A)",end=100,ERR=100) line
    end do
    read(fileunit,"(I12)",end=100,ERR=100) nrho_db
    close(fileunit)
    line=''
    return

100 stop 'STOP (d3d_iterdb): file i/o problem in d3d iterdb file'

  end subroutine get_entrysize_iterdb_d3d




  !>Read variable defined by var_name from ITER data base file
  !! \param infile string containing the file name
  !! \param var_name variable to be read
  !! \param rho GENE's rho_tor grid for profile representation
  !! \param var_out variable on GENE grid
  !! \param grad_out negative logarithmic gradient of variable on GENE grid
  !! \param declare variable as optional (no stop if not found)
  subroutine read_profiles_iterdb_d3d(infile, var_name, nx0_file, var_db, rho_db, opt_var)

    character (*), intent (in) :: infile
    character (*), intent (in) :: var_name
    integer, intent(in):: nx0_file
    real, dimension (1:nx0_file), intent (out) :: rho_db, var_db
    logical, intent(in), optional:: opt_var

    integer :: fileunit, ierr
    integer :: nrho_db
    CHARACTER(len=80) :: line=''

    call get_unit_nr(fileunit)
    open (unit=fileunit, file=trim(infile), status="old", action="read",&
         & iostat=ierr)
    if (ierr.ne.0) then
       write(*,'(3A)') "ERROR: file ",trim(infile)," not found!"
       STOP
    endif

    DO WHILE (INDEX(line,'nj').EQ.0)
       READ(fileunit, "(A)",END=100,ERR=100) line
    END DO
    READ(fileunit,"(I12)",END=100,ERR=100) nrho_db

    DO WHILE ((INDEX(line,'rho grid, meters').EQ.0).AND.&
         &(INDEX(line,'rho grid, normalized').EQ.0))
       READ(fileunit, "(A)", END=100, ERR=100) line
    END DO
    READ(fileunit,*, END=100, ERR=100) rho_db
    IF (INDEX(line,'normalized').EQ.0) &
         &rho_db = rho_db / MAXVAL(rho_db)

    DO WHILE (INDEX(line,TRIM(var_name)).EQ.0)
       READ(fileunit, "(A)", END=200, ERR=100) line
    END DO
    READ(fileunit,*, END=100, ERR=100) var_db

    close(fileunit)

    RETURN

100   STOP 'STOP (d3d_iterdb): file i/o problem in d3d iterdb file'
200   close(fileunit);

  end subroutine read_profiles_iterdb_d3d


  subroutine cast2darr1d(arr2d,arr1d,nelem)

    integer, intent(in) :: nelem
    real, dimension(1:nelem) :: arr2d, arr1d
    arr1d = arr2d

  end subroutine cast2darr1d


#ifdef WITHFUTILS
  !>Read number of radial gridpoints from Chease file
  subroutine get_entrysize_chease(NPSI)

    use futils
    integer :: hdf5_ioutgyro, npsi

    if (magn_geometry.ne.'chease') then
      print*, "Profile type 4 can be used only with CHEASE geometry input"  
      stop
    endif    
    select case(x_def)
    case('arho_t') 
    case default
      print*, "CHEASE internal profiles can be used only if x_def set to 'arho_t'"  
      stop
    end select
    if (geomfile.eq.'') geomfile='ogyropsi.h5'
    CALL openf(trim(geomdir)//'/'//trim(geomfile),hdf5_ioutgyro)
    CALL getatt(hdf5_ioutgyro,"/data", "NPSI",NPSI)
    CALL closef(hdf5_ioutgyro)

  end subroutine get_entrysize_chease


  !>Read profiles from Chease file
  subroutine read_profiles_chease(npsi,rhot_chease,Te_chease,Ti_chease,ne_chease,ni_chease)

    use futils

    integer :: npsi
    real, dimension(1:npsi), intent(out) :: rhot_chease,Te_chease
    real, dimension(1:npsi):: Ti_chease
    real, dimension(1:npsi):: ne_chease,ni_chease

    integer :: hdf5_ioutgyro

    if (magn_geometry.ne.'chease') then
       print*, "Profile type 4 can be used only with CHEASE geometry input"  
       stop
    endif    
    select case(x_def)
    case('arho_t') 
    case default
       !could be extended in the future
       print*, "CHEASE internal profiles can be used only if x_def set to 'arho_t'"  
       stop
    end select
    if (geomfile.eq.'') geomfile='ogyropsi.h5'
    CALL openf(trim(geomdir)//'/'//trim(geomfile),hdf5_ioutgyro)
    CALL getatt(hdf5_ioutgyro,"/data", "NPSI",NPSI)
    
    CALL getarr(hdf5_ioutgyro, "/data/var1d/rho_tor", rhot_chease)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/ne", ne_chease)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/ni", ni_chease)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/Te", Te_chease)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/Ti", Ti_chease)
    
    CALL closef(hdf5_ioutgyro)

  end subroutine read_profiles_chease
#endif


end module profile_io
