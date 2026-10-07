! Copyright (c) Authors:
! Lukas Hronek
! 2026
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.
!
! Every text message of the library goes through this module.
!   negf_log_unit  per-run diagnostics log (fort.12347): opened by negf_output_init on the
!                  master rank; other ranks get /dev/null unless SMEAGOL_LOG_ALL_RANKS is set
!                  (then fort.12347.<rank>)
!   negf_warn      non-fatal message: host output on the master, stderr on other ranks,
!                  log on every rank
!   negf_abort     fatal message to the host output (master) and stderr, flushed, then
!                  MPI_Abort on the SMEAGOL communicator; with collective=.true. (every rank
!                  calls it) the ranks synchronise after flushing so no buffered output is lost
! Ranks other than the master never write to the host output unit: the host program may
! have connected that unit to its output file on every rank.
module mNegfOutput
#ifdef MPI
  use mpi_siesta, only: MPI_COMM_WORLD
#endif
  implicit none
  private

  integer, parameter, public :: negf_log_unit = 12347
  integer, parameter :: stderr_unit = 0
  integer, save :: out_unit = 6
  integer, save :: comm = -1
  integer, save :: rank = 0
  logical, save :: master = .true.
  logical, save :: log_open = .false.
  logical, save :: log_fresh = .true.

  public :: negf_output_init, negf_output_finalize, negf_master, negf_out_unit
  public :: negf_warn, negf_abort, negf_log_flush

contains

  subroutine negf_output_init(mpi_comm, mpi_rank, host_unit)
    integer, intent(in) :: mpi_comm, mpi_rank
    integer, intent(in), optional :: host_unit
    character(len=16) :: env
    character(len=32) :: fname
    integer :: envlen, envstat, ios
    logical :: opened

    comm = mpi_comm
    rank = mpi_rank
    master = (mpi_rank == 0)
    if (present(host_unit)) out_unit = host_unit
    if (log_open) return

    if (master) then
      fname = 'fort.12347'
    else
      call get_environment_variable('SMEAGOL_LOG_ALL_RANKS', env, envlen, envstat)
      if (envstat == 0 .and. envlen > 0 .and. env(1:1) /= '0') then
        write (fname, '(a,i0)') 'fort.12347.', rank
      else
        fname = '/dev/null'
      end if
    end if
    inquire (unit=negf_log_unit, opened=opened)
    if (opened) close (negf_log_unit)
    if (log_fresh .and. fname /= '/dev/null') then
      open (unit=negf_log_unit, file=trim(fname), status='replace', action='write', iostat=ios)
    else
      open (unit=negf_log_unit, file=trim(fname), status='unknown', position='append', action='write', iostat=ios)
    end if
    log_open = (ios == 0)
    log_fresh = .false.
  end subroutine negf_output_init

  subroutine negf_log_flush()
! checkpoint for the diagnostics log so that a host abort loses at most the current energy point
    if (log_open) flush (negf_log_unit)
  end subroutine negf_log_flush

  subroutine negf_output_finalize()
    if (log_open) then
      flush (negf_log_unit)
      close (negf_log_unit)
      log_open = .false.
    end if
  end subroutine negf_output_finalize

  logical function negf_master()
    negf_master = master
  end function negf_master

  integer function negf_out_unit()
    negf_out_unit = out_unit
  end function negf_out_unit

  logical function is_true(flag)
    logical, intent(in), optional :: flag
    is_true = .false.
    if (present(flag)) is_true = flag
  end function is_true

  subroutine negf_warn(msg, icode, collective)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: icode
    logical, intent(in), optional :: collective
    character(len=len(msg)+24) :: text

    text = msg
    if (present(icode)) write (text, '(a,a,i0)') trim(msg), ', code ', icode
    if (master) then
      write (out_unit, '(2a)') 'SMEAGOL warning: ', trim(text)
      flush (out_unit)
    else if (.not. is_true(collective)) then
      write (stderr_unit, '(a,i0,2a)') 'SMEAGOL warning (rank ', rank, '): ', trim(text)
    end if
    write (negf_log_unit, '(2a)') 'warning: ', trim(text)
    flush (negf_log_unit)
  end subroutine negf_warn

  subroutine negf_abort(msg, icode, collective)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: icode
    logical, intent(in), optional :: collective
    character(len=len(msg)+24) :: text
    integer :: ierr

    text = msg
    if (present(icode)) write (text, '(a,a,i0)') trim(msg), ', code ', icode
    if (master) then
      write (out_unit, '(2a)') 'SMEAGOL ERROR: ', trim(text)
      write (stderr_unit, '(2a)') 'SMEAGOL ERROR: ', trim(text)
    else if (.not. is_true(collective)) then
      write (stderr_unit, '(a,i0,2a)') 'SMEAGOL ERROR (rank ', rank, '): ', trim(text)
    end if
    write (negf_log_unit, '(2a)') 'ERROR: ', trim(text)
    flush (out_unit)
    flush (stderr_unit)
    flush (negf_log_unit)
#ifdef MPI
    if (comm == -1) comm = MPI_COMM_WORLD
    if (is_true(collective)) call MPI_Barrier(comm, ierr)
    call MPI_Abort(comm, 1, ierr)
#endif
    stop 1
  end subroutine negf_abort

end module mNegfOutput
