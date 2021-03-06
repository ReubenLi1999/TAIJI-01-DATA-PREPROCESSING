module taiji_01_preprocessing
    use num_kinds
    use bspline_module

    implicit none

    type gps
        real(kind = 8)                      :: pos(3)
        real(kind = 8)                      :: vel(3)
        real(kind = 8)                      :: acc(3)
        real(kind = 8)                      :: pos_incr(3)
        real(kind = 8)                      :: vel_incr(3)
        real(kind = 8)                      :: time
        real(kind = 8)                      :: rcv_time
        integer(kind = 4)                   :: pos_qualflg
        integer(kind = 4)                   :: vel_qualflg
        integer(kind = 4)                   :: vac_qualflg
        integer(kind = 4)                   :: month
    end type gps

    type io_file
        character(len = 100)                :: name
        integer(kind = 4)                   :: unit
        integer(kind = 8)                   :: nrow
    end type io_file

    type(gps)    , allocatable              :: s(:)
    type(io_file), allocatable              :: i_f(:)
    type(io_file), allocatable              :: o_f(:)
    type(io_file), ALLOCATABLE              :: f_f(:)
    type(io_file)                           :: log_f
    
contains
    
    integer function get_file_n(ifileunit)
        !! from www.fcode.cn
        implicit none
        integer(kind = 4), intent(in   )    :: iFileUnit
        character(len = 1)                  :: cDummy
        integer(kind = 4)                   :: ierr

        get_file_n = 0
        rewind(iFileUnit)

        do
            read(iFileUnit, *, ioStat=ierr) cDummy
            if(ierr /= 0) exit
            get_file_n = get_file_n + 1
        end do

        rewind(iFileUnit)
    end function get_file_n

    subroutine gps_deglitch(me)
        implicit none

        type(gps), intent(inout)            :: me(:)

        INTEGER(kind = 8)                   :: index, i, j

        deglitch_loop: do i = 1, size(me), 1
            if (me(i)%pos_qualflg == 1) then
                me(i: size(me))%pos(1) = me(i: size(me))%pos(1) - me(i)%pos_incr(1)
                me(i: size(me))%pos(2) = me(i: size(me))%pos(2) - me(i)%pos_incr(2)
                me(i: size(me))%pos(3) = me(i: size(me))%pos(3) - me(i)%pos_incr(3)
                me(i: size(me))%vel(1) = me(i: size(me))%vel(1) - me(i)%vel_incr(1)
                me(i: size(me))%vel(2) = me(i: size(me))%vel(2) - me(i)%vel_incr(2)
                me(i: size(me))%vel(3) = me(i: size(me))%vel(3) - me(i)%vel_incr(3)
            end if
        end do deglitch_loop

        write(*, *) 'Deglitching completed...'
        
    end subroutine gps_deglitch

    subroutine diff9(time, vel, acc, sample_rate)
        !!------------------------------------------------------------------------------------------
        !! Description:
        !!      This function is designed to compute the deriavative of an array.
        !! Input:
        !!      time (real64, array): the time epoches for the input signal
        !!      vel (real64, array): the orginal signals;
        !!      acc (real64, array): the first deriavative of the mentioned original signals, and 
        !!                           originally are null
        !!      sample_rate(int32): the sample rate of the input original signal
        !! Output:
        !!      acc (real64, array): the first deriavative of the input original signal;
        !! Algorithm:
        !!      For acc[5: end - 4] part, the nine-points differentiator is applied, however, due to
        !!  the property of this differentiator, the first four points and the last four points is 
        !!  not capable to be computed, therefore, the acc[1: 4] and acc[end - 3: end] are fulfilled
        !!  by the extrapolation method.
        !!------------------------------------------------------------------------------------------
        implicit none
        INTEGER(kind = 4), PARAMETER            :: nx    = 100      !! number of input
        INTEGER(kind = 4), PARAMETER            :: order = 4        !! order in x for interpolation
        INTEGER(kind = 4), PARAMETER            :: iknot = 0        !! automatically select the knots
        INTEGER(kind = 4)                       :: idx              !! derivative of piecewise polynomial to evaluate
        INTEGER(kind = 4)                       :: iflag            !! input array examination flag
        INTEGER(kind = 4)                       :: inbvx            !! initialisation flag
        INTEGER(kind = 4)                       :: err, i

        REAL(kind = 8), INTENT(IN   )           :: sample_rate

        real(kind = 8), INTENT(IN   )           :: time(:)
        real(kind = 8), intent(in   )           :: vel(:)
        real(kind = 8), INTENT(inOUT)           :: acc(:)

        real(kind = 8), ALLOCATABLE             :: reg(:, :)
        real(kind = 8)                          :: tx_bfore(nx + order), tx_after(nx + order)
        real(kind = 8), DIMENSION(3 * order)    :: w1_1d

        LOGICAL                                 :: extrap

        idx = 0
        inbvx = 1
        
        allocate(reg(8, size(vel) - 8), stat=err)
        if (err /= 0) print *, "reg: Allocation request denied"

        reg = 0; acc = 0
        reg(1, :) = vel(1: size(vel) - 8)
        reg(2, :) = vel(2: size(vel) - 7)
        reg(3, :) = vel(3: size(vel) - 6)
        reg(4, :) = vel(4: size(vel) - 5)
        reg(5, :) = vel(6: size(vel) - 3)
        reg(6, :) = vel(7: size(vel) - 2)
        reg(7, :) = vel(8: size(vel) - 1)
        reg(8, :) = vel(9: size(vel))

        !! nine-points differentiator 
        acc(5: size(acc) - 4) = sample_rate / 280.0_wp * (reg(1, :) - 32.0_wp / 3.0_wp * reg(2, :) &
                + 56.0_wp * reg(3, :) - 224.0_wp * reg(4, :) + 224.0_wp * reg(5, :) &
                - 56.0_wp * reg(6, :) + 32.0_wp / 3.0_wp * reg(7, :) - reg(8, :))

        !! extrapolation
        !! initialise
        call db1ink(time(5: nx + 4), nx, acc(5: nx + 4), order, iknot, tx_bfore, acc(5: nx + 4), iflag)
        if (iflag /= 0) then
            write(*, *) 'Error initializing 1D spline: '//get_status_message(iflag)
            stop
        end if

        call db1ink(time(size(time) - nx - 3: size(time) - 4), nx, acc(size(acc) - nx - 3: size(acc) - 4), &
                    order, iknot, tx_after, acc(size(acc) - nx - 3: size(acc) - 4), iflag)
        if (iflag /= 0) then
            write(*, *) 'Error initializing 1D spline: '//get_status_message(iflag)
            stop
        end if

        extrapolation_loop: do i = 1, 4, 1
            extrap = .true.
            call db1val(time(i), idx, tx_bfore, nx, order, acc(5: nx + 4), acc(i), iflag, inbvx, w1_1d, extrap=extrap)
            write(*, *) acc(i)
            call db1val(time(size(time) - i + 1), idx, tx_after, nx, order, &
                        acc(size(acc) - nx - 3: size(acc) - 4), acc(size(acc) - i + 1), iflag, inbvx, w1_1d, extrap=extrap)
        end do extrapolation_loop

        if (allocated(reg)) deallocate(reg, stat=err)
        if (err /= 0) print *, "reg: Deallocation request denied"

        write(*, *) 'Accelaration computing completed...'

    end subroutine diff9

    subroutine detect_missing(time, flag)
        implicit none
        real(kind = 8)   , INTENT(INOUT)    :: time(:)
        INTEGER(kind = 4), INTENT(INOUT)    :: flag(:)

        INTEGER(kind = 8)                   :: i

        flag = 0
        detect_missing_loop: do i = 1, size(time) - 1, 1
            if (time(i + 1) - time(i) > 1.0_wp) then
                flag(i + 1) = 1
            end if
        end do detect_missing_loop

        write(*, *) 'Detecting missing point completed...'
    
    end subroutine detect_missing

    subroutine output_data2file(me, o_file_unit)
        implicit none
        type(gps), INTENT(IN   )            :: me(:)

        integer(kind = 4)                   :: o_file_unit

        integer(kind = 8)                   :: i

        output_loop: do i = 1, size(me), 1
            write(o_file_unit, '(f15.1, 1x, 10f24.12, 1x, i1)') s(i)%time, s(i)%rcv_time, s(i)%pos, s(i)%vel, s(i)%acc, s(i)%vac_qualflg
        end do output_loop

    end subroutine output_data2file

    subroutine creat_clk_file(me)
        implicit none

        TYPE(gps), intent(in   )            :: me(:)
        real(wp), allocatable               :: gps_time(:)

        INTEGER(ip)                         :: err

        allocate(gps_time(size(me)), stat=err)
        if (err /= 0) print *, "gps_time: Allocation request denied"

        gps_time = me%time - time(1)

        if (allocated(gps_time)) deallocate(gps_time, stat=err)
        if (err /= 0) print *, "gps_time: Deallocation request denied"

    end subroutine creat_clk_file
    
end module taiji_01_preprocessing

program main
    use taiji_01_preprocessing
    use num_kinds

    implicit none

    integer(kind = 4)                       :: num_input_files
    integer(kind = 4)                       :: err1, err2, ios, err3, ios1, err
    integer(kind = 4)                       :: nheader
    integer(kind = 8)                       :: index, i, j

    character(len = 300), ALLOCATABLE       :: file_header(:)
    character(len = 100)                    :: date, time, zone

    real(kind = 8)                          :: temp
    real(kind = 8)                          :: start_epoch, stop_epoch

    call CPU_TIME(start_epoch)

    num_input_files = 1
    nheader = 71
    log_f%name = '..//log.txt'
    log_f%unit = 89
    
    ! open the log file
    open(unit=log_f%unit, file=log_f%name, iostat=ios, &
         status="unknown", action="write", position="append")
    if (ios /= 0) stop "Error opening log file"

    write(log_f%unit, *) '--------------------------------------------------------------------'
    call DATE_AND_TIME(date, time, zone)
    write(log_f%unit, *) 'Date: ', date
    write(log_f%unit, *) 'Time: ', time
    write(*, *) 'Programme launching...'
    write(log_f%unit, *) 'Programme launching'

    ! allocate the list to contain the file header
    allocate(file_header(nheader), stat=err1)
    if (err1 /= 0) print *, "file_header: Allocation request denied"
    
    ! allocate the list to contain the information of input files
    allocate(i_f(num_input_files), stat=err2)
    if (err2 /= 0) print *, "i_f: Allocation request denied"

    allocate(f_f(num_input_files), stat=err)
    if (err /= 0) print *, "f_f: Allocation request denied"

    allocate(o_f(num_input_files), stat=err)
    if (err /= 0) print *, "o_f: Allocation request denied"

    i_f%name = [ &
        '..//output//taiji-01-0860-earth-fixed-system-2019-12.txt' &
    ]

    f_f%name = [ &
        '..//output//taiji-01-0860-position-velocity-flag-2019-12.txt' &
    ]
    
    o_f%name = [ &
        '..//output//taiji-01-0860-earth-fixed-system-2019-12-after-fortran.txt'&
    ]

    month_file_loop: do index = 1, num_input_files, 1
        i_f(index)%unit = 34
        f_f(index)%unit = 76
        o_f(index)%unit = 29
        
        open(unit=i_f(index)%unit, file=i_f(index)%name, iostat=ios, status="old", action="read")
        if (ios /= 0) stop "Error opening file unit i_f(index)%unit"

        open(unit=f_f(index)%unit, file=f_f(index)%name, iostat=ios, status="old", action="read")
        if (ios /= 0) stop "Error opening file unit f_f(index)%unit"

        open(unit=o_f(index)%unit, file=o_f(index)%name, iostat=ios, status="unknown", action="write")
        if (ios /= 0) stop "Error opening file unit o_f(index)%unit"
        
        ! get the number of rowsof the input files
        i_f(index)%nrow = get_file_n(i_f(index)%unit)
        f_f(index)%nrow = get_file_n(f_f(index)%unit)
        if (i_f(index)%nrow - nheader /= f_f(index)%nrow - 1) stop "The nrows of gps file and increment file do not conform"

        ! read headers
        read_header: do i = 1, nheader, 1
            read(i_f(index)%unit, '(a)') file_header(i)
            write(o_f(index)%unit, '(a)') trim(file_header(i))
        end do read_header

        read(f_f(index)%unit, '(a)') temp

        ! allocate the gps data matrix of this month
        allocate(s(i_f(index)%nrow - nheader), stat=err3)
        if (err3 /= 0) print *, "s: Allocation request denied"

        read_data_gps: do i = 1, i_f(index)%nrow - nheader, 1
            read(i_f(index)%unit, *) s(i)%time, s(i)%rcv_time, s(i)%pos, s(i)%vel
        end do read_data_gps

        read_data_increment: do i = 1, f_f(index)%nrow - 1, 1
            read(f_f(index)%unit, *) s(i)%pos_qualflg, s(i)%vel_qualflg, s(i)%pos_incr, s(i)%vel_incr
        end do read_data_increment

        call gps_deglitch(s)

        call diff9(s%time, s%vel(1), s%acc(1), 1.0_wp)
        call diff9(s%time, s%vel(2), s%acc(2), 1.0_wp)
        call diff9(s%time, s%vel(3), s%acc(3), 1.0_wp)

        call detect_missing(s%time, s%vac_qualflg)

        call output_data2file(s, o_f(index)%unit)

        call creat_clk_file(s)
        
        close(unit=i_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit i_f(index)%unit"

        close(unit=f_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit f_f(index)%unit"

        close(unit=o_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit f_f(index)%unit"
        
        if (allocated(s)) deallocate(s, stat=err3)
        if (err3 /= 0) print *, "s: Deallocation request denied"
    end do month_file_loop

    if (allocated(file_header)) deallocate(file_header, stat=err1)
    if (err1 /= 0) print *, "file_header: Deallocation request denied"

    if (allocated(i_f)) deallocate(i_f, stat=err2)
    if (err2 /= 0) print *, "i_f: Deallocation request denied"

    if (allocated(f_f)) deallocate(f_f, stat=err)
    if (err /= 0) print *, "f_f: Deallocation request denied"

    if (allocated(o_f)) deallocate(o_f, stat=err)
    if (err /= 0) print *, "o_f: Deallocation request denied"

    call CPU_TIME(stop_epoch)

    write(log_f%unit, *) 'Exited with code=0 in ', stop_epoch - start_epoch, ' seconds.'
    write(*, *) 'Exited with code=0 in ', stop_epoch - start_epoch, ' seconds.'

    close(unit=log_f%unit, iostat=ios)
    if (ios /= 0) stop "Error closing file unit log file"
end program main
