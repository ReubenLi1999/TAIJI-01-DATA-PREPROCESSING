module taiji_01_preprocessing
    use num_kinds
    use bspline_module

    implicit none

    type, public:: spacecraft
        real(kind = 8)                      :: pos_e(3)                 !! position vector in earth-fixed frame
        real(kind = 8)                      :: pos_i(3)                 !! position vector in the inertial frame
        real(kind = 8)                      :: vel_e(3)                 !! velocity vector in earth-fixed frame
        real(kind = 8)                      :: acc_e(3)                 !! accelaration vector in earth-fixed frame
        real(kind = 8)                      :: pos_e_incr(3)            !! the possible jump of position vector in earth-fixed frame
        real(kind = 8)                      :: vel_e_incr(3)            !! the possible jump of velocity vector in earth-fixed frame
        real(kind = 8)                      :: vel_e_diffed(3)          !! the differentiation of velocity vector in earth-fixed frame
        real(kind = 8)                      :: rvel_i(3)                !! the relative velocity to the atmosphere in the gcrs
        real(kind = 8)                      :: rvel_s(3)                !! the relative velocity to the atmosphere in the srf
        real(kind = 8)                      :: i2s_eul(3)               !! the euler angles representing the attitude from gcrs to srf
        real(kind = 8)                      :: i2s_m(3, 3)              !! the direction cosine matrix of the corresponding i2s_eul
        real(kind = 8)                      :: s2i_m(3, 3)              !! the direction cosine matrix from srf to gcrs
        real(kind = 8)                      :: t2s_eul(3)               !! the euler angles representing the attitude from trajectory frame to srf
        real(kind = 8)                      :: t2s_m(3, 3)              !! the direction cosine matrix of the corresponding t2s_eul
        real(kind = 8)                      :: time                     !! the corresponding gps time
        real(kind = 8)                      :: rcv_time                 !! the correspongding receiving time or satellite time
        real(kind = 8)                      :: ambient_density          !! the ambient density of this satellite
        real(kind = 8)                      :: area(6)                  !! the surface area of the satellite
        real(kind = 8)                      :: theta_v2n(6)             !! the angle of the surface normal vector and relative velocity in radian
        real(kind = 8)                      :: air_drag_i(3)            !! the air drag force in the gcrs
        real(kind = 8)                      :: air_drag_s(3)            !! the air drag force in the srf
        real(kind = 8)                      :: eul(3)                   !! the ready euler angles
        real(kind = 8)                      :: m(3, 3)                  !! the ready direction cosine matrix
        real(kind = 8)                      :: surface_vector_s(6, 3)   !! the surface normal vectors of the satellite in the srf
        real(kind = 8)                      :: surface_vector_i(6, 3)   !! the surface normal vectors of the satellite in the gcrs
        real(kind = 8)                      :: mass                     !! the mass of the whole satellite
        real(kind = 8)                      :: right_ascension          !! the right ascension of the Sun on the point
        real(kind = 8)                      :: negative_declination     !! the negative decliantion of the Sun on the point
        real(kind = 8)                      :: pos_n(3)                 !! the new position vector in the new-formed frame
        real(kind = 8)                      :: sp_c(3)                  !! the absoption, diffused and spectral coefficients
        real(kind = 8)                      :: uvector_s2sun_i(3)         !! the unit vector pointing from the satellite to the Sun
        real(kind = 8)                      :: phi_sur2sun(6)           !! the angle between surface normal and the direction to the Sun
        real(kind = 8)                      :: jdtt                     !! the julian day form of the utc for this point
        real(kind = 8)                      :: posvel_sun_i(6)          !! the position and velocity vector of the Sun in the centre of the Earth in the inertial frame
        real(kind = 8)                      :: solar_pressure_i(3)      !! the solar pressure accelaration in the inertial frame on this point
        real(kind = 8)                      :: solar_pressure_s(3)      !! the solar pressure accelaration in the srf
        real(kind = 8)                      :: vector1(3), vector2(3)
        real(kind = 8)                      :: angle
        integer(kind = 4)                   :: pos_e_qualflg            !! the quality flag for pos_eition
        integer(kind = 4)                   :: vel_e_qualflg            !! the quality flag for velocity
        integer(kind = 4)                   :: vac_qualflg              !! the flag that signals the vacuum before the epoch
        integer(kind = 4)                   :: month                    !! the month for the epoch of the spacecraft
        integer(kind = 4)                   :: dayofyear                !! the day of this year
        integer(kind = 4)                   :: secofday                 !! the seconds of this day
        integer(kind = 4)                   :: lat                      !! the latitude in the geodetic frame
        integer(kind = 4)                   :: lon                      !! the longitude in the geodetic frame
        integer(kind = 4)                   :: alt                      !! the altitude in the geodetic frame
        integer(kind = 4)                   :: year                     !! the year of this epoch
        integer(kind = 4)                   :: flag_solar_pressure      !! 1 for not in the shadow of the Earth, 0 for in the shadow of the Earth
    contains
        procedure                           :: eul2m                      => eul2m
        procedure                           :: angle2vectors              => angle2vectors
        procedure                           :: air_drag                   => air_drag
        procedure                           :: utct2jdtt                  => utct2jdtt
        procedure                           :: solar_pressure             => solar_pressure
        procedure                           :: date2radecl                => date2RADECL
        procedure                           :: check_shadow               => check_shadow

    end type spacecraft

    type io_file
        character(len = 100)                :: name
        integer(kind = 4)                   :: unit
        integer(kind = 8)                   :: nrow
        integer(kind = 4)                   :: nheader
        character(len = 300), ALLOCATABLE   :: file_header(:)
    contains
        procedure                           :: write_air_drag_file        => write_air_drag_file
    end type io_file

    class(spacecraft)    , allocatable      :: s(:)
    class(io_file)       , allocatable      :: i_f(:) !! the input file
    class(io_file)       , allocatable      :: o_f(:) !! the output file
    class(io_file)       , ALLOCATABLE      :: f_f(:) !! the position and velocity flag file
    class(io_file)       , ALLOCATABLE      :: c_f(:) !! clk file
    class(io_file)       , ALLOCATABLE      :: d_f(:) !! the air density file
    class(io_file)       , ALLOCATABLE      :: z_f(:) !! the attitude file
    class(io_file)       , ALLOCATABLE      :: v_f(:) !! the relative velocity file included in the positions and velocities file in gcrs
    class(io_file)       , ALLOCATABLE      :: a_f(:) !! the air drag accelaration file in the gcrs
    class(io_file)       , ALLOCATABLE      :: s_f(:) !! the file containing the vector pointing from the Earth to the Sun in the gcrs
    type(io_file)                           :: log_f  !! the log file

    real(wp), PARAMETER                     :: pi = acos(-1.0_wp)
    real(wp), PARAMETER                     :: radius_earth = 6371.0_wp     !! the radius of the earth in km
    real(wp), PARAMETER                     :: fe = 1370.0_wp               !! the solar flux, unit=w/m^2
    real(wp), PARAMETER                     :: c = 299792458.0_wp        !! the speed of light
    real(wp), PARAMETER                     :: au_unit = 149597871.0_wp     !! Astronomical unit
    
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

        type(spacecraft), intent(inout)     :: me(:)

        INTEGER(kind = 8)                   :: index, i, j

        deglitch_loop: do i = 1, size(me), 1
            if (me(i)%pos_e_qualflg == 1) then
                me(i: size(me))%pos_e(1) = me(i: size(me))%pos_e(1) - me(i)%pos_e_incr(1)
                me(i: size(me))%pos_e(2) = me(i: size(me))%pos_e(2) - me(i)%pos_e_incr(2)
                me(i: size(me))%pos_e(3) = me(i: size(me))%pos_e(3) - me(i)%pos_e_incr(3)
                me(i: size(me))%vel_e(1) = me(i: size(me))%vel_e(1) - me(i)%vel_e_incr(1)
                me(i: size(me))%vel_e(2) = me(i: size(me))%vel_e(2) - me(i)%vel_e_incr(2)
                me(i: size(me))%vel_e(3) = me(i: size(me))%vel_e(3) - me(i)%vel_e_incr(3)
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
        type(spacecraft), INTENT(IN   )     :: me(:)

        integer(kind = 4)                   :: o_file_unit

        integer(kind = 8)                   :: i

        write(o_file_unit, '(a)') 'gps_time, rcv_time, pos_e, vel, acc, vac_qualflg, vel_diffed, pj_qualflg, vp_qualflg'

        output_loop: do i = 1, size(me), 1
            write(o_file_unit, '(f15.1, 1x, 10f24.12, 1x, i1, 1x, 3f24.12, 1x, i1, 1x, i1)') s(i)%time, s(i)%rcv_time, &
                                    s(i)%pos_e, s(i)%vel_e, s(i)%acc_e, s(i)%vac_qualflg, s(i)%vel_e_diffed, s(i)%pos_e_qualflg, s(i)%vel_e_qualflg
        end do output_loop

    end subroutine output_data2file

    subroutine creat_clk_file(me, o_file_unit)
        implicit none

        TYPE(spacecraft)  , intent(in   )     :: me(:)
        INTEGER(ip)       , INTENT(IN   )     :: o_file_unit
        real(wp), allocatable                 :: gps_time(:)
        real(wp), ALLOCATABLE                 :: clk_drift(:)

        INTEGER(ip)                           :: err, i

        CHARACTER(len = 80)                   :: date, temp_c

        allocate(gps_time(size(me)), stat=err)
        if (err /= 0) print *, "gps_time: Allocation request denied"
        allocate(clk_drift(size(me)), stat=err)
        if (err /= 0) print *, "clk_drift: Allocation request denied"

        !! obtain the date
        call DATE_AND_TIME(date, temp_c, temp_c)

        gps_time = me%time - me(1)%time
        clk_drift = gps_time - me%rcv_time + floor(me(1)%rcv_time)
        write(o_file_unit, '(a)') 'header:'
        write(o_file_unit, '(4x, a, i10)') 'dimensions: ', size(me)
        write(o_file_unit, '(4x, a)') 'global_attributes: '
        write(o_file_unit, '(8x, a)') 'acknowledgement: TAIJI-01 is the first grativational prospecting test satellite of China, one of whose core loads is inertial sensor, therefore, the GPS data of this satellite can be used in the inversion of coefficients of earth gravity field.'
        WRITE(o_file_unit, '(8x, a)') 'creator_institution: CAS/IM'
        write(o_file_unit, '(8x, a)') 'creator_name: TAIJI-01 data product system'
        write(o_file_unit, '(8x, a)') 'creator_type: group'
        write(o_file_unit, '(8x, a, a)') 'date_created: ', date
        write(o_file_unit, '(8x, a)') 'institution: CAU'
        write(o_file_unit, '(8x, a)') 'instrument: GPS and local clock'
        write(o_file_unit, '(8x, a)') 'keywords: TAIJI-01, gravity field'
        WRITE(o_file_unit, '(8x, a)') 'processing_level: 1E'
        write(o_file_unit, '(8x, a)') 'data_product: 0860'
        write(o_file_unit, '(8x, a)') 'programme: Taiji Programme'
        WRITE(o_file_unit, '(8x, a)') 'publisher_institution: CAS/ISSI'
        write(o_file_unit, '(8x, a)') 'source: Earth-Fixed Frame trajectories for Taiji-01'
        write(o_file_unit, '(8x, a)') 'summary: 1-Hz drift between the GPS time and satellite time'
        write(o_file_unit, '(8x, a, f15.7)') 'time_coverage_start: ', me(1)%time
        write(o_file_unit, '(8x, a, f15.7)') 'time_coverage_stop: ', me(size(me))%time
        write(o_file_unit, '(8x, a)') 'title: Taiji-01 Clock Offset Data'
        write(o_file_unit, '(4x, a)') 'varaiables: '
        write(o_file_unit, '(8x, a)') '- self_rcv_time: '
        write(o_file_unit, '(12x, a)') 'comment: 1st column'
        write(o_file_unit, '(12x, a)') 'long_name: Continuous seconds past the beginning epoch of this month'
        write(o_file_unit, '(12x, a)') 'unit: second'
        write(o_file_unit, '(8x, a)') '- time_drift: '
        write(o_file_unit, '(12x, a)') 'comment: 2nd column'
        write(o_file_unit, '(12x, a)') 'long_name: Clock Offset'
        write(o_file_unit, '(12x, a)') 'unit: second'
        write(o_file_unit, '(a)') '# End of header'
        write(o_file_unit, '(a, 4x, a)') 'rcv_time', 'time_drift'

        output_clk_loop: do i = 1, size(me), 1
            write(o_file_unit, '(f20.10, 4x, f20.10)') me(i)%rcv_time, clk_drift(i)
        end do output_clk_loop

        if (allocated(gps_time)) deallocate(gps_time, stat=err)
        if (err /= 0) print *, "gps_time: Deallocation request denied"
        if (allocated(clk_drift)) deallocate(clk_drift, stat=err)
        if (err /= 0) print *, "clk_drift: Deallocation request denied"

    end subroutine creat_clk_file

    subroutine eul2m(me, d_or_r, seq)
        implicit none

        class(spacecraft) , INTENT(INout)   :: me
        CHARACTER(len = 1), INTENT(IN   )   :: d_or_r

        real(wp), PARAMETER                 :: d2r_factor = acos(-1.0_wp) / 180.0_wp

        real(wp)                            :: r(3, 3, 3)

        INTEGER(ip)                         :: i
        INTEGER(IP)                         :: nrows
        INTEGER(ip)                         :: seq(3)

        nrows = size(me%eul) / 3_ip
        me%m = 0
        !! degree to radian
        if (d_or_r == 'd') then
            me%eul = me%eul * d2r_factor
        else if ((d_or_r /= 'd') .or. (d_or_r /= 'r')) then
            stop "Error with bad degree/radian flag for computing rotation matrix from euler angles"
        end if

        ! the first rotation matrix
        r(1, 1, 1) =  1.0_wp
        r(1, 1, 2) =  0.0_wp
        r(1, 1, 3) =  0.0_wp
        r(1, 2, 1) =  0.0_wp
        r(1, 2, 2) =  cos(me%eul(1))
        r(1, 2, 3) =  sin(me%eul(1))
        r(1, 3, 1) =  0.0_wp
        r(1, 3, 2) = -sin(me%eul(1))
        r(1, 3, 3) =  cos(me%eul(1))

        ! the second rotation matrix
        r(2, 1, 1) =  cos(me%eul(2))
        r(2, 1, 2) =  0.0_wp
        r(2, 1, 3) = -sin(me%eul(2))
        r(2, 2, 1) =  0.0_wp
        r(2, 2, 2) =  1.0_wp
        r(2, 2, 3) =  0.0_wp
        r(2, 3, 1) =  sin(me%eul(2))
        r(2, 3, 2) =  0.0_wp
        r(2, 3, 3) =  cos(me%eul(2))

        ! the third rotation matrix
        r(3, 1, 1) =  cos(me%eul(3))
        r(3, 1, 2) =  sin(me%eul(3))
        r(3, 1, 3) =  0.0_wp
        r(3, 2, 1) = -sin(me%eul(3))
        r(3, 2, 2) =  cos(me%eul(3))
        r(3, 2, 3) =  0.0_wp
        r(3, 3, 1) =  0.0_wp
        r(3, 3, 2) =  0.0_wp
        r(3, 3, 3) =  1.0_wp

        ! creat the rotation matrix following the sequence
        me%m = matmul(r(seq(3), :, :), matmul(r(seq(2), :, :), r(seq(1), :, :)))

    end subroutine eul2m

    subroutine angle2vectors(me)
        implicit none

        class(spacecraft), INTENT(INOUT)            :: me

        me%angle = acos(dot_product(me%vector1, me%vector2) / (norm2(me%vector1) * norm2(me%vector2)))
    end subroutine angle2vectors

    subroutine air_drag(me, flag, j)
        implicit none
        class(spacecraft) , INTENT(INOUT)                :: me

        CHARACTER(len = 1), INTENT(IN   )                :: flag

        integer(kind = 8) , INTENT(IN   )                :: j

        integer(kind = 4)                                :: k

        if (flag == 'i') then 
            if (me%theta_v2n(j) < (pi / 2.0)) then
                    k = ceiling(j / 2.0_wp)
                    me%air_drag_i = 1.0_wp / 2.0_wp * me%ambient_density * & !! ambient air density
                                    me%area(j) * & !! area of the surface
                                    norm2(me%rvel_i)**2 * & !! magnitude of the velocty relative to the atmosphere
                                    (-2.2_wp) * cos(me%theta_v2n(j)) * me%rvel_i / norm2(me%rvel_i) / & !! ballistic coefficient vector
                                    me%mass & !! the mass of the satellite
                                    + me%air_drag_i
            end if
        else if (flag == 's') then
            if (me%theta_v2n(j) < (pi / 2.0)) then
                    k = ceiling(j / 2.0_wp)
                    me%air_drag_s = 1.0_wp / 2.0_wp * me%ambient_density * & !! ambient air density
                                    me%area(j) * & !! area of the surface
                                    norm2(me%rvel_s)**2 * & !! magnitude of the velocty relative to the atmosphere
                                    (-2.2_wp) * cos(me%theta_v2n(j)) * me%rvel_s / norm2(me%rvel_s) / & !! ballistic coefficient vector
                                    me%mass & !! the mass of the satellite
                                    + me%air_drag_s
            end if
        end if
    end subroutine air_drag

    subroutine write_air_drag_file(me)
        implicit none
        CLASS(io_file), INTENT(INOUT)               :: me

        INTEGER(ip)                                 :: err, i

        CHARACTER(len = 80)                         :: date, temp_c

        !! obtain the date
        call DATE_AND_TIME(date, temp_c, temp_c)

        write(me%unit, '(a)') 'header:'
        write(me%unit, '(4x, a)') 'global_attributes: '
        write(me%unit, '(8x, a)') 'acknowledgement: TAIJI-01 is the first grativational prospecting test satellite of China, one of whose core loads is inertial sensor, therefore, the GPS data of this satellite can be used in the inversion of coefficients of earth gravity field.'
        WRITE(me%unit, '(8x, a)') 'creator_institution: CAS/IM'
        write(me%unit, '(8x, a)') 'creator_name: TAIJI-01 data product system'
        write(me%unit, '(8x, a)') 'creator_type: group'
        write(me%unit, '(8x, a, a)') 'date_created: ', date
        write(me%unit, '(8x, a)') 'institution: CAU'
        write(me%unit, '(8x, a)') 'instrument: GPS and local clock'
        write(me%unit, '(8x, a)') 'keywords: TAIJI-01, gravity field'
        WRITE(me%unit, '(8x, a)') 'processing_level: 1E'
        write(me%unit, '(8x, a)') 'data_product: 0222'
        write(me%unit, '(8x, a)') 'programme: Taiji Programme'
        WRITE(me%unit, '(8x, a)') 'publisher_institution: CAS/ISSI'
        write(me%unit, '(8x, a)') 'source: Earth-Fixed Frame trajectories for Taiji-01'
        write(me%unit, '(8x, a)') 'summary: 1-Hz air drag accelaration in the inertial frame'
        write(me%unit, '(8x, a)') 'title: Taiji-01 Clock Offset Data'
        write(me%unit, '(4x, a)') 'varaiables: '
        write(me%unit, '(8x, a)') '- gps_time: '
        write(me%unit, '(12x, a)') 'comment: 1st column'
        write(me%unit, '(12x, a)') 'long_name: Continuous seconds past the beginning epoch of this month'
        write(me%unit, '(12x, a)') 'unit: second'
        write(me%unit, '(8x, a)') '- air_drag_acc_x: '
        write(me%unit, '(12x, a)') 'comment: 2nd column'
        write(me%unit, '(12x, a)') 'long_name: the air drag accelaration along the x axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(8x, a)') '- air_drag_acc_y: '
        write(me%unit, '(12x, a)') 'comment: 3rd column'
        write(me%unit, '(12x, a)') 'long_name: the air drag accelaration along the y axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(8x, a)') '- air_drag_acc_z: '
        write(me%unit, '(12x, a)') 'comment: 4th column'
        write(me%unit, '(12x, a)') 'long_name: the air drag accelaration along the z axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(8x, a)') '- solar_pressure_acc_x: '
        write(me%unit, '(12x, a)') 'comment: 5th column'
        write(me%unit, '(12x, a)') 'long_name: the solar pressure accelaration along the x axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(8x, a)') '- solar_pressure_acc_y: '
        write(me%unit, '(12x, a)') 'comment: 6th column'
        write(me%unit, '(12x, a)') 'long_name: the solar pressure accelaration along the y axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(8x, a)') '- solar_pressure_acc_z: '
        write(me%unit, '(12x, a)') 'comment: 7th column'
        write(me%unit, '(12x, a)') 'long_name: the solar pressure accelaration along the z axis of the gcrs'
        write(me%unit, '(12x, a)') 'unit: km/s^2'
        write(me%unit, '(a)') '# End of header'

    end subroutine write_air_drag_file

    subroutine utct2jdtt(me)
        class(spacecraft), INTENT(INOUT)                    :: me

        me%jdtt = (me%time + 32.184_wp + 37.0_wp) / 86400.0_wp + 2458118.833333333_wp

    end subroutine utct2jdtt

    subroutine solar_pressure(me, flag, j)
        class(spacecraft) , INTENT(INout)                   :: me

        CHARACTER(len = 1), INTENT(IN   )                   :: flag

        integer(kind = 8) , INTENT(IN   )                   :: j

        integer(kind = 4)                                   :: k

        if (me%flag_solar_pressure == 1_ip) then
            if (flag == 'i') then 
                if (me%phi_sur2sun(j) < (pi / 2.0)) then
                        k = ceiling(j / 2.0_wp)
                        me%solar_pressure_i = me%solar_pressure_i + &       !! the accumulation term
                                            (-1) * cos(me%phi_sur2sun(j)) * &!! the cosine of the angle between the surface normal and the direction to the Sun
                                            me%area(j) *        &         !! the area of the surface
                                            fe / c *              &         !! the constant in the equation, solar flux and speed of light
                                            ((1.0_wp - me%sp_c(1)) * me%uvector_s2sun_i + &
                                            2.0_wp * me%surface_vector_i(j, :) * &
                                            (me%sp_c(1) * cos(me%phi_sur2sun(j)) + 1.0_wp / 3.0_wp * me%sp_c(2))) / &
                                            me%mass
                end if
            end if
        else
            me%solar_pressure_i = 0.0_wp
        end if

    end subroutine solar_pressure

    subroutine date2RADECL(me)
        implicit none
        class(spacecraft), INTENT(INOUT)                    :: me

        real(kind = 8) jdstart, t
        real(kind = 8) L0, M, Ci, THETA, OMIGA, lamda, year, THETA2000, eph, s1, s2, s3, s4

        jdstart = 2451545.0_wp
        t = (me%jdtt - jdstart) / 36525.0_wp
        
        L0=280.46645_wp + 36000.76983_wp * t + 0.0003032_wp * t**2
        
        M = 357.52910_wp + 35999.05030_wp * t - 0.0001559_wp *t**2 - 0.00000048_wp * t**3
        
        Ci = (1.914600_wp - 0.004817_wp * t - 0.000014_wp * t**2) * sind(M) + &
                (0.019993_wp - 0.000101_wp * t) * sind(2.0_wp * M) + 0.000290_wp * sind(3.0_wp * M)
        
        THETA = L0 + Ci
        
        OMIGA = 125.04_wp - 1934.136_wp * t
        lamda = THETA - 0.00569_wp - 0.00478_wp * sind(OMIGA)
        
        year = floor((me%jdtt - jdstart) / 365.25_wp) + 2000.0_wp
        THETA2000 = lamda - 0.01397_wp * (year - 2000.0_wp)
        
        eph = 23.439291_wp - 0.013004_wp * t - 0.00059_wp * t / 3600.0_wp + 0.001813_wp*t**3 / 3600.0_wp + 0.00256_wp * cosd(OMIGA)
        
        s1 = cosd(eph) * sind(THETA2000)
        s2 = cosd(THETA2000)
        s3 = sind(eph)
        s4 = sind(THETA2000)
        
        me%right_ascension = atan2d(cosd(eph) * sind(THETA2000), cosd(THETA2000))
        
        me%negative_declination = asind(sind(eph) * sind(THETA2000))  
    end subroutine date2RADECL

    subroutine check_shadow(me)
        class(spacecraft), INTENT(INOUT)                        :: me
        real(wp)                                                :: r(3, 3, 3)

        ! the second rotation matrix
        r(2, 1, 1) =  cosd(-1.0_wp * me%negative_declination)
        r(2, 1, 2) =  0.0_wp
        r(2, 1, 3) = -sind(-1.0_wp * me%negative_declination)
        r(2, 2, 1) =  0.0_wp
        r(2, 2, 2) =  1.0_wp
        r(2, 2, 3) =  0.0_wp
        r(2, 3, 1) =  sind(-1.0_wp * me%negative_declination)
        r(2, 3, 2) =  0.0_wp
        r(2, 3, 3) =  cosd(-1.0_wp * me%negative_declination)

        ! the third rotation matrix
        r(3, 1, 1) =  cosd(me%right_ascension)
        r(3, 1, 2) =  sind(me%right_ascension)
        r(3, 1, 3) =  0.0_wp
        r(3, 2, 1) = -sind(me%right_ascension)
        r(3, 2, 2) =  cosd(me%right_ascension)
        r(3, 2, 3) =  0.0_wp
        r(3, 3, 1) =  0.0_wp
        r(3, 3, 2) =  0.0_wp
        r(3, 3, 3) =  1.0_wp

        me%pos_n = matmul(r(2, :, :), matmul(r(3, :, :), me%pos_i))
        
        me%flag_solar_pressure = 1_ip
        if (me%pos_n(1) < 0.0_wp) then
            if (me%pos_n(2)**2 + me%pos_n(3)**2 <= radius_earth**2) then
                me%flag_solar_pressure = 0_ip
                print *, 1
            end if
        end if
    end subroutine check_shadow
    
end module taiji_01_preprocessing

program main
    use taiji_01_preprocessing
    use num_kinds
    use Vars_wlm

    implicit none

    integer(kind = 4)                       :: num_input_files
    integer(kind = 4)                       :: err1, err2, ios, err3, ios1, err
    integer(kind = 8)                       :: index, i, j, k
    integer(kind = 4)                       :: count

    character(len = 100)                    :: date, time, zone

    real(kind = 8)                          :: temp
    real(kind = 8)                          :: temp_i2s_eul(3), temp_t2s_eul(3)
    real(kind = 8)                          :: start_epoch, stop_epoch

    INTEGER(kind = 4)                       :: air_drag_or_not
    INTEGER(kind = 4)                       :: solar_pressure_or_not

    call CPU_TIME(start_epoch)

    num_input_files = 1
    log_f%name = '..//log.txt'
    log_f%unit = 89
    
    ! open the log file
    open(unit=log_f%unit, file=log_f%name, iostat=ios, &
         status="unknown", action="write", position="append")
    if (ios /= 0) stop "Error opening log file"

    write(log_f%unit, *) '-------------------------------------------------------------------------'
    call DATE_AND_TIME(date, time, zone)
    write(log_f%unit, *) 'Date: ', date
    write(log_f%unit, *) 'Time: ', time
    write(*, *) 'Programme launching...'
    write(log_f%unit, *) 'Programme launching'
    
    !----------------------------------------------------------------------------------------------!
    ! allocate the list to contain the information of input files
    allocate(i_f(num_input_files), stat=err) !! the input positions and velocites file in itrs
    if (err /= 0) print *, "i_f: Allocation request denied"

    allocate(f_f(num_input_files), stat=err) !! the positions and velocities flag file
    if (err /= 0) print *, "f_f: Allocation request denied"

    allocate(o_f(num_input_files), stat=err) !! the output file in itrs
    if (err /= 0) print *, "o_f: Allocation request denied"

    allocate(c_f(num_input_files), stat=err) !! the clk file from satellite time to gps time
    if (err /= 0) print *, "c_f: Allocation request denied"

    ALLOCATE(d_f(num_input_files), stat=err) !! the ambient air density file
    if (err /= 0) print *, "d_f: Allocation request denied"

    ALLOCATE(z_f(num_input_files), stat=err) !! the attitude file for tf2srf and gcrs2srf
    if (err /= 0) print *, "z_f: Allocation request denied"

    ALLOCATE(v_f(num_input_files), stat=err) !! the input relative velocities in gcrs
    if (err /= 0) print *, "v_f: Allocation request denied"

    ALLOCATE(a_f(num_input_files), stat=err) !! the output air drag exerting on the satellite in the gcrs
    if (err /= 0) print *, "a_f: Allocation request denied"

    ALLOCATE(s_f(num_input_files), stat=err) !! the file containing the vector pointing the Sun from the Earth in the gcrs
    if (err /= 0) print *, "s_f: Allocation request denied"
    !----------------------------------------------------------------------------------------------！

    !----------------------------------------------------------------------------------------------!
    i_f%name = [ &
        '..//output//09//taiji-01-0860-earth-fixed-system-2019-09.txt' &
    ]

    f_f%name = [ &
        '..//output//09//taiji-01-0860-position-velocity-flag-2019-09.txt' &
    ]
    
    o_f%name = [ &
        '..//output//09//taiji-01-0860-earth-fixed-system-2019-09-after-fortran.txt'&
    ]

    c_f%name = [ &
        '..//input//09//KX09_SAT_1D_ENG_0111_20190900T000000_01_04.txt'&
    ]

    d_f%name = [ &
        '..//output//09//taiji-01-0860-air-density-2019-09.txt'&
    ]

    z_f%name = [ &
        '..//output//09//taiji-01-0811-attitude-2019-09.txt' &
    ]

    v_f%name = [ &
        '..//output//09//taiji-01-0866-gcrs-2019-09.txt' &
    ]

    a_f%name = [ &
        '..//output//09//taiji-01-0222-non-gravitational-gcrs-2019-09.txt' &
    ]

    s_f%name = [ &
        '..//output//09//taiji-01-0333-earth2sun-2019-09.txt' &
    ]
    !----------------------------------------------------------------------------------------------!

    month_file_loop: do index = 1, num_input_files, 1
        i_f(index)%unit = 34
        f_f(index)%unit = 76
        o_f(index)%unit = 29
        c_f(index)%unit = 59
        d_f(index)%unit = 39
        z_f(index)%unit = 98
        v_f(index)%unit = 68
        a_f(index)%unit = 49
        s_f(index)%unit = 123

        i_f(index)%nheader = 91
        f_f(index)%nheader = 1
        d_f(index)%nheader = 1
        z_f(index)%nheader = 44
        v_f(index)%nheader = 50

        !------------------------------------------------------------------------------------------!
        ! allocate the list to contain the file header
        allocate(i_f(index)%file_header(i_f(index)%nheader), stat=err)
        if (err /= 0) print *, "file_header: Allocation request denied"

        allocate(z_f(index)%file_header(z_f(index)%nheader), stat=err)
        if (err /= 0) print *, "file_header: Allocation request denied"

        allocate(d_f(index)%file_header(d_f(index)%nheader), stat=err)
        if (err /= 0) print *, "file_header: Allocation request denied"

        allocate(v_f(index)%file_header(v_f(index)%nheader), stat=err)
        if (err /= 0) print *, "file_header: Allocation request denied"
        !------------------------------------------------------------------------------------------!
        
        !------------------------------------------------------------------------------------------!
        open(unit=i_f(index)%unit, file=i_f(index)%name, iostat=ios, status="old", action="read")
        if (ios /= 0) stop "Error opening file unit i_f(index)%unit"

        open(unit=f_f(index)%unit, file=f_f(index)%name, iostat=ios, status="old", action="read")
        if (ios /= 0) stop "Error opening file unit f_f(index)%unit"

        open(unit=o_f(index)%unit, file=o_f(index)%name, iostat=ios, status="unknown", action="write")
        if (ios /= 0) stop "Error opening file unit o_f(index)%unit"

        open(unit=c_f(index)%unit, file=c_f(index)%name, iostat=ios, status="unknown", action="write")
        if (ios /= 0) stop "Error opening file unit c_f(index)%unit"

        open(unit=d_f(index)%unit, file=d_f(index)%name, iostat=ios, status="unknown", action='read')
        if (ios /= 0) stop "Error opening file unit d_f(index)%unit"

        open(unit=z_f(index)%unit, file=z_f(index)%name, iostat=ios, status="unknown", action='read')
        if (ios /= 0) stop "Error opening file unit z_f(index)%unit"

        open(unit=v_f(index)%unit, file=v_f(index)%name, iostat=ios, status='unknown', action='read')
        if (ios /= 0) stop "Error opening file unit v_f(index)%unit"

        open(unit=a_f(index)%unit, file=a_f(index)%name, iostat=ios, status='unknown', action='write')
        if (ios /= 0) stop "Error opening file unit a_f(index)%unit"

        open(unit=s_f(index)%unit, file=s_f(index)%name, iostat=ios, status='unknown', action='write')
        if (ios /= 0) stop "Error opening file unit a_f(index)%unit"
        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        ! get the number of rowsof the input files
        i_f(index)%nrow = get_file_n(i_f(index)%unit)
        f_f(index)%nrow = get_file_n(f_f(index)%unit)
        d_f(index)%nrow = get_file_n(d_f(index)%unit)
        v_f(index)%nrow = get_file_n(v_f(index)%unit)
        z_f(index)%nrow = get_file_n(z_f(index)%unit)

        if (i_f(index)%nrow - i_f(index)%nheader /= f_f(index)%nrow - f_f(index)%nheader) &
        stop "The nrows of gps file and increment file do not conform"

        if (i_f(index)%nrow - i_f(index)%nheader /= z_f(index)%nrow - z_f(index)%nheader) &
        stop "The nrows of gps file and attitude file do not conform"

        if (i_f(index)%nrow - i_f(index)%nheader /= v_f(index)%nrow - v_f(index)%nheader) &
        stop "The nrows of gps file and gcrs file do not conform"

        if (i_f(index)%nrow - i_f(index)%nheader /= d_f(index)%nrow - d_f(index)%nheader) &
        stop "The nrows of gps file and ambient air density file do not conform"
        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        ! read headers
        read_header_i: do i = 1, i_f(index)%nheader - 1, 1
            read(i_f(index)%unit, '(a)') i_f(index)%file_header(i)  
            write(o_f(index)%unit, '(a)') trim(i_f(index)%file_header(i))
        end do read_header_i

        read_header_z: do i = 1, z_f(index)%nheader, 1
            read(z_f(index)%unit, '(a)') z_f(index)%file_header(i)
        end do read_header_z

        read_header_v: do i = 1, v_f(index)%nheader, 1
            read(v_f(index)%unit, '(a)') v_f(index)%file_header(i)
        end do read_header_v

        read(i_f(index)%unit, '(a)') temp
        read(f_f(index)%unit, '(a)') temp
        read(d_f(index)%unit, '(a)') temp
        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        ! allocate the gps data matrix of this month
        allocate(s(i_f(index)%nrow - i_f(index)%nheader), stat=err3)
        if (err3 /= 0) print *, "s: Allocation request denied"
        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        read_data_gps: do i = 1, i_f(index)%nrow - i_f(index)%nheader, 1
            read(i_f(index)%unit, *) s(i)%time, s(i)%rcv_time, s(i)%pos_e, s(i)%vel_e
        end do read_data_gps

        read_data_increment: do i = 1, f_f(index)%nrow - 1, 1
            read(f_f(index)%unit, *) s(i)%pos_e_qualflg, s(i)%vel_e_qualflg, s(i)%pos_e_incr, s(i)%vel_e_incr
        end do read_data_increment

        k = 1
        read_data_attitude: do i = 1, z_f(index)%nrow - z_f(index)%nheader, 1
            read(z_f(index)%unit, *) temp, temp_i2s_eul, temp_t2s_eul
            if (any(s%time == temp)) then
                s(k)%i2s_eul = temp_i2s_eul; s(k)%t2s_eul = temp_t2s_eul; k = k + 1
            end if
        end do read_data_attitude

        read_data_density: do i = 1, d_f(index)%nrow - d_f(index)%nheader, 1
            read(d_f(index)%unit, *) s(i)%ambient_density
        end do read_data_density

        read_data_gcrs: do i = 1, v_f(index)%nrow - v_f(index)%nheader, 1
            read(v_f(index)%unit, *) temp, s(i)%pos_i, s(i)%rvel_i
        end do read_data_gcrs
        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        ! call gps_deglitch(s)

        call diff9(s%time, s%vel_e(1), s%acc_e(1), 1.0_wp)
        call diff9(s%time, s%vel_e(2), s%acc_e(2), 1.0_wp)
        call diff9(s%time, s%vel_e(3), s%acc_e(3), 1.0_wp)

        call diff9(s%time, s%pos_e(1), s%vel_e_diffed(1), 1.0_wp)
        call diff9(s%time, s%pos_e(2), s%vel_e_diffed(2), 1.0_wp)
        call diff9(s%time, s%pos_e(3), s%vel_e_diffed(3), 1.0_wp)

        call detect_missing(s%time, s%vac_qualflg)

        call output_data2file(s, o_f(index)%unit)
        call creat_clk_file(s, c_f(index)%unit)

        call a_f(index)%write_air_drag_file()

        print *, '----------------------------------------------------------------------------'
        print *, 'Air drag computation(1 for yes, 0 for no):'
        print *, '----------------------------------------------------------------------------'
        read(*, *) air_drag_or_not

        print *, '----------------------------------------------------------------------------'
        print *, 'Solar pressure computation(1 for yes, 0 for no):'
        print *, '----------------------------------------------------------------------------'
        read(*, *) solar_pressure_or_not
        
        if (air_drag_or_not == 1 .and. solar_pressure_or_not == 1) then
            air_drag_loop: do i = 1, i_f(index)%nrow - i_f(index)%nheader, 1
                !! mass of the spacecraft
                s(i)%mass = 183.0_wp  !! kg

                !! euler angle to direction cosine matrix
                s(i)%eul = s(i)%i2s_eul
                call s(i)%eul2m('d', [3_ip, 1_ip, 2_ip])
                s(i)%i2s_m = s(i)%m
                s(i)%s2i_m = transpose(s(i)%i2s_m)

                !! surface properties
                s(i)%area = [0.4356_wp, 0.4356_wp, 0.8448_wp, 0.8448_wp, 0.8448_wp, 0.8448_wp] ! m^2
                s(i)%surface_vector_s(1, :) = [ 1.0_wp,  0.0_wp,  0.0_wp]
                s(i)%surface_vector_s(2, :) = [-1.0_wp,  0.0_wp,  0.0_wp]
                s(i)%surface_vector_s(3, :) = [ 0.0_wp,  1.0_wp,  0.0_wp]
                s(i)%surface_vector_s(4, :) = [ 0.0_wp, -1.0_wp,  0.0_wp]
                s(i)%surface_vector_s(5, :) = [ 0.0_wp,  0.0_wp,  1.0_wp]
                s(i)%surface_vector_s(6, :) = [ 0.0_wp,  0.0_wp, -1.0_wp]

                !! utct2julian day
                call s(i)%utct2jdtt()
                !! the right ascension and the negative declination for this moment
                call s(i)%date2radecl()
                !! the position vector of the Sun, 11 for the Sun, 3 for the Earth
                !! notice: the unit of the output is a.u.
                call pleph(s(i)%jdtt, 11_ip, 3_ip, s(i)%posvel_sun_i)
                write(s_f(index)%unit, '(4f20.10)') s(i)%jdtt, s(i)%posvel_sun_i(1: 3)
                s(i)%posvel_sun_i = s(I)%posvel_sun_i * au_unit
                !! the solar pressure spectral diffused absorption coefficients
                s(i)%sp_c = [1.0_wp / 3.0_wp, 1.0_wp / 3.0_wp, 1.0_wp / 3.0_wp]
                !! the unit vector from satellite to the Sun
                s(i)%uvector_s2sun_i = (s(i)%pos_i - s(i)%posvel_sun_i(1: 3)) / norm2(s(i)%pos_i - s(i)%posvel_sun_i(1: 3))
                !! check whether the satellite is in the shadow of the Earth or not
                call s(i)%check_shadow()
                
                print *, matmul(s(i)%i2s_m, s(i)%rvel_i)

                s(i)%air_drag_i = 0.0_wp
                s(i)%solar_pressure_i = 0.0_wp
                loop_surface_vector_trans: do j = 1, 6, 1
                    !! air drag---------------------------------------------------------------------
                    !! transform the surface normal vector from srf to gcrs
                    s(i)%surface_vector_i(j, :) = matmul(s(i)%s2i_m, s(i)%surface_vector_s(j, :))
                    !! the angle between the surface normal vector and the relative velocity in radian
                    s(i)%vector1 = s(i)%rvel_i; s(i)%vector2 = s(i)%surface_vector_i(j, :)
                    call s(i)%angle2vectors()
                    s(i)%theta_v2n(j) = s(i)%angle
                    !! air drag in the inertial frame
                    call s(i)%air_drag('i', j)
                    !!------------------------------------------------------------------------------

                    !! solar pressure---------------------------------------------------------------
                    !! the angle between the surface normal and the direction to the Sun
                    s(i)%vector1 = s(i)%uvector_s2sun_i; s(i)%vector2 = s(i)%surface_vector_i(j, :)
                    call s(i)%angle2vectors()
                    s(i)%phi_sur2sun(j) = s(i)%angle
                    !! solar pressure in the inertial frame
                    call s(i)%solar_pressure('i', j)
                end do loop_surface_vector_trans

                !! air_drag  and solar pressure from gcrs into srf
                s(i)%air_drag_s = matmul(s(i)%i2s_m, s(i)%air_drag_i)
                s(i)%solar_pressure_s = matmul(s(i)%i2s_m, s(i)%solar_pressure_i)
                !! output the air drag and solar preesure accelaration in the gcrs to a_f
                write(a_f(index)%unit, '(f10.1, 4x, 6es23.15)') s(i)%time, s(i)%air_drag_s, s(i)%solar_pressure_s
            end do air_drag_loop
        end if 

        !------------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------------!
        close(unit=i_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit i_f(index)%unit"

        close(unit=f_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit f_f(index)%unit"

        close(unit=o_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit o_f(index)%unit"

        close(unit=c_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit c_f(index)%unit"
        
        close(unit=d_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit d_f(index)%unit"

        close(unit=z_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit z_f(index)%unit"

        close(unit=v_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit v_f(index)%unit"

        close(unit=a_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit a_f(index)%unit"

        close(unit=s_f(index)%unit, iostat=ios)
        if (ios /= 0) stop "Error closing file unit s_f(index)%unit"

        if (allocated(s)) deallocate(s, stat=err3)
        if (err3 /= 0) print *, "s: Deallocation request denied"

        if (allocated(i_f(index)%file_header)) deallocate(i_f(index)%file_header, stat=err)
        if (err /= 0) print *, "file_header: Deallocation request denied"

        if (allocated(z_f(index)%file_header)) deallocate(z_f(index)%file_header, stat=err)
        if (err /= 0) print *, "file_header: Deallocation request denied"

        if (allocated(d_f(index)%file_header)) deallocate(d_f(index)%file_header, stat=err)
        if (err /= 0) print *, "file_header: Deallocation request denied"

        if (allocated(v_f(index)%file_header)) deallocate(v_f(index)%file_header, stat=err)
        if (err /= 0) print *, "file_header: Deallocation request denied"
        !------------------------------------------------------------------------------------------!
    end do month_file_loop

    !----------------------------------------------------------------------------------------------!
    if (allocated(i_f)) deallocate(i_f, stat=err)
    if (err /= 0) print *, "i_f: Deallocation request denied"

    if (allocated(f_f)) deallocate(f_f, stat=err)
    if (err /= 0) print *, "f_f: Deallocation request denied"

    if (allocated(o_f)) deallocate(o_f, stat=err)
    if (err /= 0) print *, "o_f: Deallocation request denied"

    if (allocated(c_f)) deallocate(c_f, stat=err)
    if (err /= 0) print *, "c_f: Deallocation request denied"

    if (allocated(d_f)) deallocate(d_f, stat=err)
    if (err /= 0) print *, "d_f: Deallocation request denied"

    if (allocated(z_f)) deallocate(z_f, stat=err)
    if (err /= 0) print *, "z_f: Deallocation request denied"

    if (allocated(v_f)) deallocate(v_f, stat=err)
    if (err /= 0) print *, "v_f: Deallocation request denied"

    if (allocated(a_f)) deallocate(a_f, stat=err)
    if (err /= 0) print *, "a_f: Deallocation request denied"
    !----------------------------------------------------------------------------------------------!

    call CPU_TIME(stop_epoch)

    write(log_f%unit, *) 'Exited with code=0 in ', stop_epoch - start_epoch, ' seconds.'
    write(*, *) 'Exited with code=0 in ', stop_epoch - start_epoch, ' seconds.'

    close(unit=log_f%unit, iostat=ios)
    if (ios /= 0) stop "Error closing file unit log file"
end program main
