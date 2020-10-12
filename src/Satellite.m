classdef Satellite
    %   Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pos_i;                %% the position vector in the inertial frame like gcrs
        vel_i;                %% the velocity vector in the inertial frame like gcrs
        pos_e;                %% the position vector in the ecef frame like itrs
        vel_e;                %% the velocity vector in the ecef frame like itrs
        gps_time;             %% the gps time for the object
        right_ascension;      %% the right ascension for this epoch
        negative_declination  %% the negative declination for this epoch
        year;                 %% the year of the epoch in gps_time
        month;                %% the month of the epoch in gps_time
        date;                 %% the date of the epoch in gps_time
        hour;                 %% the hour of the epoch in gps_time
        minute;               %% the minute of the epoch in gps_time
        second;               %% the second of the epoch in gps_time
        time_string;          %% the time string of the epoch in gps_time
        surface_area;         %% the surface area of the satellite
        surface_normal_unit;  %% the normal unit vector for the surfaces of the satellite
        sun_i;                %% the position of the Sun in the inertial frame
        phi_s2n;              %% the angle between the surface normal and the the direction to the Sun
        solar_pressure_c;     %% the three coefficients for solar pressure
        vector_s2s;           %% the unit vector pointing from the satellite to the Sun
        i2s_eul;              %% the euler angles from the inertial frame to srf
        i2s_m;                %% the rotation matrix from the inertial frame to srf
        s2i_eul;              %% the euler angles from srf to the inertial frame
        s2i_m;                %% the rotation matrix from srf to the inertial frame
        euler_angle;          
        rotation_matrix;
    end
    
    methods
        % initiate the Satellite object
        function obj = Satellite(gps_time, pos_i)
            obj.gps_time = gps_time;
            obj.pos_i    = pos_i;
        end
        
        % append the attitude angles from the inertial frame into srf
        function obj = get_i2s_eul(obj, att)
            obj.i2s_eul = att;
        end
        
        % convert the gps time into a classic form
        function obj = gpst2yyyymmddhhmmss(gps_time, year, month, date_init)
            obj.year = year;
            obj.month = month;
            obj.date = floor((gps_time - gps_time(1)) / 86400.) + date_init;
        end
    end
end

