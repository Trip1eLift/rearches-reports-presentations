classdef SpaceCraft
    properties
        position
        velocity
        fuel_burn_duration
        v_exhaust
        mass
        dm
        dt
    end
    
    methods
        function obj = properties(obj, Time_max, dt, frame_mass, fuel_mass, fuel_burn_duration, v_exhaust)
            obj.fuel_burn_duration = fuel_burn_duration;
            obj.v_exhaust = v_exhaust;
            obj.dt = dt;
            
            Length = Time_max / dt;
            
            position_temp(1, Length) = Motion2D;
            velocity_temp(1, Length) = Motion2D;
            obj.position = position_temp;
            obj.velocity = velocity_temp;
            
            obj.dm = fuel_mass / fuel_burn_duration * dt;
            obj.mass = zeros(1, Length);
            obj.mass(1) = frame_mass + fuel_mass;
            for i = 1:Length - 1
                if (i <= fuel_burn_duration / dt)
                    obj.mass(i + 1) = obj.mass(i) - obj.dm;
                else
                    obj.mass(i + 1) = obj.mass(i);
                end
            end
        end
        
        function obj = initial_positionXY(obj, x, y)
            obj.position(1) = obj.position(1).setXY(x, y);
        end
        
        function obj = initial_positionRT(obj, r, theta)
            obj.position(1) = obj.position(1).setRT(r, theta);
        end
        
        function obj = initial_velocityXY(obj, x, y)
            obj.velocity(1) = obj.velocity(1).setXY(x, y);
        end
        
        function obj = initial_velocityRT(obj, r, theta)
            obj.velocity(1) = obj.velocity(1).setRT(r, theta);
        end
        
        function acc = upward_aceleration_at(obj, i)
            % Special case: When the rocket is lauched, it's launched
            % vertically
            acc = Motion2D;
            acc_mag = abs( -obj.v_exhaust * obj.dm / obj.mass(i) / obj.dt );
            acc = acc.setRT( acc_mag, obj.position(i).getT() );
        end
        
        function acc = forward_aceleration_at(obj, i)
            % General case: The rocket pushed itself forward
            acc = Motion2D;
            acc_mag = abs( -obj.v_exhaust * obj.dm / obj.mass(i) / obj.dt );
            v_mag = sqrt( obj.velocity(i).getX()^2 + obj.velocity(i).getY()^2 );
            acc = acc.setXY( obj.velocity(i).getX() / v_mag * acc_mag, obj.velocity(i).getY() / v_mag * acc_mag );
        end
    end
end
        