classdef Motion2D
    properties
        x %{mustBeNumeric}
        y
        r
        theta
    end
    methods(Static)
        
    end
    methods
        %%% get
        function x = getX(obj)
            x = obj.x;
        end
        function y = getY(obj)
            y = obj.y;
        end
        function r = getR(obj)
            r = obj.r;
        end
        function theta = getT(obj)
            theta = obj.theta;
        end
        function [x, y] = getXY(obj)
            x = obj.x;
            y = obj.y;
        end
        function [r, theta] = getRT(obj)
            r = obj.r;
            theta = obj.theta;
        end
        
        %%% set
        function obj = setXY(obj, x, y)
            obj.x = x;
            obj.y = y;
            [obj.r, obj.theta] = obj.cartesian2polar(x, y);
        end
        function obj = setRT(obj, r, theta)
            if (r < 0)
                obj.r = -r;
                obj.theta = theta + pi;
            else
                obj.r = r;
                obj.theta = theta;
            end
            [obj.x, obj.y] = obj.polar2cartesian(r, theta);
            [obj.r, obj.theta] = obj.cartesian2polar(obj.x, obj.y);
        end
        
        %%% utility methods
        function [r, theta] = cartesian2polar(~, x, y)
            
            if (y > 0 && x > 0)     % first quadrant
                theta = atan(y / x);
            elseif (y > 0 && x < 0) % second quadrant
                theta = atan(y / x) + pi;
            elseif (y < 0 && x < 0) % third quadrant
                theta = atan(y / x) - pi;
            else                    % forth quadrant
                theta = atan(y / x);
            end
           
            r = sqrt( x^2 + y^2 );
        end

        function [x, y] = polar2cartesian(~, r, theta)
            x = r * cos(theta);
            y = r * sin(theta);
        end
        function obj = plus(obj1, obj2)
            obj = Motion2D;
            obj = obj.setXY(obj1.getX() + obj2.getX(), obj1.getY() + obj2.getY());
        end
        function obj = minus(obj1, obj2)
            obj = Motion2D;
            obj = obj.setXY(obj1.getX() - obj2.getX(), obj1.getY() - obj2.getY());
        end
        function obj = scale(self, coeff)
            r_temp = self.getR();
            t = self.getT();
            r_temp = r_temp * coeff;
            obj = self.setRT(r_temp, t);
        end
    end
end