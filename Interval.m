classdef Interval
    properties
        l
        h
    end
    methods
        function obj = Interval(l, h)
            obj.l = l;
            obj.h = h;
        end
        function r = plus(a,b)
            if isa(a, 'Interval') && isa(b, 'Interval')
                r = Interval(a.l+b.l, a.h+b.h);
            elseif isa(a, 'Interval') && isa(b, 'double')
                r = Interval(a.l+b, a.h+b);
            elseif isa(a, 'double') && isa(b, 'Interval')
                r = Interval(a+b.l, a+b.h);
            else
                error('Interval addition not supported for %s and %s',class(a),class(b))
            end
        end
        function r = minus(a,b)
            if isa(a, 'Interval') && isa(b, 'Interval')
                r = Interval(a.l-b.h, a.h-b.l);
            elseif isa(a, 'Interval') && isa(b, 'double')
                r = Interval(a.l-b, a.h-b);
            elseif isa(a, 'double') && isa(b, 'Interval')
                r = Interval(a-b.h, a-b.l);
            else
                error('Interval subtraction not supported for %s and %s',class(a),class(b))
            end
        end
        function r = mtimes(A,b)
            A_p = A.*(A>0);
            A_m = A_p - A;
            
            r = Interval(A_p*b.l - A_m*b.h, A_p*b.h - A_m*b.l);
        end
        function r = intersect(a,b)
            r = Interval(max(a.l,b.l), min(a.h,b.h));
        end
    end
end