function y = convert_to_0_to_a(x, a)
%CONVERT_TO_0_TO_A Map the input x to the interval [0, a) using linear
%modulo. Returns x mod a.
y = x - a * floor(x / a);
end