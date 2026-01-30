function y = convert_to_plus_minus_a(x, a)
%CONVERT_TO_PLUS_MINUS_A Map the input x to the interval [-a, a) using 
%linear modulo.
y = x - 2 * a * floor((x / a + 1) / 2);
end