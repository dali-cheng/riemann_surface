function y = lorentzian_fittype(x, x0, delta_x, gamma, BL, D, num_res_loc_int)
y = ones(size(x));
for loop_index = 0 : (num_res_loc_int - 1)
    y = y - D ./ (1 + ((x - x0 - loop_index * delta_x) / gamma).^2);
end
y = y * BL;
end