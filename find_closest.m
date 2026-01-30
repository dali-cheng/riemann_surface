function [result_value, result_index] = find_closest(data_list, target_value)
%FIND_CLOSEST_INDEX Find the closest value of target_value in the data_list
%   We have a data_list and a target_value. We want to find the value in
%   the data_list that is closest to the target_value. Return such value
%   and its index in the data_list.
[~, result_index] = min(abs(data_list - target_value));
result_value = data_list(result_index);
end