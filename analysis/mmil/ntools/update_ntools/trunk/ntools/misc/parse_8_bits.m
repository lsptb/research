%
% parse_8_bits(array)
%
% takes an array of eight bits (1,0) and returns their value
% as interpreted as an 8 bit little endian value
%
% example: parse_8_bits([0 0 0 0 0 0 0 1]) == 1
%
function value=parse_8_bits(array)
    value = 0;
    value = value + 2^7*array(1);
    value = value + 2^6*array(2);
    value = value + 2^5*array(3);
    value = value + 2^4*array(4);
    value = value + 2^3*array(5);
    value = value + 2^2*array(6);
    value = value + 2^1*array(7);
    value = value + 1*array(8);
