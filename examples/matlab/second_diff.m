function K = second_diff(n)
% Input: n: a nonnegative integer 2 or greater
% Output: a second difference matrix
% K = [ 2  -1            ]
%     [-1   2  -1        ]
%     [    ..   ..  ..   ]
%     [      ..  ..    -1]
%     [             -1  2]
    K = toeplitz([2 -1 zeros(1, n - 2)]);
end