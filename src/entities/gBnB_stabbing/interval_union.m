function [lower, upper] = interval_union(left, right)
% function [lower upper] = MergeBrackets(left, right)
%
% Purpose: Interval merging
%
% Given N input closed intervals in braket form:
%   Ii := [left(i),right(i)], i = 1,2...,N (mathematical notation)
% The set union{Ii) can be written as a canonical partition by
%   intervals Jk; i.e., union{Ii) = union(Jk), where Jk are M intervals
%   (with M<=N, so the partition is minimum cardinal), and {Jk} are
%   disjoint to each other (their intersections are empty). This function
%   returns Jk = [lower(k),upper(k)], k=1,2,...M, in the ascending sorted
%   order.
%
% EXAMPLE USAGE:
%   >> [lower upper] = MergeBrackets([0 1 2 3 4],[1.5 1.6 3.5 3 5])
%   	lower =   0    2  4
%       upper = 1.6  3.5  5
%
% Algorithm complexity: O(N*log(N))
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Original: 25-May-2009

% sort the rest by left bound
[left, iorder] = sort(left);
right = right(iorder);
% Allocate, as we don't know yet the size, we assume the largest case
lower = [];
upper = [];

% Initialize
l = left(1);
u = right(1);
k = 0;
% Loop on brakets
for i=1:4
    if left(i) > u % new Jk detected
        % Stack the old one
        k = k+1;
        lower = [lower, l];
        upper = [upper, u];
        % Reset l and u
        l = left(i);
        u = right(i);
    else
        u = max(u, right(i));
    end
end
% Stack the last one

lower = [lower, l];
upper = [upper, u];

end
