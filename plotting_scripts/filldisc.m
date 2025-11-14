function y = filldisc(radius)
%FILLDISC  Returns an array with a disc full of ones.
%    Y = FILLDISC(R)  Returns a square array of size 2*ROUND(R)+1.
%    Elements less than R from the central elements are set to 1, rest to
%    0. This is done fairly efficiently, taking advantage of symmetries.

% Copyright David Young 2010

rsq = radius * radius;
R = round(radius);
c = R + 1;
t = 2*R + 1;
y = zeros(t);

% points on St George's cross
y(:, c) = 1;
y(c, :) = 1;

for rstart = 1:R
    rend = round(sqrt(rsq - rstart * rstart));
    if rend < rstart
        break
    end

    % points on St Andrew's cross
    y(c+rstart, c+rstart) = 1;
    y(c+rstart, c-rstart) = 1;
    y(c-rstart, c+rstart) = 1;
    y(c-rstart, c-rstart) = 1;

    % fill in octants
    y(c+rstart:c+rend, c+rstart) = 1;
    y(c+rstart:c+rend, c-rstart) = 1;
    y(c-rend:c-rstart, c+rstart) = 1;
    y(c-rend:c-rstart, c-rstart) = 1;
    y(c+rstart, c+rstart:c+rend) = 1;
    y(c-rstart, c+rstart:c+rend) = 1;
    y(c+rstart, c-rend:c-rstart) = 1;
    y(c-rstart, c-rend:c-rstart) = 1;
end

end

