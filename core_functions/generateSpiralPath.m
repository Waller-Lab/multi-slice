%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a spiral path of quasi-equidistant spacing, spanning values
% zero to one
%
% July 06, 2020 by Shwetadwip Chowdhury
%
% inputs:
%           N:          number of points in spiral path
%           revs:       number of revolutions the spiral path takes
%                       value in .TIF file)
%
%           (Optional inputs)
%           outCirc:    boolean value for whether to include an outer
%                       circle of points after the spiral path
%           N_out:      number of points in the outer circle
%
% outputs:
%           x:          x-coord vector of points on spiral path
%           y:          y-coord vector of points on spiral path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = generateSpiralPath(N, revs, outCirc, N_out)

if nargin <3
    outCirc = false;
end
if ~outCirc
    N_out = 0;
end

N = N-N_out;
t = linspace(0,2*pi,N);
x = t.*cos(revs*t)/2/pi;
y = t.*sin(revs*t)/2/pi;
dist = zeros(size(x));
for k = 2:length(dist)
    dist(k) = sqrt((x(k)-x(k-1)).^2+(y(k)-y(k-1)).^2)+dist(k-1);
end
coef = mean( t(N).^2./dist(N) );
distPerAcq = dist(end)/N;
tInc = sqrt((0:distPerAcq:dist(end)).*coef); tInc = tInc(1:end-1);


x = tInc.*cos(revs*tInc)/2/pi;
y = tInc.*sin(revs*tInc)/2/pi;

if outCirc
    t = linspace(0,2*pi,N_out);
    x = [x cos(t)];
    y = [y sin(t)];
end
