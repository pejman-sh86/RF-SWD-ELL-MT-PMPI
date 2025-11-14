%===========================================================================
   function matrix_plot(A_in,inorm);
%===========================================================================

% Matlab routine to plot a matrix A_in.
% Normalization: inorm=1 ... normalize matrix as a whole
%                inorm=2 ... normalize each row of the matrix individually

%---------------------------------------------------------------------------

[nrow, ncol] = size(A_in);

A = zeros(nrow,ncol);

if (inorm == 1) 
   norm = max(max(abs(A_in)));
   A = A_in/(norm*2);
end

if (inorm == 2)
   for irow=1:nrow,
      norm = max(abs(A_in(irow,:)));
      A(irow,:) = A_in(irow,:)/(norm*2);
   end
end

xx = 1:ncol;
fill([xx(1),xx,xx(ncol)],[nrow,A(1,:)+nrow,nrow],[0.8,0.8,0.8]);
axis([1 ncol 0 nrow+1])
hold on

for irow=2:nrow,
   fill([xx(1),xx,xx(ncol)],[nrow-irow+1,A(irow,:)+...
      nrow-irow+1,nrow-irow+1],[0.8,0.8,0.8]);
end

for irow=1:nrow,
   plot([1 ncol],[irow irow],':k');
end

for icol=1:ncol,
   plot([icol icol],[0 nrow+1],':k');
end

hold off

return
