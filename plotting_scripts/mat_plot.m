function mat_plot(A,col);

n_i = size(A,1);
n_j = size(A,2);

%for i = 1:n_i
%  for j = 1:n_j
%    mcor(i,j) = mcov(i,j)/sqrt(mcov(i,i)*mcov(j,j));
%  end
%end

B = A+repmat(((2*n_i)-1:-2:0)',1,n_i);

figure(1);
hold on;grid off;
for i = 1:n_i
    for j = i:n_j
        if(abs(A(i,j))>0.3)
            if(i<B(i,j))
                (n_i-i)*2+1
                plot([j j],[(n_i-i)*2+1 B(i,j)],'-b','LineWidth',3);
            else
                plot([j j],[B(i,j) (n_i-i)*2+1],'-b','LineWidth',3);
            end
            if(abs(A(i,j))~=1)
                tmp = A(i,j);
                fprintf(1,'%2i\t %2i\t %8.6f\n',i,j,tmp);
            end
            plot([1 7],[(n_i-i)*2+1 (n_i-i)*2+1],'-k');
            plot([1 7],[(n_i-i)*2 (n_i-i)*2],':k');
        end;
    end;
end;
plot([1 7],[13 1],'-b');
figure(2);
hold on;grid on;
for i = 1:n_i
    plot(B(i,:),col);
end;
set(gca,'XTick',[1:1:n_j])
set(gca,'YTick',[1:2:n_i*2])
%saveas(gca,plotfile1,'epsc2');

return

