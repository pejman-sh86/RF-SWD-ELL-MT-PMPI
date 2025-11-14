function plot_unitlag_cor(filename);

C1 = load(filename);
for i=1:size(C1,1);
    for j=1:size(C1,2);
        R1(i,j) = C1(i,j)/sqrt(C1(i,i)*C1(j,j));
    end;
end;

figure();
matrix_plot(C1,2);
figure();
matrix_plot(R1,2);

return;