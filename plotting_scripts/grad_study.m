function [] = grad_study();


load('gradient_study.mat');

for j=1:length(freq); [ref(:,j)]=ref_nlay_gfav3(thd,model,freq(j),frbw); end;

F(1).model = model_1;
for i = 2:31;

  F(i).model = [F(i-1).model(1,:) ; F(i-1).model(2,:);F(i-1).model(2:end,:)];
  F(i).model(2:i+1,1) = 0.31/i;
  for j=1:i; F(i).model(j+1,4)=1.3+(0.3/(i-1))*(j-1); end;

end;

for i = 1:31;

  for j=1:length(freq); [F(i).ref(:,j)]=ref_nlay_gfav3(thd,F(i).model,freq(j),frbw); end;
  F(i).res = F(i).ref-ref;

end;

figure(1);
hold on;
plot(thd,ref(:,4),'k');
for i = 1:31-1;
  plot(thd,F(i).ref(:,4),'--r');
end;

figure(2);
hold on;
for i = 1:31-1;
  plot(thd,F(i).res(:,4),'--r');
end;

save bla;

return;
