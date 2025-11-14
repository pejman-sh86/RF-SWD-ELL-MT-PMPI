function [] = bryan_space();

load picked_tt_int;
t2 = tt(:,3:end);
t1 = tt(:,1:end-2);
r = repmat(r,7,1);

size(r)

r2 = r(:,3:end);
r1 = r(:,1:end-2);

r3 = r(:,2:end-1);


p = (t2-t1)./(r2-r1);


NAVE = 10;
figure(1);hold on;box on;
for j = 1:2

    NRAN = length(p);
    for i=1:NRAN

        if(i <= NAVE)
            NAVE1 = 1;
        else
            NAVE1 = i-NAVE;
        end
        if(i >= NRAN-NAVE)
            NAVE2 = NRAN;
        else
            NAVE2 = i+NAVE;
        end
        p2(j,i) = mean(p(j,NAVE1:NAVE2));

    end


    plot(r3(j,:),p(j,:));
    plot(r3(j,:),p2(j,:),'k');

end


return;

