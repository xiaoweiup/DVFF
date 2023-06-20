%求轮廓系数    %silhouette factor
%v0.1.0 初代版本

n = 7;
for j=1:10
    idx = kmeans(C(:,1:1000),n);
    C(:,1001) = idx;

    s_sum = 0;
    for i=1:size(C,1)
        tic
        a_cluster = C((C(:,1001)==C(i,1001)),:);
        if size(a_cluster,1)==1
            a=1;
        else

            a=sum(dist(a_cluster(:,1:1000),C(i,1:1000)))/(size(a_cluster,1)-1);
        end

        for k=1:n
            if k==C(i,1001)
                continue
            end
            b_cluster = C((C(:,1001)==k),:);

            b2=1000000000;
            b=sum(dist(b_cluster(:,1:1000),C(i,1:1000)))/size(b_cluster,1);
            if(b<b2)
                b2=b;
            end
        end
        s_sum = s_sum+(b2-a)/max(a,b2);
        toc

    end
    s_sum = s_sum/size(C,1);
    ss_sum(j)=s_sum;
end
