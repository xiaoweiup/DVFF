%求轮廓系数    %silhouette factor
%v0.2.0.20220901       在v0.1.0的版本上增加保留每一类簇的质心位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v0.2.1.20220915       在v0.2.0的版本上修复了C原子和簇质心变量名相同的bug

n = 7;
Cluster=zeros(10*n,1001);
for j=1:10
    [idx,CC] = kmeans(O(:,1:1000),n);
    Cluster((j-1)*n+1:j*n,1)=1:n;
    Cluster((j-1)*n+1:j*n,2:1001)=CC;
    O(:,1001) = idx;

    s_sum = 0;
    for i=1:size(O,1)
        tic
        a_cluster = O((O(:,1001)==O(i,1001)),:);
        if size(a_cluster,1)==1
            a=1;
        else

            a=sum(dist(a_cluster(:,1:1000),O(i,1:1000)))/(size(a_cluster,1)-1);
        end

        for k=1:n
            if k==O(i,1001)
                continue
            end
            b_cluster = O((O(:,1001)==k),:);

            b2=1000000000;
            b=sum(dist(b_cluster(:,1:1000),O(i,1:1000)))/size(b_cluster,1);
            if(b<b2)
                b2=b;
            end
        end
        s_sum = s_sum+(b2-a)/max(a,b2);
        toc

    end
    s_sum = s_sum/size(O,1);
    ss_sum(j)=s_sum;
end
