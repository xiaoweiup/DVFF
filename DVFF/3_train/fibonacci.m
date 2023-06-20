function [Sphere_points] = fibonacci(x,y,z)

    samp_num = 1000;
    R = 0.1;

    gold_ratio = (sqrt(5)-1)/2;                                            %0.618
    N = samp_num;

    Sphere_points = zeros(1000,4);
    for n=1:N
        Sphere_points(n,3) = (2*n-1)/N-1;
        Sphere_points(n,1) = sqrt(1-Sphere_points(n,3)^2)*cos(2*pi*n*gold_ratio);
        Sphere_points(n,2) = sqrt(1-Sphere_points(n,3)^2)*sin(2*pi*n*gold_ratio);
    end

    Sphere_points(:,3) = Sphere_points(:,3)*R+z;
    Sphere_points(:,1) = Sphere_points(:,1)*R+x;
    Sphere_points(:,2) = Sphere_points(:,2)*R+y;

end

