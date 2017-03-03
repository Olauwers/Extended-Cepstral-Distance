function [inputs,outputs,clusterid] = generateelectriccircuits(N)

    n = 15;
    tolerance = 10^(-3);
    copies = 50;

    R1 = 100;
    L11 = 20;
    L12 = 60;
    C1 = 50;
    R2 = 100;
    L21 = 200;
    L22 = 160;
    C2 = 75;
    A1 = [0,1/L12,0;-1/C1,0,-1/C1;0,1/L11,-R1/L11];
    A2 = [0,1/L22,0;-1/C2,0,-1/C2;0,1/L21,-R2/L21];
    B1 = [0;1/C1;0];
    B2 = [0;1/C2;0];
    C1 = [0 1 -R1];
    C2 = [0 1 -R2];
    D = 0;
    sys1 = ss(A1,B1,C1,D);
    sys2 = ss(A2,B2,C2,D);
    sys1d = idpoly(c2d(sys1,1));
    sys2d = idpoly(c2d(sys2,1));
    
    inputs = zeros(4*copies,N);
    outputssys1 = zeros(4*copies, N);
    outputssys2 = zeros(4*copies, N);
    clusterid = zeros(4*copies,1);

    for j = 1:copies
    
        input1 = idpoly(drss(n));
        [z, p, k] = zpkdata(input1,'v');
        while any(abs(abs(z)-1)<tolerance) || any(abs(abs(p)-1)<tolerance)
            input1 = drss(n);
            input1 = idpoly(input1);
            [z, p, k] = zpkdata(input1,'v');
            disp('Retried input1 because of zero or pole close to 1')
        end

        input2 = idpoly(drss(n));
        [z, p, k] = zpkdata(input2,'v');
        while any(abs(abs(z)-1)<tolerance) || any(abs(abs(p)-1)<tolerance)
            input2 = drss(n);
            input2 = idpoly(input2);
            [z, p, k] = zpkdata(input2,'v');
            disp('Retried input2 because of zero or pole close to 1')
        end


    inputs(4*(j-1)+1,:) = 1*sim(input1,100*rand*randn(N,1));
    inputs(4*(j-1)+2,:) = 1*sim(input2,100*rand*randn(N,1));
    inputs(4*(j-1)+3,:) = 10*(idinput(N,'sine')+0.1*randn(N,1));
    inputs(4*(j-1)+4,:) = 100*rand*idinput(N,'rgs');


    outputssys1(4*(j-1)+1,:) = sim(sys1d,inputs(4*(j-1)+1,:)');
        clusterid(4*(j-1)+1,:) = 1;
    outputssys1(4*(j-1)+2,:) = sim(sys1d,inputs(4*(j-1)+2,:)');
        clusterid(4*(j-1)+2,:) = 1;
    outputssys1(4*(j-1)+3,:) = sim(sys1d,inputs(4*(j-1)+3,:)');
        clusterid(4*(j-1)+3,:) = 1;
    outputssys1(4*(j-1)+4,:) = sim(sys1d,inputs(4*(j-1)+4,:)');
        clusterid(4*(j-1)+4,:) = 1;

    outputssys2(4*(j-1)+1,:) = sim(sys2d,inputs(4*(j-1)+1,:)');
        clusterid(4*copies + 4*(j-1)+1,:) = 2;
    outputssys2(4*(j-1)+2,:) = sim(sys2d,inputs(4*(j-1)+2,:)');
        clusterid(4*copies + 4*(j-1)+2,:) = 2;
    outputssys2(4*(j-1)+3,:) = sim(sys2d,inputs(4*(j-1)+3,:)');
        clusterid(4*copies + 4*(j-1)+3,:) = 2;
    outputssys2(4*(j-1)+4,:) = sim(sys2d,inputs(4*(j-1)+4,:)');
        clusterid(4*copies + 4*(j-1)+4,:) = 2;

    disp(j)

    end

    outputs = [outputssys1;outputssys2];
    inputs = [inputs;inputs];

end