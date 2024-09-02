% Follows Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances
% by Takahiro Fukui, Yasuhiro Hatsugai and Hiroshi Suzuki
% DOI: 10.1143/JPSJ.74.1674

syms kx ky t;

% q = 3
phi = 1/3
H = @(kx,ky,t) [-2*t*cos(ky-2*pi*phi), -t, -t*exp(-1i*3*kx);
    -t, -2*t*cos(ky-4*pi*phi), -t;
    -t*exp(1i*3*kx), -t, -2*t*cos(ky-6*pi*phi)
    ]

N1 = 9;
N2 = 27;
kxs = linspace(0,2*pi/3, N1);
kys = linspace(0,2*pi, N2);

% Plot the energy bands.
energies = zeros(N1,N2,3);
for i = 1:N1
    for j = 1:N2
        [Vs,D] = eig(H(kxs(i),kys(j),1));
        energies(i,j,1) = D(1,1);
        energies(i,j,2) = D(2,2);
        energies(i,j,3) = D(3,3);
    end
end
figure;
hold on
[X, Y] = meshgrid(kxs, kys);
surf(X', Y', energies(:,:,1), 'FaceColor', 'r', 'FaceAlpha', 0.5);
surf(X', Y', energies(:,:,2), 'FaceColor', 'g', 'FaceAlpha', 0.5);  
surf(X', Y', energies(:,:,3), 'FaceColor', 'b', 'FaceAlpha', 0.5);  

hold off;                    
xlabel('kx');                 
ylabel('ky');                 
zlabel('E');                 
title('Energy Bands');   
view(3);

% Calculate Chern Number.
U1 = zeros(N1,N2);
U2 = zeros(N1,N2);
A1 = zeros(N1,N2);
A2 = zeros(N1,N2);
F_first = zeros(N1,N2);
F_second = zeros(N1,N2);
n12 = zeros(N1,N2);
for i = 1:N1
    for j = 1:N2
        [Vs,D] = eig(H(kxs(i),kys(j),1));
        [Vs_up,D] = eig(H(kxs(mod(i,N1)+1),kys(j),1));
        [Vs_right,D] = eig(H(kxs(i),kys(mod(j,N2)+1),1));

        % Gather wave-functions.
        nkl = Vs(:,2);
        nkl_up = Vs_up(:,2);
        nkl_right = Vs_right(:,2);

        % Calculate Link Variable
        U1(i,j) = dot(nkl, nkl_up);
        U1(i,j) = U1(i,j)/norm(U1(i,j));
        U2(i,j) = dot(nkl, nkl_right);
        U2(i,j) = U2(i,j)/norm(U2(i,j));

        % Calculate Gauge Potential
        A1(i,j) = log(U1(i,j));
        A2(i,j) = log(U2(i,j));
    end
end
U1;
U2;
A1;
A2;
c = 0;
for i = 1:N1
    for j = 1:N2

        % Calculate Field Strength (First Way)
        F_first(i,j) = log(U1(i,j)*U2(mod(i,N1)+1,j)/(U1(i,mod(j,N2)+1)*U2(i,j)));
        c = c + F_first(i,j);

        % Calculate Field Strength (Second Way)
        % D1A2 - D2A1
        F_second(i,j) = (A2(mod(i,N1)+1,j) - A2(i,j)) - (A1(i,mod(j,N2)+1)-A1(i,j));
        
        % Calculate n12 field by checking how many multiples of pi
        % we need to shift by.
        if imag(F_second(i,j)) < -pi
           if imag(F_second(i,j)) < -3*pi
                n12(i,j) = 2;
           else
                n12(i,j) = 1;
           end
        elseif imag(F_second(i,j)) > pi
           if imag(F_second(i,j)) > 3*pi
                n12(i,j) = -2;
           else
                n12(i,j) = -1;
           end
        end
    end
end
F_first
F_second;
n12;
c_first = (1/(2*pi*1i))*c
c_second = 0;
for i = 1:N1
    for j=1:N2
        c_second = c_second + n12(i,j);
    end
end
c_second

figure;
hold on
surf(X, Y, real(-1i*F_first)', 'FaceColor', 'y', 'FaceAlpha', 0.5);

hold off;                    
xlabel('kx');                 
ylabel('ky');                 
zlabel('-iF_12');                 
title('Field Strength');   
view(3);
