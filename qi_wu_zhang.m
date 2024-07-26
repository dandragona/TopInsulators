syms kx ky t;


H = @(kx,ky,u) [u+cos(kx)+cos(ky), sin(kx)-1i*sin(ky); sin(kx)+1i*sin(ky), -(u+cos(kx)+cos(ky))];

u = 2.1
N1 = 15;
N2 = 15;
kxs = linspace(0,2*pi, N1);
kys = linspace(0,2*pi, N2);

% Plot the energy bands.
energies = zeros(N1,N2,2);
for i = 1:N1
    for j = 1:N2
        [Vs,D] = eig(H(kxs(i),kys(j),u));
        energies(i,j,1) = D(1,1);
        energies(i,j,2) = D(2,2);
    end
end
figure;
hold on
[X, Y] = meshgrid(kys, kxs);
surf(X, Y, energies(:,:,1), 'FaceColor', 'r', 'FaceAlpha', 0.5);
surf(X, Y, energies(:,:,2), 'FaceColor', 'g', 'FaceAlpha', 0.5);  

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
        [Vs,D] = eig(H(kxs(i),kys(j),u));
        [Vs_up,D] = eig(H(kxs(mod(i,N1)+1),kys(j),u));
        [Vs_right,D] = eig(H(kxs(i),kys(mod(j,N2)+1),u));

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
F_first;
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