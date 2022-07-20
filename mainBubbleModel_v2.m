clc, clear all, close all

D = 0.025 ;
inclinacao = 0 ;
inclinacaor = inclinacao * pi/180 ;

rugTub = 0 ;
g = 9.81 ;

rhoL = 999.0 ;
miL = 0.000855 ;

rhoG = 1.2 ;
miG = 17.0e-6;


ene = 0.0 ;
tensup = 0.07 ;

Co = 1.12 ;

jG = 1;
jL = 10;

while jG <=10

    
RS = 1;
J = jG + jL ;

%RLS = fun_rlsGregory(J);
RLS = 1;
f = fun_freqSchulkes(rhoL, miL, D,  inclinacaor, g, jL, J); %Schulkes(2011)

%global rhoL Co RLS miL rhoG miG tensup f D g rugTub inclinacaor jL jG J RS ene


contlinhas = 0 ;
linhas = 1;
while (contlinhas < linhas)
    contlinhas = contlinhas+1 ;
    if (inclinacao >= 50)
        %call Vertical (jGk(contlinhas),jLk(contlinhas),RSk(contlinhas))
        [res1,res2,LB, LS, RGB] = vertical(Co, RLS, rhoL, miL, rhoG, miG, tensup, ...
            D, g, inclinacaor, jL, J, f, ene);
    else
        
        [res1,res2,LB, LS, RGB] = horizontal(Co, RLS, rhoL, miL, rhoG, miG, ...
            tensup, D, g, inclinacaor, jL, J, f, ene);
    end
    jG = jG + 1

end
figure(01)
subplot(2,1,1)
plot(res1,res2)
hold on
grid on
title('Modelo de Bolha para JL=10m/s e JG=[1:10]m/s (Comprimento)')
legend('JL=1.0','JL=2.0','JL=3.0','JL=4.0','JL=5.0','JL=6.0','JL=7.0','JL=8.0','JL=9.0','JL=10.0')
xlabel('Comprimento (m)')
ylabel('Fração de Vazio')

subplot(2,1,2)
hold on
plot(res1/LB+1,res2)
grid on
title('Modelo de Bolha para JL=10m/s e JG=[1:10]m/s (Adimensional)')
legend('JL=1.0','JL=2.0','JL=3.0','JL=4.0','JL=5.0','JL=6.0','JL=7.0','JL=8.0','JL=9.0','JL=10.0')
xlabel('Comprimento Adimensional')
ylabel('Fração de Vazio')

end



% 11 : perfil.txt
% 12 : resultados.txt