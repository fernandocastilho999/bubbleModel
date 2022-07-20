function [res1,res2,LB, LS, RGB] = vertical(Co, RLS, rhoL, miL, rhoG, miG, tensup, ...
    D, g, inclinacaor, jL, J, f, ene)

%global Co RLS rhoL miL rhoG miG tensup D g inclinacaor jL J ene


dh = D/1e3 ;
primeira = 0 ;

%perfil(1,:) = [0 0.5] ;


%%%%%%%%%%%%%%%%%%%%
%RLS = fun_rlsGregory(J) ;

RLBi = RLS ;
alturaFilme = D/2*(1-sqrt(1-RLBi));

Vdrift = 1.54*(tensup*g*(rhoL-rhoG)/rhoL^2 )^0.25*RLS^ene*sin(inclinacaor) ; ...
    %alterado para 9.81 -> g EQ 56 E 57


UGS = 1*J + Vdrift ; %Conferir
ULS = J - UGS*(1-RLS) ;
VD = 0.54*sqrt(g*D)*cos(inclinacaor)+0.35*sqrt(g*D)*sin(inclinacaor) ; %Bendiksen 1984
UT = Co*J + VD ;
altB = alturaFilme ;

somaULB = 0. ;
voltotal = 0. ;
invIncDB = 0. ;
somaRLB = 0. ;
compB = 0. ;
loop = 1 ; 
incdBolha = 0. ;

k = 0 ;
condicao = 1 ;
ii = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (condicao == 1) %Loop de integração da bolha
        
    if (RLS >= 1)
        UGS = 0 ;
    else
        Vdrift = 1.54*(tensup*g*(rhoL-rhoG)/rhoL^2 )^0.25 ...
            *RLS^ene*sin(inclinacaor) ;
        UGS = J + Vdrift ;
    end

    RLB = 1 -(1 -2 *altB/D)^2 ;  %Vertical

    %Cálculo interativo dos parâmetros associados ao desenho da bolha
    ULS =  J - UGS*(1-RLS) ;

    if (loop == 1)
        ULBo = UT-((UT-ULS)*RLS/RLB) ;
        ALo = RLB*pi*D^2 /4 ;  %Vertical
    else
        ULBo = ULB ;
        ALo = AL ;
    end
    ULB = UT-((UT-ULS)*RLS/RLB) ; %Velocidade do líquido na bolha

    if (RLB == 1)
        UGB = 0 ;
    else
        UGB =  UT - ((UT-UGS)*(1 -RLS)/(1 -RLB)) ; %Velocidade do gás na bolha
    end

    %Líquido
    SL = pi*D ; %Vertical
    AL = RLB*pi*D^2 /4 ;  %Vertical
    DhL = 4 *AL/SL ; %Diametro hidráulico de liquido
    ReL = rhoL*ULB*DhL/miL ; %Reynolds do liquido
    fL = fun_fBlasius(ReL) ;
    %fl = fator_atrito_hall(Rel,dhl)
    talL = rhoL*fL*abs(ULB)*ULB/2  ;

    %Gas
    SG = 0 ; %Vertical - sem atrito gás e parede
    AG = (1-RLB)*pi*D^2 /4 ;  %Vertical
    SI = pi*(D-2 *altB) ; %Vertical, alterado 19/01/2017 de pi*(D-altB)
    DhG = 4 *AG/(SG+SI) ;
    ReG = rhoG*UGB*DhG/miG ;
    fG = 0 ; %Vertical
    talG = rhoG*fG*abs(UGB)*UGB/2 ;  %Vai zerar

    %Interfacial
    fI = 0.005 *(1 +300 *altB/D) ; %Wallis 1969 - recomendada por Tailtel e Barnea (1990)
    talI = rhoG*fI*abs(UGB-ULB)*(UGB-ULB)/2 ; 

    dR = 4 *(1 -2 *altB/D)/D ; %Derivada dRLB/dh
    termoFilmeL = talL*SL/AL ;
    if (RLB == 1)   %sE RLB = 1, não existe gás.
        termoFilmeG = 0 ;
    else
        termoFilmeG = talG*SG/AG ;
    end
    termoFilmeI = talI*SI/AL ;
    termoGravit = (rhoL - rhoG ) * g * sin(inclinacaor) ;
    inerciaL = rhoL * (UT-ULB)^2  * dR /  RLB ;
    if (RLB == 1) 
        inerciaG = 0 ;
    else
        inerciaG = rhoG  * (UT-UGB)^2  * dR/ (1-RLB) ;
    end
    forHidro = (rhoL - rhoG ) * g * cos(inclinacaor) ; %Alterado 19/01/2017 de 9.81 para g - caso queira mudar a gravidade
    numerador = termoFilmeL - termoFilmeG - termoFilmeI + termoGravit ; %Bateu
    denominador = forHidro - inerciaG - inerciaL ;
    incdBolha   =   numerador / denominador ;

    %Ponta da bolha quadrada
    if(incdBolha > 0  && primeira == 0)  
        ii = ii + 1 ;
        res1(ii)=compB;
        res2(ii)=altB/D ;
        altB = altB - dh ;
        somaULB = 0  ;
        voltotal = 0  ;
        invIncDB = 0  ;
        somaRLB = 0 ; 
        compB = 0  ;
        incdBolha = 0 ;        
        loop = loop + 1 ;

    else
        primeira = 1 ; %variável para não voltar para a parte "quadrada" da ponta da bolha
        if(incdBolha < 0) 
            verificador =  abs(atand(numerador/denominador)) ; 
            dCompB      =   invIncDB * dh ; %Ponto anterior
            invIncDB    =   1 /incdBolha ; %Novo ponto - Apenas usado no próximo loop

            %Integração Dos Parametros
            prodVelFilm =   ((ULB+ULBo)/2 )*((AL+ALo)/2 )*dCompB ;
            somaULB     =   somaULB + prodVelFilm ;
            volsecao    =   ((AL+ALo)/2 *dCompB) ;
            voltotal = voltotal + volsecao ;
            dCompRLB    =  dCompB * RLB ;
            somaRLB        =   somaRLB + abs(dCompRLB) ;
            LB = abs(compB) ;
            LS = UT/f-LB ;
            
            if (voltotal == 0)  
                ULBm = 0 ;
                RLBm = 0 ;

            else
                ULBm = abs(somaULB/voltotal) ;
                RLBm = abs(somaRLB / LB) ;
            end

            betam = 1 -(jL-ULBm*RLBm)/UT/RLS/(1 -RLBm/RLS) ;
            beta = LB/(LB+LS) ;
            ebeta = (betam-beta)*100 /beta ; %ebeta = (beta-betam)*100/beta

            %Critério de Parada
            if (abs(ebeta) > 0.1 )  
                condicao = 1 ;
            else
                condicao = 0 ;
            end

            if (ebeta < 0)   % if ebeta > 0
                if(k == 1)  
                    compB = compB - dCompB ;
                    altB = altB + dh ;
                    somaULB = somaULB - prodVelFilm ; %VOLTAR O PONTO
                    voltotal = voltotal - volsecao ;
                    somaRLB    = somaRLB - abs(dCompRLB) ;
                    k = 0 ;
                end
                dh = dh/10 ;
            else
                loop = loop + 1 ;
                compB = compB + invIncDB* dh ;
                altB  = altB - dh ;
                k = 1 ;
            end               
        else
            if(k == 1)
                compB = compB - invIncDB* dh ;
                altB = altB + dh ;
                somaULB = somaULB - prodVelFilm ;
                voltotal = voltotal - volsecao ;
                somaRLB    = somaRLB - abs(dCompRLB) ;
                loop = loop - 1 ;
                k = 0 ;
            end
            dh = dh/10  ;             
        end
    end
end
RGB = 1 - RLBm ;



