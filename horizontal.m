function [res1,res2,LB, LS, RGB] = horizontal(Co, RLS, rhoL, miL, rhoG, miG, ...
    tensup, D, g, inclinacaor, jL, J, f, ene)

%global Co RLS rhoL miL rhoG miG tensup D g inclinacaor jL jG J f ene

% DADOS DE ENTRADA %
%    inclinacaor = inclinacao*(PI)/180 

%Passo
dh = D/1e3 ;
primeira = 0 ;

%Velocidades Superficiais



alturaFilme = D*0.9  ;  %0.9 altura inicial
Vdrift = 1.54 *(tensup*g*(rhoL-rhoG)/rhoL^2  )^0.25 *RLS^ene*sin(inclinacaor) ;
UGS = J + Vdrift ;
ULS = J - UGS*(1-RLS) ;
%Co = CoM
VD = 0.54 *sqrt(g*D)*cos(inclinacaor)+0.35 *sqrt(g*D)*sin(inclinacaor) ; %Bendiksen 1984 
UT = Co*J + VD ;
%    CInf = CInfM
%    VD = CInf * sqrt(g*D)
%    UT = Co*J + VD
altB = alturaFilme ;
somaULB     = 0  ;
voltotal     = 0 ; 
invIncDB    = 0 ;
somaRLB        = 0 ;
compB       = 0 ;
loop        = 1   ;
incdBolha = 0 ;
k = 0 ;
ii=0;

condicao = 1 ; %Condição para entrada do primeiro loop
while (condicao == 1)
   

    if (RLS >= 1 ) 
        UGS = 0 ;
    else
        Vdrift = 1.54 *(tensup*g*(rhoL-rhoG)/rhoL^2  )^0.25 *RLS^ene*sin(inclinacaor);
        UGS = J + Vdrift ;
    end
    angInterno = 2 *acos(1 -2 *altB/D) ;
    RLB = (angInterno - sin(angInterno))/2 /pi ;

    %Cálculo interativo dos parâmetros associados ao desenho da bolha
    ULS =  J - UGS*(1-RLS) ;
    if (loop == 1) 
        ULBo = UT-((UT-ULS)*RLS/RLB) ;
        ALo = D^2 *(angInterno - sin(angInterno))/8  ;
    else
        ULBo = ULB;
        ALo = AL;
    end
    ULB = UT-((UT-ULS)*RLS/RLB) ; %Velocidade do líquido na bolha
    UGB =  UT - ((UT-UGS)*(1 -RLS)/(1 -RLB)) ; %Velocidade do gás na bolha

    %Liquido
    SL = angInterno*D/2 ; %Perimetro molhado do líquido
    AL = D^2 *(angInterno - sin(angInterno))/8  ; %Área do líquido
    DhL = 4 *AL/SL ; %Diametro hidráulico de liquido
    ReL = rhoL*ULB*DhL/miL ; %Reynolds do liquido
    fL = fun_fBlasius(ReL) ;
    %fl = fator_atrito_hall(Rel,dhl)       
    talL = rhoL*fL*abs(ULB)*ULB/2 ;

    %Gas
    SG = (2 *pi-angInterno)*D/2 ;  %Perimetro do Gas
    AG = D^2 *(2 *pi- angInterno + sin(angInterno))/8 ; %Area do gás
    SI = D*sin(angInterno/2 ) ; %Pertimetro interfacial
    DhG = 4 *AG/(SG+SI) ; %Diametro Hidráulico de gas
    ReG = rhoG*UGB*DhG/miG ; %Reynolds do gas
    fG = fun_fBlasius(ReG)  ;
    %fg = fator_atrito_hall(ReG,dhg)        
    talG = rhoG*fG*abs(UGB)*UGB/2 ;

    %Interfacial
    fI = 0.014 ; %Fator de atrito interfacial - Sugerido por: Cohen e Hanratty(1968); Shoham e Tailtel(1984) - Citado no artigo de Tailtel e Barnea(1990)
    talI = rhoG*fI*abs(UGB-ULB)*(UGB-ULB)/2 ;
    dR = 4 *sin(angInterno/2 )/(pi*D) ; %Derivada dRLB/dh
    termoFilmeL = talL*SL/AL ;
    termoFilmeG = talG*SG/AG ;
    termoFilmeI = talI*SI/((AL*AG/(AL+AG))) ;
    termoGravit = (rhoL - rhoG ) * g * sin(inclinacaor) ;
    inerciaL = rhoL * (UT-ULB)^2  * dR /  RLB   ;
    inerciaG = rhoG  * (UT-UGB)^2  * dR/ (1-RLB) ;
    forHidro = (rhoL - rhoG ) * g * cos(inclinacaor) ;
    numerador = termoFilmeL - termoFilmeG - termoFilmeI + termoGravit ;
    denominador = forHidro - inerciaG - inerciaL ;
    incdBolha   =   numerador / denominador ; %dh/dz
    %write(*,*),"Inclinacao",rls,rlsg,rlsb
    %pause


    %PARA CHEGAR NO NARIZ DA BOLHA   
    if(incdBolha > 0)  
        altB = altB - dh  ;
        somaULB     = 0  ;
        voltotal     = 0  ;
        invIncDB    = 0  ;
        somaRLB        = 0 ; 
        compB       = 0  ;
        loop        = 1   ;
        incdBolha = 0  ;
        %clear perfil

    %DH AJUSTÁVEL 
    else

        if(incdBolha < 0 ) 
            ii = ii +1 ;
            verificador =  abs(atand(numerador/denominador)) ;
            %write(11,*),compB,altB/D  %Plota o perfil
            res1(ii) = compB ;
            res2(ii) = altB/D  ;
            dCompB      =   invIncDB * dh  ; %Ponto anterior
            invIncDB    =   1 /incdBolha  ; %Novo ponto - Apenas usado no próximo loop

            %Integração Dos Parametros
            prodVelFilm =   ((ULB+ULBo)/2 )*((AL+ALo)/2 )*dCompB  ;
            somaULB     =   somaULB + prodVelFilm  ;
            volsecao    =   ((AL+ALo)/2 *dCompB) ;
            voltotal = voltotal + volsecao ;
            dCompRLB    =  dCompB * RLB ;
            somaRLB        =   somaRLB + abs(dCompRLB)                  ;
            LB = abs(compB) ;
            LS = UT/f-LB ;
            if (voltotal == 0 )   
                ULBm = 0  ;
                RLBm = 0  ;

            else
                ULBm = abs(somaULB/voltotal) ;
                RLBm = abs(somaRLB / LB) ;
            end

            %Cálculo do erro
            betam = 1 -(jL-ULBm*RLBm)/UT/RLS/(1 -RLBm/RLS) ;
            beta = LB/(LB+LS) ;
            ebeta = (betam-beta)*100 /beta  ; %ebeta = (beta-betam)*100/beta


            %Critério de Parada
            if (abs(ebeta) > 0.1 )  
                condicao = 1 ;
            else
                condicao = 0 ;
            end


            if (ebeta < 0 )    %Se erro < 0 não existe solução real
                if(k == 1)   %Verifica se precisa voltar 1 ponto
                    %Ponto novo = Ponto antigo
                    compB = compB - dCompB ;
                    altB = altB + dh ;
                    somaULB = somaULB - prodVelFilm  ;
                    voltotal = voltotal - volsecao ;
                    somaRLB    = somaRLB - abs(dCompRLB) ;
                    k = 0 ;
                end
                dh = dh/10 ; %Refina o passo
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
           % pause


        end
    end


end

RGB = 1 - RLBm     ;