function rls = rls_Gregory
    
    global rhoL miL rhoG tensup D g J

    ReM = rhoL* J * D/miL 
    if (ReM < 2300)
        fM = 16D0/ReM ;
    else
        fM = 0.3164/4./ReM.^0.25 ;
    end
    rls = 1-0.058*(2*(0.4*tensup/(rhoL-rhoG)/g)^0.5*...
        ((2*fM/D)*J^3)^(2/5)*(rhoL/tensup)^(3/5)-0.725)^2 ;
end
