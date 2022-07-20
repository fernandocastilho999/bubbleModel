function f = fator_atrito_blasius(Re)
    
    %global rhoL miL rhoG tensup D g
   
    if (Re < 2300)
        f = 16/Re ;
    else
        f = 0.3164/4/Re^0.25 ;
    end
end
