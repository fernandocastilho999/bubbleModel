function freq = freq_schulkes(rhoL, miL, D,  inclinacaor, g, jL, J) 
   
    %global rhoL miL D  inclinacaor g jL J

    alpha = jL./J ;
    p1 = 0.016*alpha*(2+3*alpha) ;
    Rel = rhoL*jL*D/miL ;
    
    if (Rel >= 4000D0)
        p2 = 1 ;
    else
        p2 = 12.1*Rel^-0.37 ;
    end
    Fr = jL./sqrt(D*g*cos(inclinacaor)) ;
    
    if (abs(inclinacaor) <= 0.17)
        p3 = 1+2*sign(inclinacaor).*sqrt(abs(inclinacaor))./Fr ;
    else
        p3 = 1.8*(0.6+2*inclinacaor-inclinacaor^2)/Fr ;
    end
    
    freq = J*p1*p2*p3/D ;     
end
