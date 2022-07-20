function f = fator_atrito_hall(Re,Dh)
    
    global rugTub
 
    resultado = 0.001375*(1+((20000*rugTub/Dh)+...
        (1000000/abs(Re)))^(0.333333333333333)) ;
 
end
