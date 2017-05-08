function petk07()
#
# M. Petkovic testing polynomials, page 147
#
    y = [1*[1,1,1],-2+i, -2-i, 5*i,5*i, -5*i, -5*i];
    p = poly(y);
    z = [1, 3; -2+i,1;  -2-i,1;  5*i,2; -5*i, 2];
    

    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    
   
p,z
end