function henrici()
#
# M. Petkovic testing polynomials, page 123
#
    p = poly([-4.1,-3.8,-2.05,-1.85,1.95,2.15,3.9,4.05]);
    z = [[-4.1,-3.8,-2.05,-1.85,1.95,2.15,3.9,4.05]',ones(8,1)];

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
endfunction henrici()
#
# M. Petkovic testing polynomials, page 123
#
    p = poly([-4.1,-3.8,-2.05,-1.85,1.95,2.15,3.9,4.05]);
    z = [[-4.1,-3.8,-2.05,-1.85,1.95,2.15,3.9,4.05]',ones(8,1)];

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