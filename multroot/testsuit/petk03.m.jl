function petk03()
#
# M. Petkovic testing polynomials, p134
#
    p1 = [1,-1-99*i/70];
    p2 = [1,2,3];
    p2 = conv(p2,p2);
    p = conv(p1,p2);
    p = conv(p,[1,1]);
    z = [-1.00000000000000 + 1.41421356237309i, 2;    -1.00000000000000 - 1.41421356237309i, 2;     -1-99*i/70,  1; -1, 1];

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