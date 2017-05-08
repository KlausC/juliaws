function fibsq04()
#
#  test polynomial suggested by Goedecker
#
    n = 4;
    p = [-1,ones(1,n)];
    p = conv(p,p);
    z = [-.7748041132154339,         -.7637893113374573e-1-.8147036471703865*i,         -.7637893113374573e-1+.8147036471703865*i,         1.927561975482925];
    z = [z',2*ones(n,1)];
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
endfunction fibsq04()
#
#  test polynomial suggested by Goedecker
#
    n = 4;
    p = [-1,ones(1,n)];
    p = conv(p,p);
    z = [-.7748041132154339,         -.7637893113374573e-1-.8147036471703865*i,         -.7637893113374573e-1+.8147036471703865*i,         1.927561975482925];
    z = [z',2*ones(n,1)];
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