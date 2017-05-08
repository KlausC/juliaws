function jt10b()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 10^6;
    p = poly([a,1,1/a]);
    z = [[a,1,1/a]',ones(3,1)];
    @printf("\n");
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.8f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;     
   p,z
endfunction jt10b()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 10^6;
    p = poly([a,1,1/a]);
    z = [[a,1,1/a]',ones(3,1)];
    @printf("\n");
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.8f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;     
   p,z
end