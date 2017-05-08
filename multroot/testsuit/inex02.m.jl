function inex02()
#
#  test polynomial suggested by Goedecker
#
    p = poly([(10/11)*[1,1,1,1,1],(20/11)*[1,1,1],            (30/11)*[1,1]]);
    p = round(10^6*p)/10^6;
    z = [10/11, 5; 20/11, 3; 30/11, 2];
    z = z(3:-1:1,:);
    @printf("\n");
    @printf(" Coefficients are rounded up at the 6-th digits after\n");
    @printf("    decimal point. Originally, \n");
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
endfunction inex02()
#
#  test polynomial suggested by Goedecker
#
    p = poly([(10/11)*[1,1,1,1,1],(20/11)*[1,1,1],            (30/11)*[1,1]]);
    p = round(10^6*p)/10^6;
    z = [10/11, 5; 20/11, 3; 30/11, 2];
    z = z(3:-1:1,:);
    @printf("\n");
    @printf(" Coefficients are rounded up at the 6-th digits after\n");
    @printf("    decimal point. Originally, \n");
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