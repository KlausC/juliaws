function jt01a()
#
#  test polynomial suggested by Jenkins and Traub
#
   a = 10^(10);
   p = [1,-1,-a^2,a^2];
   z = [a,1;-a,1;1,1];
   if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.5f \t \t \t %3g \n", z');
   else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       
    
    
   p,z
endfunction jt01a()
#
#  test polynomial suggested by Jenkins and Traub
#
   a = 10^(10);
   p = [1,-1,-a^2,a^2];
   z = [a,1;-a,1;1,1];
   if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.5f \t \t \t %3g \n", z');
   else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       
    
    
   p,z
end