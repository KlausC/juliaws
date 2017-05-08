function jt04()
#
#  test polynomial suggested by Jenkins and Traub
#
   p = poly([0.1,0.1,0.1,0.5,0.6,0.7]);
   z = [0.5, 1; 0.6, 1; 0.7, 1; 0.1, 3];
   @printf("\n");
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
endfunction jt04()
#
#  test polynomial suggested by Jenkins and Traub
#
   p = poly([0.1,0.1,0.1,0.5,0.6,0.7]);
   z = [0.5, 1; 0.6, 1; 0.7, 1; 0.1, 3];
   @printf("\n");
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