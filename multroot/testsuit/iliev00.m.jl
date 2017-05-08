function iliev00()
#
#  generalization of Iliev example:
#
#    (x-1) (x-2)^2 (x-3)^3
#  
   k = 0;
   p = poly([1,2,2,3,3,3]);
   z = [1,1; 2,2; 3,3];
   if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
   else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       p,z
endfunction iliev00()
#
#  generalization of Iliev example:
#
#    (x-1) (x-2)^2 (x-3)^3
#  
   k = 0;
   p = poly([1,2,2,3,3,3]);
   z = [1,1; 2,2; 3,3];
   if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
   else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       p,z
end