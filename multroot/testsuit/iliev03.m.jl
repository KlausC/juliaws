function iliev03()
#
#  generalization of Iliev example:
#
#    (x-1) (x-2)^2 (x-3)^3
#  
   k = 8;
   l = 1*k; m = 2*k; n = 3*k;
   p = poly([ones(1,l),2*ones(1,m),3*ones(1,n)]);
   z = [1,l; 2,m; 3,n];
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
endfunction iliev03()
#
#  generalization of Iliev example:
#
#    (x-1) (x-2)^2 (x-3)^3
#  
   k = 8;
   l = 1*k; m = 2*k; n = 3*k;
   p = poly([ones(1,l),2*ones(1,m),3*ones(1,n)]);
   z = [1,l; 2,m; 3,n];
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