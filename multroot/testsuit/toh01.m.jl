function toh01()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
   k = [-10:9]';
   z = 2*(k+0.5)/19 + i*sin(2*pi*(k+0.5)/19);
   p = poly(z);
   z = [z,ones(20,1)];
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
endfunction toh01()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
   k = [-10:9]';
   z = 2*(k+0.5)/19 + i*sin(2*pi*(k+0.5)/19);
   p = poly(z);
   z = [z,ones(20,1)];
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