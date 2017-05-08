function igyp03()
#
#  generalization of Igarash and Ypma
# 
   p = poly([10*(1+i)*[1,1,1],1,i,2,2i,3,4i,5]);
   z = [10*(1+i),3; 1,1; i,1; 2,1; 2i,1; 3,1; 4i,1; 5,1];
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
endfunction igyp03()
#
#  generalization of Igarash and Ypma
# 
   p = poly([10*(1+i)*[1,1,1],1,i,2,2i,3,4i,5]);
   z = [10*(1+i),3; 1,1; i,1; 2,1; 2i,1; 3,1; 4i,1; 5,1];
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