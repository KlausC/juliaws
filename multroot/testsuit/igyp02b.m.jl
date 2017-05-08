function igyp02b()
#
#  generalization of Igarash and Ypma
# 
   m = 7;
   p = poly([10*(1+i)*ones(1,m),-1*ones(1,10-m)]);
   z = [10*(1+i),m; -1, 10-m];
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
endfunction igyp02b()
#
#  generalization of Igarash and Ypma
# 
   m = 7;
   p = poly([10*(1+i)*ones(1,m),-1*ones(1,10-m)]);
   z = [10*(1+i),m; -1, 10-m];
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