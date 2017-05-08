function jt06()
#
#  test polynomial suggested by Jenkins and Traub
#
   p = poly([.1,1.001,.998,1.00002,.99999]);
   z = [.1,1.001,.998,1.00002,.99999];
   z = [z',ones(5,1)];
   @printf("\n");
   @printf(" Ill-conditioned polynomial \n");
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
endfunction jt06()
#
#  test polynomial suggested by Jenkins and Traub
#
   p = poly([.1,1.001,.998,1.00002,.99999]);
   z = [.1,1.001,.998,1.00002,.99999];
   z = [z',ones(5,1)];
   @printf("\n");
   @printf(" Ill-conditioned polynomial \n");
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