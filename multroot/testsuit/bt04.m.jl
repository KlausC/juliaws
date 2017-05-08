function bt04()
#
# Brugnano and Trigiante
#
    p = poly([1,1,1, -1*[1,1,1,1], (.5+i)*[1,1,1], (.5-i)*[1,1,1],             .5*(1+i)*[1,1], .5*(1-i)*[1,1]]);
    z = [1,3; -1,4; .5+i,3; .5-i,3; .5*(1+i),2; .5*(1-i),2];
    
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
endfunction bt04()
#
# Brugnano and Trigiante
#
    p = poly([1,1,1, -1*[1,1,1,1], (.5+i)*[1,1,1], (.5-i)*[1,1,1],             .5*(1+i)*[1,1], .5*(1-i)*[1,1]]);
    z = [1,3; -1,4; .5+i,3; .5-i,3; .5*(1+i),2; .5*(1-i),2];
    
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