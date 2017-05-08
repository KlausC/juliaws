function bt03()
#
# Brugnano and Trigiante
#
    p = poly([i*ones(1,5),-i*ones(1,5),0.5i*ones(1,4),            -0.5i*ones(1,4),0.75i,-0.75i]);
    z = [i,5; -i, 5; 0.5i, 4; -0.5i, 4; .75i, 1; -0.75i, 1];
    
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
endfunction bt03()
#
# Brugnano and Trigiante
#
    p = poly([i*ones(1,5),-i*ones(1,5),0.5i*ones(1,4),            -0.5i*ones(1,4),0.75i,-0.75i]);
    z = [i,5; -i, 5; 0.5i, 4; -0.5i, 4; .75i, 1; -0.75i, 1];
    
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
endfunction bt03()
#
# Brugnano and Trigiante
#
    p = poly([i*ones(1,5),-i*ones(1,5),0.5i*ones(1,4),            -0.5i*ones(1,4),0.75i,-0.75i]);
    z = [i,5; -i, 5; 0.5i, 4; -0.5i, 4; .75i, 1; -0.75i, 1];
    
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