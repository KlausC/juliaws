function miyak02()
#
#  test polynomial suggested by Goedecker
#  square of Fibocacci polynomial
#
    p = poly([(1.1+1.1*i)*ones(1,4),(3.2+2.3*i)*[1,1],2.1+1.5*i]);
    p = conv(p,p);
    z = [[(1.1+1.1*i),(3.2+2.3*i),2.1+1.5*i]',2*[4,2,1]'];
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    p,z
end