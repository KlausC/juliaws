function fib05()
#
#  test polynomial suggested by Goedecker
#
    n = 5;
    p = [-1,ones(1,n)];
    z = [-.6783507129699967-.4585361872731445*i,  1;          -.6783507129699967+.4585361872731445*i,  1;          .19537659464725405-.8488536405462456*i,  1;          .19537659464725405+.8488536405462456*i,  1;          1.9659482366454853,                      1];
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;p,z
endfunction fib05()
#
#  test polynomial suggested by Goedecker
#
    n = 5;
    p = [-1,ones(1,n)];
    z = [-.6783507129699967-.4585361872731445*i,  1;          -.6783507129699967+.4585361872731445*i,  1;          .19537659464725405-.8488536405462456*i,  1;          .19537659464725405+.8488536405462456*i,  1;          1.9659482366454853,                      1];
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;p,z
end