function triple03()
#
#  test polynomial suggested by Goedecker
#
    p = triple(18,10,16);
    z = [[0.9,1,1.1]',[18,10,16]'];
    @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');p,z
end