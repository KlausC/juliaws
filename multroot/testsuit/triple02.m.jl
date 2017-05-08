function triple02()
#
#  test polynomial suggested by Goedecker
#
    p = triple(10,10,10);
    z = [[0.9,1,1.1]',[10,10,10]'];
    @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');p,z
end