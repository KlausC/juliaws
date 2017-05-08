function triple01()
#
#  test polynomial suggested by Goedecker
#
    p = triple(5,5,5);
    z = [[0.9,1,1.1]',[5,5,5]'];
    @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');p,z
end