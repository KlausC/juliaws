
#include <julia.h>

int main(int argc, char *argv[]) {
    jl_init(JULIA_INIT_DIR);
    (void)jl_eval_string("println(sqrt(2.0))");
    jl_atexit_hook(0);
    return 0;
}
