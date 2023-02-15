import subprocess

#"python3 3_mat_least_square.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 10000 -tao_ls_type more-thuente -m 'motion_mesh1.msh' -o 'test3' -er 1.0 -es 1.0e-2 -lr 2.0 -ls 0.125 -vr 0.3 -vs 0.5 -k 5.0e-3 -e 4.0e-3 -p 2.0 -q 1.0"]

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 motion.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 5000 -tao_ls_type unit -m 'motion.msh' -o 'test1' -er 1.0 -es 1.0e-2 -lr 2.5 -ls 0.125 -vr 0.3 -vs 0.5 -k 5.0e-3 -e 4.0e-3 -p 2.0 -q 2.0"]



i = 1
for program in program_list:
    print("------------------------------------------------------------------------------")
    print("")
    print("Running test #{}".format(i))
    print("")
    print(program)
    print("")
    subprocess.run(program, shell = True)
    i = i + 1
