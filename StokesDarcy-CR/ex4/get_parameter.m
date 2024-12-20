function para = get_parameter()

para.box = getBox();
para.int_ps = int_ps();
para.int_pd = int_pd();

para.us1 = @fun_us1;
para.us2 = @fun_us2;
para.ps = @fun_ps;
para.us1_vec = @fun_us1_vec;
para.us2_vec = @fun_us2_vec;

para.us1_x = @fun_us1_x;
para.us1_y = @fun_us1_y;
para.us2_x = @fun_us2_x;
para.us2_y = @fun_us2_y;

para.ud1 = @fun_ud1;
para.ud2 = @fun_ud2;
para.pd = @fun_pd;

para.gs = @fun_gs;
para.fs1 = @fun_fs1;
para.fs2 = @fun_fs2;
para.gb1 = @fun_gb1;
para.gb2 = @fun_gb2;
para.gbp = @fun_gbp;

para.gd = @fun_gd;
para.fd1 = @fun_fd1;
para.fd2 = @fun_fd2;
para.eta1 = @fun_eta1;
para.eta2 = @fun_eta2;

para.nu = @fun_nu;
para.K = @fun_K;
para.alpha_BJS = @fun_alpha_BJS;
para.negativeone = @fun_negativeone;
para.one = @fun_one;