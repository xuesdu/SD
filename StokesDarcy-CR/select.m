function out = select(cc,ii)
a.end = 0;
a.pic = 1;

if cc == 1 
    a.T = 0.5;
    switch ii
        case 0
            a.dir = 'ex1'; a.m = 2; a.n = 2;
        case 1
            a.dir = 'ex1'; a.m = 4; a.n = 4;
        case 2
            a.dir = 'ex1'; a.m = 8; a.n = 8;
        case 3
            a.dir = 'ex1'; a.m = 16; a.n = 16;
        case 4
            a.dir = 'ex1'; a.m = 32; a.n = 32;
        case 5 
            a.dir = 'ex1'; a.m = 64; a.n = 64;
        case 6
            a.dir = 'ex1'; a.m = 128; a.n = 128;
        case 7
            a.dir = 'ex1'; a.m = 256; a.n = 256;
        otherwise
            a.end = 1;
    end
elseif cc == 2
    switch ii
        case 0
            a.dir = 'ex2'; a.m = 2; a.n = 2;
        case 1
            a.dir = 'ex2'; a.m = 4; a.n = 4;
        case 2
            a.dir = 'ex2'; a.m = 8; a.n = 8;
        case 3
            a.dir = 'ex2'; a.m = 16; a.n = 16;
        case 4
            a.dir = 'ex2'; a.m = 32; a.n = 32;
        case 5 
            a.dir = 'ex2'; a.m = 64; a.n = 64;
        case 6
            a.dir = 'ex2'; a.m = 128; a.n = 128;
        case 7
            a.dir = 'ex2'; a.m = 256; a.n = 256;
        otherwise
            a.end = 1;
    end
elseif cc == 3
    switch ii
        case 0
            a.dir = 'ex3'; a.m = 2; a.n = 2;
        case 1
            a.dir = 'ex3'; a.m = 4; a.n = 4;
        case 2
            a.dir = 'ex3'; a.m = 8; a.n = 8;
        case 3
            a.dir = 'ex3'; a.m = 16; a.n = 16;
        case 4
            a.dir = 'ex3'; a.m = 32; a.n = 32;
        case 5 
            a.dir = 'ex3'; a.m = 64; a.n = 64;
        case 6
            a.dir = 'ex3'; a.m = 256; a.n = 256;
        case 7
            a.dir = 'ex3'; a.m = 512; a.n = 512;
        otherwise
            a.end = 1;
    end
elseif cc == 4
    switch ii
        case 0
            a.dir = 'ex4'; a.m = 4; a.n = 4;
        case 1
            a.dir = 'ex4'; a.m = 8; a.n = 8;
        case 2
            a.dir = 'ex4'; a.m = 16; a.n = 16;
        case 3
            a.dir = 'ex4'; a.m = 32; a.n = 32;
        case 4
            a.dir = 'ex4'; a.m = 64; a.n = 64;
        case 5 
            a.dir = 'ex4'; a.m = 128; a.n = 128;
        case 6
            a.dir = 'ex4'; a.m = 256; a.n = 256;
        case 7
            a.dir = 'ex4'; a.m = 512; a.n = 512;
        otherwise
            a.end = 1;
    end
elseif cc == 5
    switch ii
        case 0
            a.dir = 'ex5'; a.m = 4; a.n = 4;
        case 1
            a.dir = 'ex5'; a.m = 8; a.n = 8;
        case 2
            a.dir = 'ex5'; a.m = 16; a.n = 16;
        case 3
            a.dir = 'ex5'; a.m = 32; a.n = 32;
        case 4
            a.dir = 'ex5'; a.m = 64; a.n = 64;
        case 5 
            a.dir = 'ex5'; a.m = 128; a.n = 128;
        case 6
            a.dir = 'ex5'; a.m = 256; a.n = 256;
        case 7
            a.dir = 'ex5'; a.m = 512; a.n = 512;
        otherwise
            a.end = 1;
    end
elseif cc == 6
    switch ii
        case 0
            a.dir = 'ex6'; a.m = 2; a.n = 2;
        case 1
            a.dir = 'ex6'; a.m = 4; a.n = 4;
        case 2
            a.dir = 'ex6'; a.m = 8; a.n = 8;
        case 3
            a.dir = 'ex6'; a.m = 16; a.n = 16;
        case 4
            a.dir = 'ex6'; a.m = 32; a.n = 32;
        case 5 
            a.dir = 'ex6'; a.m = 128; a.n = 128;
        case 6
            a.dir = 'ex6'; a.m = 256; a.n = 256;
        case 7
            a.dir = 'ex6'; a.m = 512; a.n = 512;
        otherwise
            a.end = 1;
    end
end
if a.n > 60
    a.pic = 0;
end
out = a;