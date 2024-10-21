function ax = skew_matrix(a)
    ax = [0 -a(3) a(2);
          a(3) 0 -a(1);
          -a(2) a(1) 0];
end