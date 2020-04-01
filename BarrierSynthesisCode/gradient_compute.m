function [h, x_d, y_d] = gradient_compute(x_state, dfx, dfy, g)

    h = double(g(x_state(1), x_state(2)));                                  % Can be verified with the Value in D-
    x_d = double(dfx(x_state(1), x_state(2)));                              % Can be verified with vale of [U,V] in fig(6)
    y_d = double(dfy(x_state(1), x_state(2)));                              % Can be verified with vale of [U,V] in fig(6)

end
