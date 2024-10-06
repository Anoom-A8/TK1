function test_block_lu_strategy()
    Ns = [16, 64, 128, 256, 512, 1000];
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 4), floor(N / 4); floor(N / 2), floor(N / 2)];
        
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            A = create_banded_matrix(N, p, q);
            b = rand(N, 1);
            
            printf("\nTesting N = %d, p = %d, q = %d\n", N, p, q);
            [x, L, U, exec_time] = solve_lu_banded_block(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - L*U) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            printf("Execution Time: %.6f seconds\n", exec_time);
        end
    end
end

function [L, U] = block_lu_factorization(A)
    n = size(A, 1);
    if n == 1
        L = eye(1);
        U = A;
        return;
    end

    mid = floor(n / 2);

    A11 = A(1:mid, 1:mid);
    A12 = A(1:mid, mid+1:end);
    A21 = A(mid+1:end, 1:mid);
    A22 = A(mid+1:end, mid+1:end);

    [L11, U11] = block_lu_factorization(A11);

    U12 = L11 \ A12;
    L21 = A21 / U11;

    S = A22 - L21 * U12;

    [L22, U22] = block_lu_factorization(S);

    L = [L11, zeros(size(L11, 1), size(L22, 2)); L21, L22];
    U = [U11, U12; zeros(size(U22, 1), size(U11, 2)), U22];
end

function x = forward_substitution(L, b)
    N = length(b);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i);
    end
end

function x = backward_substitution(U, b)
    N = length(b);
    x = zeros(N, 1);
    for i = N:-1:1
        x(i) = (b(i) - U(i, i+1:N) * x(i+1:N)) / U(i, i);
    end
end

function [x, L, U, exec_time] = solve_lu_banded_block(A, b)
    tic;
    
    [L, U] = block_lu_factorization(A);
    
    y = forward_substitution(L, b);
    x = backward_substitution(U, y);
    
    exec_time = toc;
end

function A = create_banded_matrix(N, p, q)
    A = zeros(N);
    for i = 1:N
        for j = max(1, i - q):min(N, i + p)
            A(i, j) = rand();
        end
    end
end

% Run the test
test_block_lu_strategy();
