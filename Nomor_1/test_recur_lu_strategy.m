function test_recur_lu_strategy()
    Ns = [16, 64, 128, 256, 512, 1000];
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 2), floor(N / 2)];
        
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            A = create_banded_matrix(N, p, q);
            b = rand(N, 1);
            
            printf("\nTesting N = %d, p = %d, q = %d\n", N, p, q);
            [x, L, U, exec_time] = solve_lu_banded_manual(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - L*U) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            
            printf("Execution Time: %.6f seconds\n", exec_time);
        end
    end
end

function [L, U] = recursive_lu_factorization(A)
    n = size(A, 1);
    
    if n == 1
        % Base case for 1x1 matrix
        L = 1;
        U = A;
        return;
    end
    
    % Partition the matrix A
    a11 = A(1, 1);
    bT = A(1, 2:end);
    c = A(2:end, 1);
    D = A(2:end, 2:end);
    
    if a11 != 0
        c_div_a11 = c / a11;
        schur_complement = D - c_div_a11 * bT;
    else
        error("Matrix is singular, cannot proceed with LU factorization");
    end
    
    % Recursively factorize the Schur complement
    [L_schur, U_schur] = recursive_lu_factorization(schur_complement);
    
    % Construct the full L and U matrices
    L = [1, zeros(1, n - 1); c_div_a11, L_schur];
    U = [a11, bT; zeros(n - 1, 1), U_schur];
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

function [x, L, U, exec_time] = solve_lu_banded_manual(A, b)
    tic;
    
    % Perform recursive LU decomposition
    [L, U] = recursive_lu_factorization(A);
    
    % Forward substitution to solve Ly = b
    y = forward_substitution(L, b);
    
    % Back substitution to solve Ux = y
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
test_recursive_lu_banded();