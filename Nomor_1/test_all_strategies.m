function test_all_strategies()
    matrices = generate_matrices();
    test_block_lu_strategy(matrices);
    test_lu_pivot_strategy(matrices);
    test_recur_lu_strategy(matrices);
end

% Function to generate and store matrices for reuse
function matrices = generate_matrices()
    Ns = [16, 64, 128, 256, 512, 1000];
    matrices = struct();  % Initialize a structure to store matrices
    
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 4), floor(N / 4); floor(N / 2), floor(N / 2)];
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            % Store each matrix in the structure
            matrices.(['N' num2str(N)]).(['p' num2str(p) '_q' num2str(q)]) = create_banded_matrix(N, p, q);
        end
    end
end

% Update each test strategy function to take the matrices as an argument
function test_block_lu_strategy(matrices)
    Ns = [16, 64, 128, 256, 512, 1000];
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 4), floor(N / 4); floor(N / 2), floor(N / 2)];
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            % Retrieve the precomputed matrix
            A = matrices.(['N' num2str(N)]).(['p' num2str(p) '_q' num2str(q)]);
            b = rand(N, 1);
            printf("\nTesting Block LU - N = %d, p = %d, q = %d\n", N, p, q);
            [x, L, U, exec_time] = solve_lu_banded_block(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - L*U) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            printf("Execution Time: %.6f seconds\n", exec_time);
        end
    end
end

% Repeat similarly for Pivot LU and Recursive LU
function test_lu_pivot_strategy(matrices)
    Ns = [16, 64, 128, 256, 512, 1000];
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 4), floor(N / 4); floor(N / 2), floor(N / 2)];
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            A = matrices.(['N' num2str(N)]).(['p' num2str(p) '_q' num2str(q)]);
            b = rand(N, 1);
            printf("\nTesting Pivot LU - N = %d, p = %d, q = %d\n", N, p, q);
            [x, L, U, P, time_taken] = solve_lu_with_pivoting(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - L*U) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            printf("Execution Time: %.6f seconds\n", time_taken);
        end
    end
end

function test_recur_lu_strategy(matrices)
    Ns = [16, 64, 128, 256, 512, 1000];
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N / 4), floor(N / 4); floor(N / 2), floor(N / 2)];
        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            A = matrices.(['N' num2str(N)]).(['p' num2str(p) '_q' num2str(q)]);
            b = rand(N, 1);
            printf("\nTesting Recursive LU - N = %d, p = %d, q = %d\n", N, p, q);
            [x, L, U, exec_time] = solve_lu_banded_manual(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - L*U) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            printf("Execution Time: %.6f seconds\n", exec_time);
        end
    end
end

function A = create_banded_matrix(N, p, q)
    A = zeros(N);
    for i = 1:N
        for j = max(1, i - q):min(N, i + p)
            A(i, j) = rand();
        end
        % Strengthen diagonal elements to improve conditioning
        A(i, i) = A(i, i) + N;  
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

function [L, U, P] = lu_factorization_with_pivoting(A)
    N = size(A, 1);
    L = eye(N);
    U = A;
    P = (1:N)';  % Initialize permutation vector

    for k = 1:N-1
        % Pivoting
        [~, t] = max(abs(U(k:N, k)));
        t = t + k - 1;

        % Swap rows in U
        U([k, t], :) = U([t, k], :);
        
        % Swap rows in L (only above k)
        if k > 1
            L([k, t], 1:k-1) = L([t, k], 1:k-1);
        end

        % Update permutation vector
        P([k, t]) = P([t, k]);

        % Eliminate below the pivot
        for i = k+1:N
            L(i, k) = U(i, k) / U(k, k);
            U(i, :) = U(i, :) - L(i, k) * U(k, :);
        end
    end
end

function [L, U] = recursive_lu_factorization(A)
    n = size(A, 1);

    % Increase threshold to minimize recursion depth
    threshold = 1024;
    
    if n <= threshold
        [L, U, ~] = lu(A);
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
    
    % Use the built-in LU factorization for larger matrices to avoid excessive recursion
    if size(schur_complement, 1) <= threshold
        [L_schur, U_schur] = lu(schur_complement);
    else
        % Avoid further recursion on large matrices
        [L_schur, U_schur] = recursive_lu_factorization(schur_complement);
    end
    
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

% Solver Functions
function [x, L, U, exec_time] = solve_lu_banded_block(A, b)
    tic;
    [L, U] = block_lu_factorization(A);
    y = forward_substitution(L, b);
    x = backward_substitution(U, y);
    exec_time = toc;
end

function [x, L, U, P, time_taken] = solve_lu_with_pivoting(A, b)
    tic;
    [L, U, P] = lu_factorization_with_pivoting(A);
    b_permuted = b(P);
    y = forward_substitution(L, b_permuted);
    x = backward_substitution(U, y);
    time_taken = toc;
end

function [x, L, U, exec_time] = solve_lu_banded_manual(A, b)
    tic;
    [L, U] = recursive_lu_factorization(A);
    y = forward_substitution(L, b);
    x = backward_substitution(U, y);
    exec_time = toc;
end

% Run the tests
test_all_strategies();
