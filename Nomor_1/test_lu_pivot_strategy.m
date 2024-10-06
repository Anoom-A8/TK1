function test_lu_pivot_strategy()
    Ns = [16, 64, 128, 256, 512, 1000];  % Set of matrix sizes
    for N = Ns
        p_q_pairs = [1, 1; 2, 1; 3, 4; floor(N/4), floor(N/4); floor(N/2), floor(N/2)];

        for i = 1:size(p_q_pairs, 1)
            p = p_q_pairs(i, 1);
            q = p_q_pairs(i, 2);
            A = create_banded_matrix(N, p, q);  % Generate the matrix based on N, p, q
            b = rand(N, 1);  % Random vector b for testing

            printf("\nTesting N = %d, p = %d, q = %d\n", N, p, q);
            
            % Solve using LU factorization with pivoting
            [xA, LA, UA, PA, timeA] = solve_lu_with_pivoting(A, b);
            printf("LU Factorization Residual Error: %.6e\n", norm(A - LA*UA) / norm(A));
            printf("Matrix Norm (Frobenius): %.4f\n", norm(A, 'fro'));
            printf("Condition Number: %.4f\n", cond(A));
            
            printf("\nExecution Time: %.6f seconds\n", timeA);
        end
    end
end

function A = create_banded_matrix(N, p, q)
    % Create an N x N matrix with bandwidth (p, q) based on the specified rules
    A = zeros(N);
    for i = 1:N
        for j = max(1, i-p):min(N, i+q)
            A(i, j) = rand();  % Populate within the band with random values
        end
    end
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

% Forward substitution (Ly = b)
function y = forward_substitution(L, b)
  N = length(b);
  y = zeros(N, 1);
  for i = 1:N
    y(i) = b(i) - L(i, 1:i-1) * y(1:i-1);
  end
end

% Backward substitution (Ux = y)
function x = backward_substitution(U, y)
  N = length(y);
  x = zeros(N, 1);
  for i = N:-1:1
    x(i) = (y(i) - U(i, i+1:N) * x(i+1:N)) / U(i, i);
  end
end

% Solve using LU factorization (Strategy A)
function [x, L, U, P, time_taken] = solve_lu_with_pivoting(A, b)
  tic;
  [L, U, P] = lu_factorization_with_pivoting(A);
  b_permuted = b(P);
  y = forward_substitution(L, b_permuted);
  x = backward_substitution(U, y);
  time_taken = toc;
end

% Run the test
test_all_strategies();