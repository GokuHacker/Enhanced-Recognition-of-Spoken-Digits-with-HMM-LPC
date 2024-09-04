// Function for the 2nd problem's solution: Viterbi algorithm
void viterbi_algo() {
    for (int i = 1; i <= N; i++) {
        delta[1][i] = pi[i] * B[i][O[1]];
        psi[1][i] = 0;
    }
    for (int t = 2; t <= T; t++) {
        for (int j = 1; j <= N; j++) {
            long double max = DBL_MIN;
            int index = 0;
            for (int i = 1; i <= N; i++) {
                if (delta[t - 1][i] * A[i][j] > max) {
                    max = delta[t - 1][i] * A[i][j];
                    index = i;
                }
            }
            delta[t][j] = max * B[j][O[t]];
            psi[t][j] = index;
        }
    }
    P_star = DBL_MIN;
    for (int i = 1; i <= N; i++) {
        if (delta[T][i] > P_star) {
            P_star = delta[T][i];
            q_star[T] = i;
        }
    }
    for (int t = T - 1; t >= 1; t--) {
        q_star[t] = psi[t + 1][q_star[t + 1]];
    }
    // Write observed sequence and most likely state sequence
    FILE *fp = fopen("predicted_seq_viterbi.txt", "w");
    for (int t = 1; t <= T; t++) {
        fprintf(fp, "%4d\t", O[t]);
    }
    fprintf(fp, "\n");
    for (int t = 1; t <= T; t++) {
        fprintf(fp, "%4d\t", q_star[t]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

// Function to calculate xi values
void calculate_xi() {
    for (int t = 1; t < T; t++) {
        long double denominator = 0.0;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                denominator += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
            }
        }
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j]) / denominator;
            }
        }
    }
}
