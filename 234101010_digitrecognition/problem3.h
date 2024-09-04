// Function for the 3rd problem's solution: Re-estimation of HMM parameters
void re_estimation() {
    for (int i = 1; i <= N; i++) {
        pi_bar[i] = gamma[1][i];
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            int mi = 0;
            long double maxValue = DBL_MIN;
            long double adjustSum = 0;
            for (int j = 1; j <= N; j++) {
                long double numerator = 0.0, denominator = 0.0;
                for (int t = 1; t <= T - 1; t++) {
                    numerator += xi[t][i][j];
                    denominator += gamma[t][i];
                }
                A_bar[i][j] = (numerator / denominator);
                if (A_bar[i][j] > maxValue) {
                    maxValue = A_bar[i][j];
                    mi = j;
                }
                adjustSum += A_bar[i][j];
            }
            A_bar[i][mi] += (1 - adjustSum);
        }
    }
    for (int j = 1; j <= N; j++) {
        int mi = 0;
        long double maxValue = DBL_MIN;
        long double adjustSum = 0;
        for (int k = 1; k <= M; k++) {
            long double numerator = 0.0, denominator = 0.0;
            for (int t = 1; t <= T; t++) {
                if (O[t] == k) {
                    numerator += gamma[t][j];
                }
                denominator += gamma[t][j];
            }
            B_bar[j][k] = (numerator / denominator);
            if (B_bar[j][k] > maxValue) {
                maxValue = B_bar[j][k];
                mi = k;
            }
            if (B_bar[j][k] < 1.00e-030) {
                B_bar[j][k] = 1.00e-030;
            }
            adjustSum += B_bar[j][k];
        }
        B_bar[j][mi] += (1 - adjustSum);
    }
    for (int i = 1; i <= N; i++) {
        pi[i] = pi_bar[i];
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            A[i][j] = A_bar[i][j];
        }
    }
    for (int j = 1; j <= N; j++) {
        for (int k = 1; k <= M; k++) {
            B[j][k] = B_bar[j][k];
        }
    }
}
