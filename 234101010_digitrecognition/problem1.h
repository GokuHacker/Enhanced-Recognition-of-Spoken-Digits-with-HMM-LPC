// Function for the forward procedure (alpha calculation)
void alphaCalc() {
    for (int i = 1; i <= N; i++) {
        alpha[1][i] = pi[i] * B[i][O[1]];
    }
    for (int t = 1; t < T; t++) {
        for (int j = 1; j <= N; j++) {
            long double sum = 0;
            for (int i = 1; i <= N; i++) {
                sum += alpha[t][i] * A[i][j];
            }
            alpha[t + 1][j] = sum * B[j][O[t + 1]];
        }
    }
    // Write the alpha values to a file
    FILE *fp = fopen("alpha.txt", "w");
    for (int t = 1; t <= T; t++) {
        for (int j = 1; j <= N; j++) {
            fprintf(fp, "%e\t", alpha[t][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// Function to calculate the model's probability (1st problem's solution)
long double calculate_score() {
    long double probability = 0;
    for (int i = 1; i <= N; i++) {
        probability += alpha[T][i];
    }
    return probability;
}

// Function for the backward procedure (beta calculation)
void calculate_beta() {
    for (int i = 1; i <= N; i++) {
        beta[T][i] = 1;
    }
    for (int t = T - 1; t >= 1; t--) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
            }
        }
    }
    // Write the beta values to a file
    FILE *fp = fopen("beta.txt", "w");
    for (int t = 1; t < T; t++) {
        for (int j = 1; j <= N; j++) {
            fprintf(fp, "%e\t", beta[t][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// Function to predict the most likely state sequence using gamma values
void predict_state_sequence() {
    for (int t = 1; t <= T; t++) {
        long double max = 0;
        int index = 0;
        for (int j = 1; j <= N; j++) {
            if (gamma[t][j] > max) {
                max = gamma[t][j];
                index = j;
            }
        }
        q[t] = index;
    }
    // Write observed sequence and predicted sequence
    FILE *fp = fopen("predicted_seq_gamma.txt", "w");
    for (int t = 1; t <= T; t++) {
        fprintf(fp, "%4d\t", O[t]);
    }
    fprintf(fp, "\n");
    for (int t = 1; t <= T; t++) {
        fprintf(fp, "%4d\t", q[t]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

// Function to calculate gamma values
void calculate_gamma() {
    for (int t = 1; t <= T; t++) {
        long double sum = 0;
        for (int i = 1; i <= N; i++) {
            sum += alpha[t][i] * beta[t][i];
        }
        for (int i = 1; i <= N; i++) {
            gamma[t][i] = alpha[t][i] * beta[t][i] / sum;
        }
    }
    // Write it to a file
    FILE *fp = fopen("gamma.txt", "w");
    for (int t = 1; t <= T; t++) {
        for (int j = 1; j <= N; j++) {
            fprintf(fp, "%.16e\t", gamma[t][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    predict_state_sequence();
}
