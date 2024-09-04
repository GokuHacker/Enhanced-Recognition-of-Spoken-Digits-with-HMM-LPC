//Normalize the data
void normalize_data(char file[100]){
	// Open input file
	FILE* input_file = fopen(file, "r");
	if (input_file == NULL) {
		printf("Error in Opening File!\n");
		return;
	}

	int amplitude = 0, average_amplitude = 0;
	int i = 0;
	int sample_count = 0;
	int min_amplitude = INT_MAX;
	int max_amplitude = INT_MIN;
	while (!feof(input_file)) {
		fscanf(input_file, "%d", &amplitude);
		average_amplitude += amplitude;
		min_amplitude = (amplitude < min_amplitude) ? amplitude : min_amplitude;
		max_amplitude = (amplitude > max_amplitude) ? amplitude : max_amplitude;
		sample_count++;
	}

	average_amplitude /= sample_count;
	T = (sample_count - FS) / 80 + 1;

	if (T > T_) {
		T = T_;
	}

	min_amplitude -= average_amplitude;
	max_amplitude -= average_amplitude;

	fseek(input_file, 0, SEEK_SET);

	while (!feof(input_file)) {
		fscanf(input_file, "%d", &amplitude);

		if (min_amplitude == max_amplitude) {
			amplitude = 0;
		} else {
			amplitude -= average_amplitude;
			amplitude = (amplitude * CLIP) / ((max_amplitude > min_amplitude) ? max_amplitude : (-1) * min_amplitude);
			// Store normalized data
			samples[i++] = amplitude;
		}
	}

	fclose(input_file);

}

void calculate_energy_of_frame(int frame_no) {
    int sample_start_index = frame_no * 80;
    energy[frame_no] = 0;
    int i = 0;
    while (i < FS) {
        energy[frame_no] += samples[i + sample_start_index] * samples[i + sample_start_index];
        energy[frame_no] /= FS;
        i++;
    }
}

long double calculate_max_energy() {
    int nf = T;
    long double max_energy = DBL_MIN;
    int f = 0;
    while (f < nf) {
        if (energy[f] > max_energy) {
            max_energy = energy[f];
        }
        f++;
    }
    return max_energy;
}

long double calculate_avg_energy(){
	int nf=T;
	long double avg_energy=0.0;
	int f=0;
	while(f<nf){
		avg_energy+=energy[f];
		f++;
	}
	return avg_energy/nf;
}

//mark starting and ending of speech activity
void mark_checkpoints(){
	int nf = T;
    int f = 0;
    while (f < nf) {
        calculate_energy_of_frame(f);
        f++;
    }
	long double threshold_energy = calculate_avg_energy() / 10;
    int isAboveThresholdStart = 1;
    int isAboveThresholdEnd = 1;
    start_frame = 0;
    end_frame = nf - 1;

	f = 0;
    while (f < nf - 5) {
        int i = 0;
        while (i < 5) {
            isAboveThresholdStart *= (energy[f + i] > threshold_energy);
            i++;
        }
        if (isAboveThresholdStart) {
            start_frame = ((f - 5) > 0) ? (f - 5) : 0;
            break;
        }
        isAboveThresholdStart = 1;
        f++;
    }
	f = nf - 1;
    while (f > 4) {
        int i = 0;
        while (i < 5) {
            isAboveThresholdEnd *= (energy[f - i] > threshold_energy);
            i++;
        }
        if (isAboveThresholdEnd) {
            end_frame = ((f + 5) < nf) ? (f + 5) : (nf - 1);
            break;
        }
        isAboveThresholdEnd = 1;
        f--;
    }
}