void record_digit()
{
	int choice;
	do
	{
	system("Live_testing\\Recording_Module.exe 3 Live_testing\\input_file.wav Live_testing\\input_file.txt");
	printf("Playing sound:\n");
	PlaySound(TEXT("Live_testing\\input_file.wav"),NULL,SND_SYNC);
	printf("Press 1:Proceed 2:Re-record digit\n");
	scanf("%d",&choice);
	printf("Recording completed\n");
	}while(choice==2);
	printf("record function done\n");
}
//function to find the short term energy of the given frame
lld short_term_energy(int l, int r,lld arr[])
{
	lld energy = 0;
	for(int i = l; i<=r; i++)
		energy += arr[i]*arr[i];
	return energy/(r-l+1);
}

//function to trim the live testing data and extract the useful information from it
void trim_file()
{
	FILE *fptr;
	fptr=fopen("Live_testing\\input_file.txt","r");
	if(fptr==NULL)
	{
		printf("Error opening file\n");
		return;
	}
	int n = 0;
	lld arr[100000],a;
	while(fscanf(fptr,"%Lf",&a) != EOF)
		arr[++n] = a;
	fclose(fptr);
	lld threshold_energy;
	threshold_energy= short_term_energy(1,320,arr);	//function call to find the threshold energy of the silence part
	printf("threshold done\n");
	int l = 0,r = n,s,f = 0;
	lld frame_energy[1000],max=0;
	//find the frame with largest short term energy
	//We move left of this frame to find the starting point of useful data and move right to find the ending point of useful data
	for(int i = 1; i+320<=n; i += 320)
	{
		frame_energy[++f] = short_term_energy(i,i+319,arr);
		if(max < frame_energy[f])
		{
			max = frame_energy[f];
			s = f;
		}
	}
	for(int i = s; i>0; i--)
	{
		if(frame_energy[i] < 3*threshold_energy)
		{
			l = i*320 + 1;
			break;
		}
	}
	for(int i = s; i<=f; i++)
	{
		if(frame_energy[i] < 3*threshold_energy)
		{
			r = (i-1)*320 + 1;
			break;
		}
	}
	fptr = fopen("Live_testing\\input_file.txt","w");
	if(fptr==NULL)
	{
		printf("Error opening file\n");
		return;
	}
	//printing the useful extracted data into corresponding text file
	for(int i = l; i<= r; i++)
		fprintf(fptr,"%Lf\n",arr[i]);
	fclose(fptr);
}
int find_digit()
{
	lld probability=0,calc_probability=0;
	int index=0;
	for(int i=0;i<=9;i++)
	{
		read_lambda_values(i);
		calc_probability=forward_procedure();
		if(calc_probability>probability)
		{
			probability=calc_probability;
			index=i;
		}
	}
	return index;
}
void predict_digit()
{
	char filename[40];
	lld arr[60000],frame_arr[framesize+1];
	int index,framecount=0,predicted_digit;
	lld count;
	lld R[p+1],A[p+1],C[p+1];
	printf("Codebook loading starts.... \n");
	load_codebook();
	printf("Codebook loading eends.... \n");
	sprintf(filename,"Live_testing/input_file.txt");
	FILE *fp=fopen(filename,"r");
	
	if(fp==NULL)
	{
		printf("unable to access Live_testing file\n");
		return;
	}
	int k=0;
	//printf("file read starts.... \n");
	while(fscanf(fp,"%Lf",&arr[k])!=EOF)
		k++;
	fclose(fp);
	printf("recog starts\n");
	dc_shift(arr);
	normalisation(arr);
	framecount=0;
	for(int x=0;x+framesize<k;x+=80)
	{
		index=0;
		for(int y=x;y<x+framesize;y++)
		{
			frame_arr[index]=arr[y];
			index++;
		}
		
		hammingWindow(frame_arr);
		calculateRi(frame_arr,R);
		calculateAi(R,A);
		calculateCi(R,A,C);
		printf("ci done\n");
		raisedsinewindow(C);
		framecount++;
		O[framecount]=tokhura_distance_index(C);
	}
	T=framecount;
	predicted_digit=find_digit();
	printf("Predicted digit is: %d\n",predicted_digit);
}
void live_testing()
{
	record_digit();
	trim_file();
	predict_digit();
}