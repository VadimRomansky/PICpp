void progon(double* a, double* c, double* b, int N, double* f, double* x, double* alpha, double* beta){
	for(int i = 0; i <N; ++i){
		alpha[i] = 0;
		beta[i] = 0;
	}
	alpha[0] = -b[0]/c[0];
	beta[0] = f[0]/c[0];

	for(int j = 1; j < N; ++j){
		alpha[j] = (-b[j])/(a[j-1]*alpha[j-1] + c[j]);      
		beta[j] = (f[j] - a[j-1]*beta[j-1])/(a[j-1]*alpha[j-1] + c[j]);
	}

	x[N] = (f[N] - a[N-1]*beta[N-1])/(a[N-1]*alpha[N-1] + c[N]);
	for(int j = 0; j < N; ++j){
		x[N-1-j] = alpha[N-j-1]*x[N-j]+beta[N-j-1];
	}
}
