// Fits1.cpp : Defines the entry point for the console application.
//

#include "../cfitsio/include/fitsio.h"
#include <cstdio>

int main() {
	const int Nx = 4001;
	const int Ny = 101;
	double** image = new double*[Nx];
	for (int i = 0; i < Nx; ++i) {
		image[i] = new double[Ny];
	}

	FILE* input = fopen("../B.dat", "r");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			double a;
			fscanf(input, "%lf", &a);
			image[i][j] = a;
		}
	}
	fclose(input);

	fitsfile* fptr; /* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements, exposure;
	long naxes[2] = {300, 200}; /* image is 300 pixels wide by 200 rows */
	short array[200][300];
	status = 0; /* initialize status before calling fitsio routines */
	fits_create_file(&fptr, "../B.fits", &status); /* create new file */
	/* Create the primary array image (16-bit short integer pixels */
	fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
	/* Write a keyword; must pass the ADDRESS of the value */
	exposure = 1500.;
	fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);
	/* Initialize the values in the image with a linear ramp function */
	for (jj = 0; jj < naxes[1]; jj++)
		for (ii = 0; ii < naxes[0]; ii++)
			array[jj][ii] = ii + jj;

	nelements = naxes[0] * naxes[1]; /* number of pixels to write*/
	/* Write the array of integers to the image */
	fits_write_img(fptr, TSHORT, fpixel, nelements, array[0], &status);
	fits_close_file(fptr, &status); /* close the file */
	fits_report_error(stderr, status); /* print out any errormessages */
	return (status);

	for (int i = 0; i < Nx; ++i) {
		delete[] image[i];
	}
	delete[] image;
	return 0;
}

