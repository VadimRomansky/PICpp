// Fits1.cpp : Defines the entry point for the console application.
//

#include "../../cfitsio/include/fitsio.h"
#include <cstdio>

int main() {
	const int Nx = 5001;
	const int Ny = 1001;
	/*double** image = new double*[Nx];
	for (int i = 0; i < Nx; ++i) {
		image[i] = new double[Ny];
	}*/
    printf("aa\n");
    double* image2 = new double[Nx*Ny];

	FILE* input = fopen("../B.dat", "r");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			double a;
			fscanf(input, "%lf", &a);
			image2[i*Ny + j] = a;
            char c;
            fscanf(input, "%c", &c);
		}
	}
	fclose(input);

	fitsfile* fptr; /* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements, exposure;
	long naxes[2] = {Ny, Nx}; /* image is 300 pixels wide by 200 rows */
	//int array[naxes[1]][Ny];
	status = 0; /* initialize status before calling fitsio routines */
	fits_create_file(&fptr, "!../B.fits", &status); /* create new file */
	/* Create the primary array image (16-bit short integer pixels */
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	/* Write a keyword; must pass the ADDRESS of the value */
	exposure = 1500.;
	fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);
	/* Initialize the values in the image with a linear ramp function */
	//for (jj = 0; jj < naxes[1]; jj++) {
        //for (ii = 0; ii < naxes[0]; ii++) {
            //array[jj][ii] = ii + jj;
            //image[jj][ii] = array[jj][ii];
            //image2[ii + jj*naxes[0]] = image[jj][ii];
        //}
    //}

	nelements = naxes[0] * naxes[1]; /* number of pixels to write*/
	/* Write the array of integers to the image */
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, image2, &status);
	fits_close_file(fptr, &status); /* close the file */
	fits_report_error(stderr, status); /* print out any errormessages */

	/*for (int i = 0; i < Nx; ++i) {
		delete[] image[i];
	}
	delete[] image;*/
    //delete[] image2;
    return (status);
}

