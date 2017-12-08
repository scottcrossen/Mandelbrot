#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <sys/time.h>
#include <string.h>


double current_time()
{
  struct timeval time;
  gettimeofday(&time, NULL);
  return ((double) time.tv_sec + (double) time.tv_usec * 1e-6);
}

int main(int argc, char* argv[])
{
  /* Parse the command line arguments. */
  if (argc != 8) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* The window in the plane. */
  const double xmin = atof(argv[1]);
  const double xmax = atof(argv[2]);
  const double ymin = atof(argv[3]);
  const double ymax = atof(argv[4]);

  /* Maximum number of iterations, at most 65535. */
  const uint16_t maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  const int xres = atoi(argv[6]);
  const int yres = (xres*(ymax-ymin))/(xmax-xmin);

  /* The output file name */
  const char* filename = argv[7];

  /* Open the file and write the header. */
  FILE * fp = fopen(filename,"wb");
  //char *comment="# Mandelbrot set";/* comment should start with # */

  /*write ASCII header to the file*/
  fprintf(fp,
          "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
          xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));

  /* Precompute pixel width and height. */
  double dx=(xmax-xmin)/xres;
  double dy=(ymax-ymin)/yres;

  double x, y; /* Coordinates of the current point in the complex plane. */
  //double u, v; /* Coordinates of the iterated point. */
  int i,j; /* Pixel counters */
  int k; /* Iteration counter */

  // Initialize results
  unsigned char*** result = (unsigned char***) malloc(sizeof(unsigned char**) * xres); // x column
  for (int i = 0; i < xres; i++) {
      result[i] = (unsigned char**) malloc(sizeof(unsigned char*) * yres); // y column
      for(int j = 0; j < yres; j++) {
          result[i][j] = (unsigned char*) malloc(sizeof(unsigned char) * 6); // pixel
      }
  }

  double start_time = current_time();
  #pragma omp parallel for private(i,j,k,y,x,u,v)
  for (j = 0; j < yres; j++) {
    y = ymax - j * dy;
    for(i = 0; i < xres; i++) {
      double u = 0.0;
      double v= 0.0;
      double u2 = u * u;
      double v2 = v*v;
      x = xmin + i * dx;
      /* iterate the point */
      for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
            v = 2 * u * v + y;
            u = u2 - v2 + x;
            u2 = u * u;
            v2 = v * v;
      };
      /* compute  pixel color and write it to file */
      if (k >= maxiter) {
        /* interior */
        const unsigned char black[] = {0, 0, 0, 0, 0, 0};
        memcpy (result[i][j], black, 6);
      }
      else {
        /* exterior */
        unsigned char color[6];
        color[0] = k >> 8;
        color[1] = k & 255;
        color[2] = k >> 8;
        color[3] = k & 255;
        color[4] = k >> 8;
        color[5] = k & 255;
        memcpy (result[i][j], color, 6);
      };
    }
  }
  double stop_time = current_time();
  printf("Elapsed time: %f\n", stop_time - start_time);

  // Write results to file
  for(int j = 0; j < yres; j++) {
      for (int i = 0; i < xres; i++) {
          fwrite(result[i][j], 6, 1, fp);
      }
  }

  // Free results
  for (int i = 0; i < xres; i++) {
      for(int j = 0; j < yres; j++) {
        free(result[i][j]);
      }
      free(result[i]);
  }
  free(result);

  fclose(fp);
  return 0;
}
