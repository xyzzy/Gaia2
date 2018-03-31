/* 
 * This program merges the 5 separate title letters into a single frame. Frames 600 to 856.
 *
 * ffmpeg  -r 16 -f image2 -i  title-text/text_1_%05d.png -vcodec libx264 -crf 31 -pix_fmt yuv420p title-text.mp4 -y
 */

/*
 *	This file is part of Gaia2, Colour particle simulation
 *	Copyright (C) 2018, xyzzy@rockingship.net
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <getopt.h>
#include <gd.h>
#include <unistd.h>
#include <signal.h>

#define MAXW (1920*1)
#define MAXH (1080*1)

typedef struct
{
  int16_t r, g, b;
  int next;
} P;

P grid[MAXW*MAXH];
P mean[MAXW*MAXH];
int evalseq[MAXW*MAXH];
int numEval;
int rmap[256*256*256];
char *opt_mmapName;
int mmapfd = -1;
int xyok[MAXW*MAXH];
int root[256*256*256];
int version[MAXW*MAXH];
int gVersion;

char *opt_output;
char *opt_input;
char *arg_image;
int opt_framenr = 0;
int opt_verbose = 0;
int opt_timer = 1;
volatile int tick = 1;
uint64_t gCntCalcArea;

int todoIn, todoOut;
uint32_t todoRGB[256*256*256];
uint32_t todoVersion[256*256*256];

void sigAlarm(int) {
	if (opt_timer) {
		tick++;
		alarm(opt_timer);
	}
}

const char * timeAsString(void) {
	static char tstr[128];

	time_t t = time(0);
	struct tm *tm = localtime(&t);
	strftime(tstr, sizeof(tstr), "%F %T", tm);

	return tstr;
}

void save(int framenr)
{
	char fname[128];
	sprintf(fname, opt_output, framenr);
	FILE *fil = fopen(fname, "wb");
	if (fil == NULL) {
		fprintf(stderr, "Could not open output file\n");
		exit(1);
	}

	gdImagePtr im = gdImageCreateTrueColor(MAXW, MAXH);

	for(int y=0; y<MAXH; y++) {
		for (int x=0; x<MAXW; x++) {

			int v = grid[y*MAXW+x].r<<16 | grid[y*MAXW+x].g<<8 | grid[y*MAXW+x].b;

			int r = ((v>>16) & 0xFF);
			int g = ((v>>8) & 0xFF);
			int b = (v & 0xFF);

			int c = gdImageColorAllocate(im, r, g, b);
			gdImageSetPixel(im, x, y, c);

		}
	}
	gdImagePng(im, fil);
	gdImageDestroy(im);
	fclose(fil);
}

int load(const char *fname)
{
gdImagePtr im = NULL;
FILE *fil;
unsigned char c[2];

		/* open source */
		fil = fopen(fname, "rb");

		if (fil == NULL) {
			fprintf(stderr, "Could not open %s: %m\n", fname);
			return -1;
		}
		if (fread(c, 2, 1, fil) == 1) {
			rewind(fil);
			if (c[0] == 0x89 && c[1] == 0x50)
				im = gdImageCreateFromPng(fil);
			if (c[0] == 0x47 && c[1] == 0x49)
				im = gdImageCreateFromGif(fil);
			if (c[0] == 0xff && c[1] == 0xd8)
				im = gdImageCreateFromJpeg(fil);
		}
		if (im == NULL) {
			fprintf(stderr, "Could not load source image %x %x\n", c[0], c[1]);
			return -1;
		}
		fclose(fil);

		for (int y=0; y<MAXH; y++) {
		for (int x=0; x<MAXW; x++) {
			uint32_t v = gdImageGetTrueColorPixel(im, x, y);

			unsigned r = ((v >> 16) & 0xFF);
			unsigned g = ((v >> 8) & 0xFF);
			unsigned b = (v & 0xFF);

			// save pixel
			assert(r < 256);
			assert(g < 256);
			assert(b < 256);

			if (r+g+b != 0) {
				int yx = y*MAXW+x;
				grid[yx].r = r;
				grid[yx].g = g;
				grid[yx].b = b;
			}
		}}

	return 0;
}

void usage(const char *argv0, int verbose)
{
	printf("Usage: %s [options] <output image>\n", argv0);
	if (verbose)
		return;

	printf("\n\
options:\n\
	--startframe=n	starting frame number\
	-v --verbose	show progress\n");
}

int main(int argc, char *argv[])
{
gdImagePtr im = NULL;
FILE *fil;
unsigned char c[2];

	// randomize
	srand(time(NULL));

	for (;;) {
		int option_index = 0;
		enum {	LO_HELP=1, LO_VERBOSE='v', LO_STARTFRAME};
		static struct option long_options[] = {
			/* name, has_arg, flag, val */
			{"help",       0, 0, LO_HELP},
			{"verbose",    0, 0, LO_VERBOSE},
			{"startframe", 1, 0, LO_STARTFRAME},
			{NULL,         0, 0, 0}
		};

		char optstring[128], *cp;
		cp = optstring;
		for (int i=0; long_options[i].name; i++) {
			if (isalpha(long_options[i].val)) {
				*cp++ = long_options[i].val;
				if (long_options[i].has_arg)
					*cp++ = ':';
			}
		}
		*cp++ = '\0';

		int c = getopt_long (argc, argv, optstring, long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case LO_VERBOSE:
			opt_verbose++;
			break;
		case LO_STARTFRAME:
			opt_framenr = strtol(optarg, NULL, 10);
			break;
		case LO_HELP:
			usage(argv[0], 0);
			exit(0);
			break;
		case '?':
			fprintf(stderr,"Try `%s --help' for more information.\n", argv[0]);
			exit(1);
			break;
		default:
			fprintf (stderr, "getopt returned character code %d\n", c);
			exit(1);
		}
	 }

	if (argc-optind < 2) {
		usage(argv[0], 1);
		exit(1);
	}
	opt_output = argv[optind++];
	opt_input = argv[optind++];
	if (optind < argc)
		arg_image = argv[optind++];

	if (opt_timer) {
		signal(SIGALRM, sigAlarm);
		alarm(opt_timer);
	}

	int inframe = 800;
	for (int inframe=600; inframe<856; inframe++) {
		char fname[1024];

		if (tick) {
			fprintf(stderr, "\r\e[K[%s] %.2f%%", timeAsString(), (100.0 * (inframe-600)) / 256);
			tick = 0;
		}

		for (int y=0; y<MAXH; y++) {
			for (int x=0; x<MAXW; x++) {
				int yx = y*MAXW+x;
				grid[yx].r = grid[yx].g = grid[yx].b = 0;
			}
		}

		sprintf(fname, "title-1/y-%05d.png", inframe);
		load(fname);
		sprintf(fname, "title-2/y-%05d.png", inframe);
		load(fname);
		sprintf(fname, "title-3/y-%05d.png", inframe);
		load(fname);
		sprintf(fname, "title-4/y-%05d.png", inframe);
		load(fname);
		sprintf(fname, "title-5/y-%05d.png", inframe);
		load(fname);
		
		save(opt_framenr++);
	}
	fprintf(stderr, "\r\e[K[%s]\n", timeAsString());

	return 0;
}
