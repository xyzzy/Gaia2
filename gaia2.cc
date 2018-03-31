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
#include <sys/types.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <gd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include "jansson.h"

#if 0

#define MAXW (1920*2)
#define MAXH (1080*2)
#define NUMSLOT 240
//#define FIVExFIVE

#elif 1

#define MAXW (1920*1)
#define MAXH (1080*1)
#define NUMSLOT 230
//#define FIVExFIVE

#else

#define MAXW 320
#define MAXH 180

#endif

typedef struct 
{
	uint32_t	rgb;
	int		next;
} P;

P grid[MAXW*MAXH];
P mean[MAXW*MAXH];
int evalseq[MAXW*MAXH];
int numEval;
int rmap[256*256*256];
uint32_t (*nearestRGB) [256*256*256][NUMSLOT] = (uint32_t(*)[256*256*256][NUMSLOT])MAP_FAILED;
uint32_t (*nearestDist)[256*256*256][NUMSLOT];
char *opt_mmapName;
int mmapfd = -1;
int xyok[MAXW*MAXH];
int root[256*256*256];
int version[MAXW*MAXH];
int gVersion;

char *arg_output;
char *arg_json;
char *arg_image;
int opt_framenr = 0;
int opt_verbose = 0;
int opt_timer = 1;
volatile int tick = 1;
uint64_t gCntCalcArea;

int opt_lead_frames;
double opt_lead_fps_start;
double opt_lead_fps_end;


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
	char fname[1024];
	sprintf(fname, arg_output, framenr);
	FILE *fil = fopen(fname, "wb");
	if (fil == NULL) {
		fprintf(stderr, "Could not open output file\n");
		exit(1);
	}

	gdImagePtr im = gdImageCreateTrueColor(MAXW, MAXH);

	for(int y=0; y<MAXH; y++) {
		for (int x=0; x<MAXW; x++) {
			int yx = y*MAXW+x;

			int c = gdImageColorAllocate(im, (grid[yx].rgb>>16)&0xff, (grid[yx].rgb>>8)&0xff, (grid[yx].rgb>>0)&0xff);
			gdImageSetPixel(im, x, y, c);

		}
	}
	gdImagePng(im, fil);
	gdImageDestroy(im);
	fclose(fil);
}

static inline void calcmean(int x, int y)
{
	int r=0, g=0, b=0, n = 0;

	// weighted square without center pixel
	// 121
	// 242
	// 121

	int yx;

	if (x < 0 || y < 0 || x >= MAXW || y >= MAXH || !xyok[y*MAXW+x])
		return;

	if (version[y * MAXW + x] == gVersion)
		return;
	version[y * MAXW + x] = gVersion;

	gCntCalcArea++;

#ifdef FIVExFIVE
	/*
	 *  1  4  6  4  1
	 *  4 16 24 16  4
	 *  6 24 36 24  6
	 *  4 16 24 16  4
	 *  1  4  6  4  1
	 */

	yx = (y-2)*MAXW+(x-2);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=4; r += 4*(uint8_t)(grid[yx].rgb>>16); g += 4*(uint8_t)(grid[yx].rgb>>8); b += 4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=6; r += 6*(uint8_t)(grid[yx].rgb>>16); g += 6*(uint8_t)(grid[yx].rgb>>8); b += 6*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=4; r += 4*(uint8_t)(grid[yx].rgb>>16); g += 4*(uint8_t)(grid[yx].rgb>>8); b += 4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)(grid[yx].rgb); }

	yx = (y-1)*MAXW+(x-2);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 4; r +=  4*(uint8_t)(grid[yx].rgb>>16); g +=  4*(uint8_t)(grid[yx].rgb>>8); b +=  4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=16; r += 16*(uint8_t)(grid[yx].rgb>>16); g += 16*(uint8_t)(grid[yx].rgb>>8); b += 16*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=24; r += 24*(uint8_t)(grid[yx].rgb>>16); g += 24*(uint8_t)(grid[yx].rgb>>8); b += 24*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=16; r += 16*(uint8_t)(grid[yx].rgb>>16); g += 16*(uint8_t)(grid[yx].rgb>>8); b += 16*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 4; r +=  4*(uint8_t)(grid[yx].rgb>>16); g +=  4*(uint8_t)(grid[yx].rgb>>8); b +=  4*(uint8_t)(grid[yx].rgb); }

	yx = (y  )*MAXW+(x-2);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 6; r +=  6*(uint8_t)(grid[yx].rgb>>16); g +=  6*(uint8_t)(grid[yx].rgb>>8); b +=  6*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=24; r += 24*(uint8_t)(grid[yx].rgb>>16); g += 24*(uint8_t)(grid[yx].rgb>>8); b += 24*(uint8_t)(grid[yx].rgb); }
	yx++;
	// skip center with weitht 36
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=24; r += 24*(uint8_t)(grid[yx].rgb>>16); g += 24*(uint8_t)(grid[yx].rgb>>8); b += 24*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 6; r +=  6*(uint8_t)(grid[yx].rgb>>16); g +=  6*(uint8_t)(grid[yx].rgb>>8); b +=  6*(uint8_t)(grid[yx].rgb); }

	yx = (y+1)*MAXW+(x-2);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 4; r +=  4*(uint8_t)(grid[yx].rgb>>16); g +=  4*(uint8_t)(grid[yx].rgb>>8); b +=  4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=16; r += 16*(uint8_t)(grid[yx].rgb>>16); g += 16*(uint8_t)(grid[yx].rgb>>8); b += 16*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=24; r += 24*(uint8_t)(grid[yx].rgb>>16); g += 24*(uint8_t)(grid[yx].rgb>>8); b += 24*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=16; r += 16*(uint8_t)(grid[yx].rgb>>16); g += 16*(uint8_t)(grid[yx].rgb>>8); b += 16*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+= 4; r +=  4*(uint8_t)(grid[yx].rgb>>16); g +=  4*(uint8_t)(grid[yx].rgb>>8); b +=  4*(uint8_t)(grid[yx].rgb); }


	yx = (y+2)*MAXW+(x-2);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=4; r += 4*(uint8_t)(grid[yx].rgb>>16); g += 4*(uint8_t)(grid[yx].rgb>>8); b += 4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=6; r += 6*(uint8_t)(grid[yx].rgb>>16); g += 6*(uint8_t)(grid[yx].rgb>>8); b += 6*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=4; r += 4*(uint8_t)(grid[yx].rgb>>16); g += 4*(uint8_t)(grid[yx].rgb>>8); b += 4*(uint8_t)(grid[yx].rgb); }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)(grid[yx].rgb); }
#else
	/*
	 *  1 2 1
	 *  2 4 2
	 *  1 2 1
	 */

	yx = (y-1)*MAXW+(x-1);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)grid[yx].rgb; }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=2; r += 2*(uint8_t)(grid[yx].rgb>>16); g += 2*(uint8_t)(grid[yx].rgb>>8); b += 2*(uint8_t)grid[yx].rgb; }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)grid[yx].rgb; }

	yx = (y  )*MAXW+(x-1);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=2; r += 2*(uint8_t)(grid[yx].rgb>>16); g += 2*(uint8_t)(grid[yx].rgb>>8); b += 2*(uint8_t)grid[yx].rgb; }
	yx++;
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=2; r += 2*(uint8_t)(grid[yx].rgb>>16); g += 2*(uint8_t)(grid[yx].rgb>>8); b += 2*(uint8_t)grid[yx].rgb; }

	yx = (y+1)*MAXW+(x-1);
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)grid[yx].rgb; }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=2; r += 2*(uint8_t)(grid[yx].rgb>>16); g += 2*(uint8_t)(grid[yx].rgb>>8); b += 2*(uint8_t)grid[yx].rgb; }
	yx++;
	if (yx>=0 && yx< MAXW*MAXH && xyok[yx]) { n+=1; r += 1*(uint8_t)(grid[yx].rgb>>16); g += 1*(uint8_t)(grid[yx].rgb>>8); b += 1*(uint8_t)grid[yx].rgb; }
#endif

	// deliberate mix of float and ints
	//
	// the << floors the float

	yx = y * MAXW + x;
	uint32_t rgb = mean[yx].rgb;

	// unlink
	if (root[rgb] == yx) {
		root[rgb] = mean[yx].next;
	} else {
		for (int z=root[rgb]; z!=-1; z=mean[z].next) {
			if (mean[z].next == yx) {
				mean[z].next = mean[yx].next;
				break;
			}
		}
	}


	rgb = (r / n) << 16 | (g / n) << 8 | (b / n);
	rgb = (*nearestRGB)[rgb][0]; // find nearest colour

	mean[yx].rgb = rgb;
	mean[yx].next = root[rgb];
	root[rgb] = yx;
}

void calcarea(int x, int y)
{
#ifdef FIVExFIVE
	calcmean(x-2,y-2);
	calcmean(x-2,y-1);
	calcmean(x-2,y  );
	calcmean(x-2,y+1);
	calcmean(x-2,y+2);

	calcmean(x-1,y-2);
	calcmean(x-1,y-1);
	calcmean(x-1,y  );
	calcmean(x-1,y+1);
	calcmean(x-1,y+2);

	calcmean(x  ,y-2);
	calcmean(x  ,y-1);
	calcmean(x  ,y  );
	calcmean(x  ,y+1);
	calcmean(x  ,y+2);

	calcmean(x+1,y-2);
	calcmean(x+1,y-1);
	calcmean(x+1,y  );
	calcmean(x+1,y+1);
	calcmean(x+1,y+2);

	calcmean(x+2,y-2);
	calcmean(x+2,y-1);
	calcmean(x+2,y  );
	calcmean(x+2,y+1);
	calcmean(x+2,y+2);
#else
	calcmean(x-1,y-1);
	calcmean(x-1,y  );
	calcmean(x-1,y+1);

	calcmean(x  ,y-1);
	calcmean(x  ,y  );
	calcmean(x  ,y+1);

	calcmean(x+1,y-1);
	calcmean(x+1,y  );
	calcmean(x+1,y+1);
#endif
}

void rotate(int worse_x, int worse_y, int better_x, int better_y)
{

#if 0
	int x, y, dx, dy, lastyx, yx;
	x = better_x; y = better_y;
	lastyx = y*MAXW+x;
	P t = grid[lastyx];
        do {
		dx = worse_x-x; if (dx<0) dx = -1; if (dx>1) dx = 1;
		dy = worse_y-y; if (dy<0) dy = -1; if (dy>1) dy = 1;
		x += dx; y += dy;

		yx = y*MAXW+x;
		if (xyok[yx]) {
			grid[lastyx] = grid[yx];
			lastyx = yx;
		}
	} while (dx || dy);
	grid[lastyx] = t;

	gVersion++;

	x = better_x; y = better_y;
        do {
		dx = worse_x-x; if (dx<0) dx = -1; if (dx>1) dx = 1;
		dy = worse_y-y; if (dy<0) dy = -1; if (dy>1) dy = 1;
		calcarea(x, y);
		x += dx; y += dy;
	} while (dx || dy);
	calcarea(x, y);
#elif 1
	int x = better_x;
	int y = better_y;
	int dx =  abs(worse_x-x);
	int sx = x<worse_x ? 1 : -1;
	int dy = -abs(worse_y-y);
	int sy = y<worse_y ? 1 : -1;
	int err = dx+dy, e2; /* error value e_xy */

	int yx, lastyx = y*MAXW+x;
	P t = grid[lastyx];

	for(;;){  /* loop */
		if (x==worse_x && y==worse_y) break;
		e2 = 2*err;
		if (e2 >= dy) { err += dy; x += sx; } /* e_xy+e_x > 0 */
		if (e2 <= dx) { err += dx; y += sy; } /* e_xy+e_y < 0 */

		yx = y*MAXW+x;
		if (xyok[yx]) {
			grid[lastyx] = grid[yx];
			lastyx = yx;
		}
	}
	grid[lastyx] = t;

	gVersion++;

	x = better_x;
	y = better_y;
	dx =  abs(worse_x-x);
	sx = x<worse_x ? 1 : -1;
	dy = -abs(worse_y-y);
	sy = y<worse_y ? 1 : -1;
	err = dx+dy; /* error value e_xy */

	calcarea(x, y);

	for(;;){  /* loop */
		if (x==worse_x && y==worse_y) break;
		e2 = 2*err;
		if (e2 >= dy) { err += dy; x += sx; } /* e_xy+e_x > 0 */
		if (e2 <= dx) { err += dx; y += sy; } /* e_xy+e_y < 0 */

		calcarea(x, y);
	}
#else

#endif
}

int comparCenterDist(const void *_a, const void *_b)
{
	const int *a = (const int *) _a;
	const int *b = (const int *) _b;

	int x = *a%MAXW;
	int y = *a/MAXW;

	int da = (x-MAXW/2)*(x-MAXW/2) + (y-MAXH/2)*(y-MAXH/2);

	x = *b%MAXW;
	y = *b/MAXW;
	int db = (x-MAXW/2)*(x-MAXW/2) + (y-MAXH/2)*(y-MAXH/2);

	if (da < db)
		return -1; //!
	else if (da > db)
		return +1; //!
	else
		return 0;
}

int comparCenterDistArg(const void *_a, const void *_b, void *_arg)
{
	const int *a = (const int *) _a;
	const int *b = (const int *) _b;
	int *arg = (int *) _arg;
	int cx = arg[0];
	int cy = arg[1];

	int x = *a%MAXW;
	int y = *a/MAXW;

	int da = (x-cx)*(x-cx) + (y-cy)*(y-cy);

	x = *b%MAXW;
	y = *b/MAXW;
	int db = (x-cx)*(x-cx) + (y-cy)*(y-cy);

	if (da < db)
		return -1; //!
	else if (da > db)
		return +1; //!
	else
		return 0;
}

int comparPolar(const void *_a, const void *_b)
{
	const int *a = (const int *) _a;
	const int *b = (const int *) _b;

	int xa = (*a%MAXW) - (MAXW/2);
	int ya = (*a/MAXW) - (MAXH/2);

	double angleA = -atan2(ya,xa);
//	angleA -= ((xa)*(xa) + (ya)*(ya)) / 1e-10;


	int xb = (*b%MAXW) - (MAXW/2);
	int yb = (*b/MAXW) - (MAXH/2);

	double angleB = -atan2(yb,xb);
//	angleB -= ((xb)*(xb) + (yb)*(yb)) / 1e-10;

	angleA = fmod(angleA+M_PI, M_PI/2);
	angleB = fmod(angleB+M_PI, M_PI/2);

	if (angleA < angleB)
		return -1; //!
	else if (angleA > angleB)
		return +1; //!
	else
		return 0;
}

void usage(const char *argv0, int verbose)
{
	printf("Usage: %s [options] <output image> <input image> <numframes>\n", argv0);
	if (verbose)
		return;

	printf("\n\
options:\n\
	--seed=n	starting seed of random generator\n\
	--startframe=n	starting frame number\
	-v		\n\
	--verbose	show progress\n");
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
		enum {	LO_HELP=1, LO_VERBOSE='v', LO_SEED='s', LO_STARTFRAME};
		static struct option long_options[] = {
			/* name, has_arg, flag, val */
			{"help",       0, 0, LO_HELP},
			{"verbose",    0, 0, LO_VERBOSE},
			{"seed",       1, 0, LO_SEED},
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
		case LO_SEED:
			srand(strtol(optarg, NULL, 10));
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
	arg_output = argv[optind++];
	arg_json = argv[optind++];
	if (optind < argc)
		arg_image = argv[optind++];

	if (opt_timer) {
		signal(SIGALRM, sigAlarm);
		alarm(opt_timer);
	}

	/*
	** Load json
	*/
	json_error_t jerror;

	// load json
	FILE *f = fopen(arg_json, "r");
	if (!f) {
		fprintf(stderr, "fopen(%s) returned: %m\n", arg_json);
		exit(1);
	}
	json_t *jInput = json_loadf(f, 0, &jerror);
	if (jInput == 0) {
		fprintf(stderr, "Failed to load/decode %s line %d: %s\n", arg_json, jerror.line, jerror.text);
		exit(1);
	}
	fclose(f);

	// extract metrics
	opt_lead_frames = json_integer_value(json_object_get(jInput, "lead_frames"));
	opt_lead_fps_start = json_real_value(json_object_get(jInput, "lead_fps_start"));
	opt_lead_fps_end = json_real_value(json_object_get(jInput, "lead_fps_end"));

	/*
	 * Determine scan order
	 */
	{
		/*
		 * Circular grid. Horizontal center aligned.
		 *      1
		 *     2 3
		 *    4   5
		 *     6 7
		 *      8
		 */
		if (1) {
			// use rmap[] as intermediate
			for (int i = 0; i < MAXW * MAXH; i++)
				xyok[i] = rmap[i] = 0;
			numEval = 0;
			for (int d = 0; d < MAXH / 2; d++) {

				if (tick) {
					fprintf(stderr, "\r\e[K[%s] %.2f%% evalseq", timeAsString(), (100.0 * d) / (MAXH / 2));
					tick = 0;
				}
				for (int i = 0; i < MAXW * MAXH; i++) {
					int x = i % MAXW;
					int y = i / MAXW;

					if ((x - MAXW / 2) * (x - MAXW / 2) + (y - MAXH / 2) * (y - MAXH / 2) <= d * d && xyok[i] == 0) {
						xyok[i]++;
						evalseq[numEval++] = i;
					}
				}
			}

			/*
			 * Order above to polar (radar)
			 */
			if (0) {
				qsort(evalseq, numEval, sizeof(evalseq[0]), comparPolar);
			}

			/*
			 * Order above to random
			 */
			if (0) {
				for (uint32_t ix = 0; ix < numEval; ix++) {
					// get random index
					uint32_t other = (uint32_t) rand() % numEval;
					// swap
					int t = evalseq[ix];
					evalseq[ix] = evalseq[other];
					evalseq[other] = t;
				}
			}

			fprintf(stderr, "\r\e[K[%s]\n", timeAsString());
		}

		/*
		 * Rectangle, ordered starting from center moving towards edges.
		 */
		if (0) {
			for (int i = 0; i < MAXW * MAXH; i++) {
				xyok[i]++;
				evalseq[numEval++] = i;
			}
			qsort(evalseq, MAXW * MAXH, sizeof(evalseq[0]), comparCenterDist);
		}

		/*
		 * Rectangle, random
		 */
		if (0) {
			for (int i = 0; i < MAXW * MAXH; i++) {
				xyok[i]++;
				evalseq[numEval++] = i;
			}
			for (uint32_t ix = 0; ix < numEval; ix++) {
				// get random index
				uint32_t other = (uint32_t) rand() % numEval;
				// swap
				int t = evalseq[ix];
				evalseq[ix] = evalseq[other];
				evalseq[other] = t;
			}
		}

		/*
		 * Letter outlines
		 */
		if (0) {
			/* open source */
			fil = fopen("title.png", "rb");

			if (fil == NULL) {
				fprintf(stderr, "Could not open title.png\n");
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

			int xlo, xhi, ylo, yhi;

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

#if 1
					// for title letters
					// also assume fps_start=0.4
					// don't frget, first frame is hardcode loaded from title.png
//					if (r != 255 || g != 0 || b != 0)
//					if (r != 0 || g != 255 || b != 0)
//					if (r != 0 || g != 0 || b != 255)
//					if (r != 255 || g != 0 || b != 255)
					if (r != 0 || g != 255 || b != 255)
						continue;
#endif

					if (numEval == 0) {
						xlo=xhi = x;
						ylo=yhi = y;
					} else {
						if (x<xlo) xlo=x; else if (x>xhi) xhi=x;
						if (y<ylo) ylo=y; else if (y>yhi) yhi=y;
					}

					int yx = y*MAXW+x;
					xyok[yx]++;
					evalseq[numEval++] = yx;
				}
			}

			int arg[2];
			arg[0] = (xlo+xhi)/2; 
			arg[1] = (ylo+yhi)/2; 
			qsort_r(evalseq, numEval, sizeof(evalseq[0]), comparCenterDistArg, arg);
		}
		fprintf(stderr,"numEval:%d\n", numEval);
	}

	/*
	 * Open RGB colour mapping
	 */
	{
		asprintf(&opt_mmapName, "rgb-%d.%d.rgb", numEval, NUMSLOT);
		mmapfd = open(opt_mmapName, O_RDONLY, 0);
		if (mmapfd >= 0)
			nearestRGB = (uint32_t(*)[256*256*256][NUMSLOT]) mmap(NULL, sizeof(uint32_t) * 256 * 256 * 256 * NUMSLOT, PROT_READ, MAP_PRIVATE, mmapfd, 0);

		if (nearestRGB == MAP_FAILED) {
			nearestRGB = (uint32_t(*)[256*256*256][NUMSLOT]) calloc((size_t)256 * 256 * 256 * NUMSLOT, sizeof(uint32_t));
			nearestDist = (uint32_t(*)[256*256*256][NUMSLOT]) calloc((size_t)256 * 256 * 256 * NUMSLOT, sizeof(uint32_t));

			for (int i = 0; i < 256 * 256 * 256; i++) {
				for (int j=0; j<NUMSLOT; j++) {
					(*nearestRGB) [i][j] = 0xffffffff;
					(*nearestDist)[i][j] = 0xffffffff;
				}
			}
		}
	}

	for (int i=0; i<256*256*256; i++) {
		rmap[i] = -1;
		root[i] = -1;
	}
	for (int i=0; i<MAXW*MAXH; i++)
		mean[i].next = -1;

	if (arg_image) {
		/* open source */
		fil = fopen(arg_image, "rb");
		if (fil == NULL) {
			fprintf(stderr, "Could not open source image\n");
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

		for (uint32_t ix=0; ix<numEval; ix++) {
			int yx = evalseq[ix];
			int x = yx % MAXW;
			int y = yx / MAXW;
			assert(x < MAXW);
			assert(y < MAXH);

			uint32_t rgb = (uint32_t)gdImageGetTrueColorPixel(im, x, y);

			// save pixel

			if (mmapfd < 0) {
				// nearest RGB
				(*nearestRGB)[rgb][0] = rgb;
				(*nearestDist)[rgb][0] = 0;
			}

			if (rmap[rgb] != -1)
				fprintf(stderr,"dupcol: %06x\n", rgb);

			grid[yx].rgb = rgb;
			rmap[rgb] = yx;
		}
	} else {

		/*
		 * Generate pixels
		 * 1920*1080=2073600
		 * (2^7)x(2^7)x(2^7) = 2097152
		 */
		uint32_t err = 0;
		uint32_t r=0,g=0,b=0;
		for (uint32_t ix=0; ix<numEval; ix++) {
			int yx = evalseq[ix];

			// save pixel
			assert(r<256);
			assert(g<256);
			assert(b<256);

			uint32_t rgb = r << 16 | g << 8 | b;

			if (mmapfd < 0) {
				// nearest RGB
				(*nearestRGB)[rgb][0] = rgb;
				(*nearestDist)[rgb][0] = 0;
			}

			grid[yx].rgb = rgb;
			rmap[rgb] = yx;

			// increment to next value
			err += 256*256*256;
			while (err >= numEval) {
				b += 1;
				if (b >= 256) {
					g++;
					b -= 256;
					if (g >= 256) {
						r++;
						g -= 256;
					}
				}
				err -= numEval;
			}
		}

		fprintf(stderr, "nextRGB: r:%d g:%d b:%d err:%d\n", r,g,b, err);

		for (int i=0; i<numEval; i++) {
			// get random index
			uint32_t other = (uint32_t) rand() % numEval;
			// swap
			P t = grid[evalseq[i]];
			grid[evalseq[i]] = grid[evalseq[other]];
			grid[evalseq[other]] = t;
		}
	}

	if (mmapfd < 0) {
		// create if not loaded
		for (int margin = 1; ; margin++) {
			int numUpdate = 0;
			int numSkipped = 0;
			for (int r = 0; r < 256; r++) {
				for (int g = 0; g < 256; g++) {
					for (int b = 0; b < 256; b++) {
						uint32_t rgb = r << 16 | g << 8 | b;

						if (tick) {
							fprintf(stderr, "\r\e[K[%s] %.2f%% %d %d", timeAsString(), (100.0 * (r * 256 * 256 + g * 256 + b)) / (256 * 256 * 256), numUpdate, numSkipped);
							tick = 0;
						}

						if ((*nearestRGB)[rgb][0] != rgb && (*nearestDist)[rgb][0] != 0xffffffff)
							continue;
						if (root[rgb] >= 3 && (*nearestDist)[rgb][NUMSLOT-1] != 0xffffffff) {
							numSkipped++;
							continue;
						}

						int changed = 0;
						for (int tr = r - margin; tr <= r + margin; tr += 1) {
							for (int tg = g - margin; tg <= g + margin; tg += 1) {
								for (int tb = b - margin; tb <= b + margin; tb += 1) {
									if (tr == r - margin || tr == r + margin || tg == g - margin || tg == g + margin || tb == b - margin || tb == b + margin)
										if (tr >= 0 && tr < 256 && tg >= 0 && tg < 256 && tb >= 0 && tb < 256) {
											uint32_t trgb = tr << 16 | tg << 8 | tb;
											if ((*nearestDist)[trgb][0] != 0)
												continue; // must be a valid rgb combo

											// Taken from www.compuphase.com/cmetric.htm
											// dist = 3*dr*dr + 4*dg*dg + 2*db*db [for max red]
											// dist = 2*dr*dr + 4*dg*dg + 3*db*db [for zero red]
											// thus:
											// dist = (2+rmean/255)*dr*dr + 4*dg*dg + (3-rmean/255)*db*db
											//
											// dist*255 = 255*(2+rmean/255)*dr*dr + 255*4*dg*dg + 255*(3-rmean/255)*db*db
											// dist*255 = (255*2+rmean)*dr*dr + 255*4*dg*dg + (255*3-rmean)*db*db

											// get colour distance test pixel
											int rmean = (r + tr) / 2;
											int dr = r - tr;
											int dg = g - tg;
											int db = b - tb;
											int dist = (255 * 2 + rmean) * dr * dr + 255 * 4 * dg * dg + (255 * 3 - rmean) * db * db;
											assert(dist >= 0);

											int s = NUMSLOT - 2;
											while (s > 0 && dist < (*nearestDist)[rgb][s]) {
												(*nearestRGB)[rgb][s+1] = (*nearestRGB)[rgb][s];
												(*nearestDist)[rgb][s+1] = (*nearestDist)[rgb][s];
												s--;
												changed++;
											}
											if (dist < (*nearestDist)[rgb][s]) {
												assert(s == 0);
												(*nearestRGB)[rgb][s] = trgb;
												(*nearestDist)[rgb][s] = dist;
												changed++;
											} else if (s + 1 < NUMSLOT && dist < (*nearestDist)[rgb][s + 1]) {
												(*nearestRGB)[rgb][s + 1] = trgb;
												(*nearestDist)[rgb][s + 1] = dist;
												changed++;
											}
										}
								}
							}
						}

						if (changed) {
							numUpdate++;
							root[rgb] = 0;
						} else {
							root[rgb]++;
						}
					}
				}
			}

			fprintf(stderr, "\r\e[K[%s] margin:%d numUpdate:%d numSkipped:%d\n", timeAsString(), margin, numUpdate, numSkipped);
			if (numUpdate == 0)
				break;
		}
		for (uint32_t i = 0; i < 256 * 256 * 256; i++) {
			if ((*nearestRGB)[i][0] == i) {
				for (int s = 0; s < NUMSLOT; s++) {
					assert((*nearestDist)[i][s] != 0xffffffff);
				}
			}
		}
	}

	if (mmapfd < 0) {
		// create if not loaded
		int fd = open(opt_mmapName, O_WRONLY|O_CREAT|O_TRUNC, 0444);
		if (fd >= 0) {
			size_t left = sizeof(uint32_t) * 256 * 256 * 256 * NUMSLOT;
			char *p = (char*)nearestRGB;
			while (left >= 0x01000000) {
				ssize_t r = write(fd, p, 0x01000000);
				assert(r == 0x01000000);
				p += 0x01000000;
				left -= 0x01000000;
			}
			if (left > 0) {
				ssize_t r = write(fd, p, left);
				assert(r == left);
			}

			close(fd);
			fprintf(stderr, "\r\e[K[%s] Written %s\n", timeAsString(), opt_mmapName);
		}
	}

	tick = 1;
	for (int i=0; i<256*256*256; i++) {
		assert((*nearestRGB)[i][0] >= 0 && (*nearestRGB)[i][0] <= 0xffffff);
		if ((*nearestRGB)[i][0] == i) {
			if (tick) {
				fprintf(stderr, "\r\e[K[%s] %.2f%% tstnearest", timeAsString(), (100.0 * i) / (256*256*256));
				tick = 0;
			}

			for (int k = 0; k < NUMSLOT; k++) {
				assert((*nearestRGB)[i][k] >= 0 && (*nearestRGB)[i][k] <= 0xffffff);
			}
		}
	}
	fprintf(stderr, "\r\e[K[%s]\n", timeAsString());

	// calculate means
	for (int i=0; i<256*256*256; i++)
		root[i] = -1;

	gVersion++;
	//		for (int w=0; w<numEval; w++)
	for (int w=numEval-1; w>=0; w--)
	{
		int yx = evalseq[w];
		int x = yx % MAXW;
		int y = yx / MAXW;

		calcmean(x,y);
	}

	/*
	** Start
	*/

	double log_fps_start = log(opt_lead_fps_start);
	double log_fps_end = log(opt_lead_fps_end);
	uint64_t totalEval = 0;
	uint64_t nextSave = numEval;

	// video is running 16fps but simulates 0.1 fps.
	// time increments of 1/16 sec
	double fps;
	{
		double pos = (double)totalEval / (numEval * opt_lead_frames);
		double e = (log_fps_end-log_fps_start) * pos + log_fps_start;
		fps = exp(e);
		nextSave = numEval * fps / 16.0;
	}
	if (opt_framenr != 0)
		fps = 16.0;

	// save initial (not for a continuation)
	if (opt_framenr == 0)
		save(opt_framenr++);

	uint64_t lo = 0xffffffffffffffff;
	int locnt = 0;
	int oddeven = 0; // 0=O2I, 1=I2O
	for (;;) {

		int cntYes=0, cntNo = 0, cntRotate=0;
		gCntCalcArea=0;
		static uint64_t stats[NUMSLOT], statsOther;

		for (int i=0; i<NUMSLOT; i++)
			stats[i] = 0;
		statsOther = 0;

		for (int w=0; w<numEval; w++) {
			int worse = oddeven ? evalseq[w] : evalseq[numEval-1-w];
			int wx = worse % MAXW;
			int wy = worse / MAXW;

			if (tick) {
				fprintf(stderr, "\r\e[K[%s] %.2f%% %d %ld %lu [%lu %lu]", timeAsString(),
					100.0*w/numEval, cntRotate, gCntCalcArea, statsOther, totalEval, nextSave);
				tick = 0;
			}

			int found = -1;
			uint32_t rgb = grid[worse].rgb;
			for (int s = 0; s < NUMSLOT; s++) {

				int rgbworse = (*nearestRGB)[rgb][s];
				assert(rgbworse != 0xffffffff);
				if (root[rgbworse] != -1) {
					cntYes++;

					int better = root[rgbworse];
					int bx = better % MAXW;
					int by = better / MAXW;
					int betterDist = (bx-wx)*(bx-wx) + (by-wy)*(by-wy);
					for (int b=root[rgbworse]; b!=-1; b=mean[b].next) {
						bx = b % MAXW;
						by = b / MAXW;
						int d = (bx-wx)*(bx-wx) + (by-wy)*(by-wy);
						if (d <= betterDist) {
							// in case of draw, get the oldest
							better = b;
							betterDist = d;
						}
					}

					assert(better >= 0 && better < MAXW*MAXH);

					if (better != worse) {
						rotate(better % MAXW, better / MAXW, worse % MAXW, worse / MAXW);    // swap/shift
						cntRotate++;
					}
					stats[s]++;
					found = s;

					break;
				}
			}
			if (found==-1)
				statsOther++;

			if (fps < 16 && totalEval == nextSave) {
				save(opt_framenr++);

				double pos = (double)totalEval / (numEval * opt_lead_frames);
				double e = (log_fps_end-log_fps_start) * pos + log_fps_start;
				fps = exp(e);
				nextSave += numEval * fps / 16.0;

				fprintf(stderr," SAVE %f\n", fps);
			}

			totalEval++;
		}

		int cnt = 0;
		for (int ii=0; ii<numEval; ii++) {
			int i = evalseq[ii];
			for (int j = root[mean[i].rgb]; j != -1; j = mean[j].next) {
				if (i==j)
					cnt++;
			}
			assert(cnt>0);
		}
		fprintf(stderr,"\r\e[K[%s] %d [%d %d %d %d %lu] ", timeAsString(), opt_framenr, cntYes, cntNo, cnt, cntRotate, gCntCalcArea);
		for (int i=0; i<NUMSLOT; i++)
			fprintf(stderr,"%ld ", stats[i]);
		fprintf(stderr," [%ld] ", statsOther);

		if (gCntCalcArea < lo) {
			lo = gCntCalcArea;
			locnt = 0;
			fprintf(stderr," %d:%ld ### ", locnt, lo);
		} else {
			locnt++;
			fprintf(stderr," %d:%ld+%ld", locnt, lo, gCntCalcArea-lo);
		}
		fprintf(stderr,"\n");

		if (fps >= 16)
			save(opt_framenr++);

		oddeven = 1-oddeven;
	}

	return 0;
}
