#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include<limits.h>
#include<vector>
#include<time.h>
#include<algorithm>

typedef int BOOL;
#define TRUE  1
#define FALSE 0

/*
#define NO_VERTICAL   0
#define VERTICAL      1
#define MORE_THAN     2
#define LESS_THAN     3
*/
#define PI 3.1415926

struct Line {
	float a1, a2, b;
	float lslope;	
};

struct Pair {
	float index,line1, line2, point;
	double x, y;
};

// Object Function
struct Obfunc {
	// xd = c1x + c2y + c
	float c1, c2, c;
};

// Establish VERTEX structure
struct Vertex {
	float x, y;
};
// some struct about the questions
typedef struct Line Line;
typedef struct Vertex Vertex;

// find some intersection
BOOL intersection(struct Line l1[], struct Line l2[], struct Vertex &v1)
{
	if ((l1->a1 * l2->a2 - l2->a1 * l1->a2) < 1e-10)
	{
		return FALSE;
	}
	v1.x = (l1->b * l2->a2 - l2->b * l1->a2) / (l1->a1 * l2->a2 - l2->a1 * l1->a2);
	v1.y = (l1->b * l2->a1 - l2->b * l1->a1) / (l1->a2 * l2->a1 - l1->a1 * l2->a2);
	return TRUE;
}

// Slope line
BOOL slope(struct Line &l)
{
	if (fabs(l.a2 - 0.0) < 1e-10)
	{
		if ((l.a1 > 0 && l.a2 < 0) || (l.a1 < 0 && l.a2 > 0))
		{
			l.lslope = FLT_MAX;
		}
		else if ((l.a1 < 0 && l.a2 < 0) || (l.a1 > 0 && l.a2 > 0))
		{
			l.lslope = -FLT_MAX;
		}
		else
		{
			l.lslope = -l.a1 / l.a2;
		}
		return FALSE;
	}
	l.lslope = -l.a1 / l.a2;
	return TRUE;
}

//trans lines to x-aixs
int transformation(struct Line lines[], struct Obfunc object, int index, int *num_g, int *num_h) {
	double the_arc = atan(-object.c1 / object.c2);
	double the_dec = atan(-object.c1 / object.c2) * 180 / PI;

	int i;
	double a1_t, a2_t, b_t;

	*num_g = 0;
	*num_h = 0;

	for (i = 0; i < index; i++) {
		a1_t = lines[i].a1;
		a2_t = lines[i].a2;
		b_t = lines[i].b;

		lines[i].a1 = cos(the_arc)*a1_t + sin(the_arc)*a2_t;
		lines[i].a2 = cos(the_arc)*a2_t - sin(the_arc)*a1_t;

		if (lines[i].a2 > 0) {
			*num_g++;
		}
		else if (lines[i].a2 < 0) {
			*num_h++;
		}
		else {
			return false;
		}
		slope(lines[i]);
	}
	if (*num_g + *num_h != index)
	{
		printf("error!\n");
		exit(-1);
	}
}

// compare funcuntion
int cmp(const void *a, const void *b)
{
	struct Line *aa = (struct Line *)a;
	struct Line *bb = (struct Line *)b;
	return ((aa->lslope > bb->lslope) ? 1 : -1);
}

double make_Pair(struct Line L1[], struct Line L2[], struct Pair p1[], struct Pair p2[], int num_g, int num_h, int *index) {
	int i, j;
	*index = 0;
	for (i = 0; i < num_g; i += 2) {
		if (L1[i].lslope - L1[i + 1].lslope < 1e-10) {
			L1[i].lslope = NULL;
		}
	}
	p1[i / 2].index = *index++;
	p1[i / 2].line1 = *index++;
	p1[i / 2].line2 = *index++;
	intersection(&L1[i], &L1[i + 1], &p1[i / 2].point);
return true;
}

bool median_find(struct Line L1[], struct Line L2[], struct Line lines[], int num_g, int num_h)
{
	int i = num_g + num_h;
	int h = 0, n = 0;
	int j = 0;
	for (; j < i; j++) {
		if (lines[j].a2>0) {
			L1[h].a1 = -lines[j].a1 / lines[j].a2;
			L1[h].a2 = 1;
			L1[h].b = lines[j].b / lines[j].a2;
			h++;
		}
		else if (lines[j].a2 < 0) {
			L1[n].a1 = -lines[j].a1 / lines[j].a2;
			L1[n].a2 = 1;
			L1[n].b = lines[j].b / lines[j].a2;
		}
		else {
			return false;
		}
	}
	return true;
}

bool test(struct Pair p1[], struct Pair p2[], struct Line l1[], struct Line l2[],int index,
			int num_g , int num_h , int num_d)
{
	/*srand((unsigned int)time(NULL));*/
	float x_g = p1[index].x;
	float y_g = p1[index].y;

	struct Line *sg = (l1[index].a1)?&


}
//main function
int main(void) {

	return 0;
}

