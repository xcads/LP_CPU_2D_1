// cpu2d_VS.cpp : 定义控制台应用程序的入口点。
//
#define _CRT_SECURE_NO_WARNINGS
// #include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <float.h>

using namespace std;

#define EQUAL_NOT_VERTICAL  0
#define EQUAL_VERTICAL      1
#define GREATER_THAN        2
#define LESS_THAN           3

#define PI                  3.1415926535

#define TRUE 1
#define FALSE 0

typedef int BOOL;

// LINE structure : Constraints
struct Line {
	// a1x + a2y >= b
	double a1, a2, b;
	double lslope;
	bool beingUsed;
};

// Object Function
struct Objfunc {
	// xd = c1x + c2y
	double c1, c2;
};

// VERTEX structure
struct Vertex {
	double x, y;
};

// PAIR structure
struct Pair {
	int index;
	Line line1, line2;
	Vertex point;
	bool beingUsed;
};

typedef struct Line Line;
typedef struct Objfunc Objfunc;
typedef struct Vertex Vertex;

vector<struct Line> originalConstraints;




// Intersection vertex
bool Intersection(struct Line *l1, struct Line *l2, struct Vertex *v1)
{
	if (abs(l1->a1 * l2->a2 - l2->a1 * l1->a2) < 1e-6)
	{
		v1 = NULL;
		return false;
	}
	v1->x = -(l1->b * l2->a2 - l2->b * l1->a2) / (l1->a1 * l2->a2 - l2->a1 * l1->a2);
	v1->y = (l1->b * l2->a1 - l2->b * l1->a1) / (l1->a2 * l2->a1 - l1->a1 * l2->a2);
	printf("Intersection: (%lf, %lf)\n", v1->x, v1->y);
	return true;
}

// Slope line
bool Slope(struct Line *l)
{
	if (fabs(l->a2 - 0.0) < 1e-6)
	{
		if ((l->a1 > 0 && l->a2 < 0) || (l->a1 < 0 && l->a2 > 0))
		{
			l->lslope = FLT_MAX;
		}
		else if ((l->a1 < 0 && l->a2 < 0) || (l->a1 > 0 && l->a2 > 0))
		{
			l->lslope = -FLT_MAX;
		}
		else
		{
			l->lslope = -l->a1 / l->a2;
		}
		return false;
	}
	l->lslope = -l->a1 / l->a2;
	return true;
}

// Compare
int cmp(const void *a, const void *b)
{
	struct Line *aa = (struct Line *)a;
	struct Line *bb = (struct Line *)b;
	return ((aa->lslope > bb->lslope) ? 1 : -1);
}

// Rotation - O(n)
bool Rotation(struct Line lines[], struct Objfunc object, int index, int *numG, int *numH)
{
	double thetaArc = atan(-object.c1 / object.c2);
	double thetaDec = atan(-object.c1 / object.c2) * 180 / PI;

	int i;
	double a1Temp, a2Temp, bTemp;

	(*numG) = 0;
	(*numH) = 0;

	for (i = 0; i < index; i += 1) {
		a1Temp = originalConstraints[i].a1;
		a2Temp = originalConstraints[i].a2;
		bTemp = originalConstraints[i].b;

		lines[i].a1 = cos(thetaArc) * a1Temp + sin(thetaArc) * a2Temp;
		lines[i].a2 = cos(thetaArc) * a2Temp - sin(thetaArc) * a1Temp;
		lines[i].b = bTemp;

		if (lines[i].a2 > 0) {
			(*numG)++;
			printf("%d", (*numG));
		}
		else if (lines[i].a2 < 0) {
			(*numH)++;
			printf("%d", (*numH));
		}
		else {
			return false;
		}

		Slope(&lines[i]);
		lines[i].beingUsed = true;

	}

	if ((*numG) + (*numH) != index) {
		printf("Fatal Error at Rotation()!\n");
		exit(-1);
	}

	return true;
}

// Separation - O(n)
bool Segmentation(struct Line I1[], struct Line I2[], struct Line lines[], int numG, int numH)
{
	int index = numG + numH;
	int i, g = 0, h = 0;
	for (i = 0; i < index; i++) {
		if (lines[i].a2 > 0) {
			I1[g].a1 = -lines[i].a1 / lines[i].a2;
			I1[g].a2 = 1;
			I1[g].b = lines[i].b / lines[i].a2;
			Slope(&I1[g]);
			I1[g].lslope = -I1[g].lslope;
			I1[g].beingUsed = true;
			g++;
		}
		else if (lines[i].a2 < 0) {
			I2[h].a1 = -lines[i].a1 / lines[i].a2;
			I2[h].a2 = 1;
			I2[h].b = lines[i].b / lines[i].a2;
			Slope(&I2[h]);
			I2[h].lslope = -I2[h].lslope;
			I2[h].beingUsed = true;
			h++;
		}
		else {
			return false;
		}
	}
	return true;
}

// New Make Pairs
bool NewMakePairs(struct Line I1[], struct Line I2[],
	struct Pair pair,
	int numG, int numH, int *index,
	double leftBound, double rightBound)
{
	return true;
}


// Make pairs
bool MakePairs(struct Line I1[], struct Line I2[],
	struct Pair pairsG[], struct Pair pairsH[],
	int numG, int numH, int *index,
	double leftBound, double rightBound)
{
	int g, h, gtemp;
	*index = 0;
	for (g = 0; g < numG; g += 1) {
		// drop
		if (I1[g].beingUsed == false) {
			continue;
		}
		for (gtemp = g + 1; gtemp < numG; gtemp++) {
			if (I1[gtemp].beingUsed == true) {
				break;
			}
		}
		if (gtemp == numG) break;

		if (abs(I1[g].lslope - I1[gtemp].lslope) < 1e-6) {
			if (I1[g].b > I1[gtemp].b) {
				I1[gtemp].beingUsed = false;
			}
			else {
				I1[g].beingUsed = false;
			}
			g = gtemp - 1;
			continue;
		}
		struct Vertex *p = (struct Vertex *)malloc(sizeof(struct Vertex));
		Intersection(&I1[g], &I1[gtemp], p);
		if (p->x < leftBound || p->x > rightBound) {
			if (abs(I1[g].lslope) > abs(I1[gtemp].lslope)) {
				I1[gtemp].beingUsed = false;
			}
			else if (abs(I1[gtemp].lslope) < abs(I1[gtemp].lslope)) {
				I1[g].beingUsed = false;
			}
			g = gtemp - 1;
			continue;
		}
		pairsG[(*index)].index = (*index);
		pairsG[(*index)].line1 = I1[g];
		pairsG[(*index)].line2 = I1[gtemp];
		pairsG[(*index)].point.x = p->x; pairsG[(*index)].point.y = p->y;
		printf("Intersection2: (%lf, %lf)\n", p->x, p->y);
		printf("Value: %lf, %lf\n", I1[g].a1 * p->x + I1[g].b, I1[gtemp].a1 * p->x + I1[gtemp].b);
		(*index)++;
	}
	/*
	for (h = 0; h < numH; h += 2) {
	// drop

	pairsH[h / 2].index = *index++;
	pairsH[h / 2].line1 = I2[h];
	pairsH[h / 2].line2 = I2[h + 1];
	Intersection(&I2[h], &I2[h + 1], &pairsH[h / 2].point);
	}*/

	return true;
}

// sg, Sg, sh, Sh
struct Vertex *TestingLine(struct Pair pairsG[], struct Pair pairsH[],
	struct Line I1[], struct Line I2[],
	int numG, int numH, int numDot)
{
	// Randomly choose a point
	srand((unsigned int)time(NULL));
	int index = rand() % numDot;
	double xPrimeG = pairsG[index].point.x;   // x' - xPrime
	double yPrimeG = pairsG[index].point.y;
	double yPrimeH;

	struct Line *sg = (&pairsG[index].line1.a1 < &pairsG[index].line2.a1) ? &pairsG[index].line1 : &pairsG[index].line2;
	struct Line *Sg = (&pairsG[index].line1.a1 >= &pairsG[index].line2.a1) ? &pairsG[index].line1 : &pairsG[index].line2;
	struct Line *sh = NULL;
	struct Line *Sh = NULL;
	// struct Line *sh = (&pairsH[index].line1.a1 < &pairsH[index].line2.a1) ? &pairsH[index].line1 : &pairsH[index].line2;
	// struct Line *Sh = (&pairsH[index].line1.a1 < &pairsH[index].line2.a1) ? &pairsH[index].line1 : &pairsH[index].line2;

	vector<int> linesG;
	vector<int> linesH;

	// Finding g(x') and H(x')
	for (int i = 0; i < numG; i++) {
		if (I1[i].beingUsed == true) {
			if (abs(yPrimeG - (I1[i].a1 * xPrimeG + I1[i].b)) > 1e-6 && yPrimeG < (I1[i].a1 * xPrimeG + I1[i].b)) {
				printf("xPrime yPrime ???: %lf %lf %lf\n", xPrimeG, yPrimeG, (I1[i].a1 * xPrimeG + I1[i].b));
				yPrimeG = I1[i].a1 * xPrimeG + I1[i].b;
				sg = &I1[i];
				Sg = &I1[i];
			}
		}
	}
	for (int i = 0; i < numH; i++) {
		if (I2[i].beingUsed == true) {
			if (sh == NULL || Sh == NULL) {
				sh = &I2[i];
				Sh = &I2[i];
				yPrimeH = I2[i].a1 * xPrimeG + I2[i].b;
			}
			else if (abs(yPrimeH - (I2[i].a1 * xPrimeG + I2[i].b)) > 1e-6 && yPrimeH > (I2[i].a1 * xPrimeG + I2[i].b)) {
				yPrimeH = I2[i].a1 * xPrimeG + I2[i].b;
				sh = &I2[i];
				Sh = &I2[i];
			}
		}
	}
	if (numH == 0) {
		yPrimeH = yPrimeG + 1000.0;
	}

	// Finding sg - min g(x') && Finding Sg - max g(x')
	/*
	struct Line *sg = &pairsG[0].line1;
	struct Line *Sg = &pairsG[0].line1;
	struct Line *sh = &pairsH[0].line1;
	struct Line *Sh = &pairsH[0].line1;
	*/
	for (int i = 0; i < numG; i++) {
		double currentLineValueG = I1[i].a1 * xPrimeG + I1[i].b;
		if (I1[i].beingUsed == false || abs(currentLineValueG - yPrimeG) >= 1e-6) {
			continue;
		}

		if (I1[i].a1 < sg->a1) {
			sg = &I1[i];
		}
		if (I1[i].a1 > Sg->a1) {
			Sg = &I1[i];
		}
	}
	// Finding sh - min h(x') && Finding Sh - max h(x')
	for (int i = 0; i < numH; i++) {
		double currentLineValueH = I2[i].a1 * xPrimeG + I2[i].b;
		if (I2[i].beingUsed == false || abs(currentLineValueH - yPrimeH) >= 1e-6) {
			continue;
		}

		if (I2[i].a1 < sh->a1) {
			sh = &I2[i];
		}
		if (I2[i].a1 > Sh->a1) {
			Sh = &I2[i];
		}
	}

	// Is feasible
	if (abs(yPrimeG - yPrimeH) < 1e-6) {
		if (sg->a1 > 0 && sg->a1 >= Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->beingUsed = false;
			}
			Sg->beingUsed = false;
			return NULL;
		}
		else if (Sg->a1 < 0 && Sg->a1 <= sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->beingUsed = false;
			}
			sg->beingUsed = false;
			return NULL;
		}
		else {
			// x* = x'
			return &(pairsG[index].point);
		}
	}
	else if (yPrimeG > yPrimeH) {   // infeasible
		if (sg->a1 > Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->beingUsed = false;
			}
			Sg->beingUsed = false;
			return NULL;
		}
		else if (Sg->a1 < sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->beingUsed = false;
			}
			sg->beingUsed = false;
			return NULL;
		}
		else if ((sg->a1 - Sh->a1) <= 0 && 0 <= (Sg->a1 - sh->a1)) {
			// no feasible
			printf("No feasible solution!\n");
			exit(0);
			return NULL;
		}
	}
	else if (yPrimeG < yPrimeH) {   // feasible
		if (sg->a1 > 0) {
			// x* < x'
			Sg->beingUsed = false;
			return NULL;
		}
		else if (Sg->a1 < 0) {
			// x* > x'
			sg->beingUsed = false;
			return NULL;
		}
		else if (sg->a1 <= 0 && 0 <= Sg->a1) {
			// x* = x'
			return &(pairsG[index].point);
		}
	}
}

void LinearProgramming(void)
{
	int indexRecord = 0;
	int numGRecord;
	int numHRecord;
	int indexPair;
	double leftBound, rightBound;
	double aTemp, bTemp, cTemp;
	bool judge = false;
	struct Objfunc object;

	while (1) {
		scanf("%lf%lf%lf", &aTemp, &bTemp, &cTemp);
		if (aTemp == 0.0 && bTemp == 0.0 && cTemp == 0.0) {
			break;
		}
		struct Line lineTemp;
		lineTemp.a1 = aTemp;
		lineTemp.a2 = bTemp;
		lineTemp.b = cTemp;
		originalConstraints.push_back(lineTemp);
		indexRecord++;
	}
	scanf("%lf%lf", &object.c1, &object.c2);
	scanf("%lf%lf", &leftBound, &rightBound);

	struct Line *lines = (struct Line *)malloc(indexRecord * sizeof(struct Line));
	struct Line *I1 = (struct Line *)malloc(indexRecord * sizeof(struct Line));
	struct Line *I2 = (struct Line *)malloc(indexRecord * sizeof(struct Line));
	struct Pair *pairG = (struct Pair *)malloc(indexRecord * sizeof(struct Pair) / 2);
	struct Pair *pairH = (struct Pair *)malloc(indexRecord * sizeof(struct Pair) / 2);
	struct Vertex *sln = NULL;

	judge = Rotation(lines, object, indexRecord, &numGRecord, &numHRecord);
	if (judge == false) {
		printf("Fatal Error at LinearProgramming() - Rotation()!\n");
		exit(-1);
	}

	judge = Segmentation(I1, I2, lines, numGRecord, numHRecord);
	if (judge == false) {
		printf("Fatal Error at LinearProgramming() - Segmentation()!\n");
		exit(-1);
	}

	while (1) {
		judge = MakePairs(I1, I2, pairG, pairH, numGRecord, numHRecord, &indexPair, leftBound, rightBound);
		if (judge == false) {
			printf("Fatal Error at LinearProgramming() - MakePairs()!\n");
			exit(-1);
		}

		sln = TestingLine(pairG, pairH, I1, I2, numGRecord, numHRecord, indexPair);
		if (sln != NULL) {
			break;
		}
	}

	printf("sln: %lf %lf", sln->x, sln->y);

	return;

}


int main()
{

	LinearProgramming();

	return 0;
}

