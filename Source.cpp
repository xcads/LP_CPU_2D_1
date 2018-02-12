/*
Author	: @Johanan Luo
Time	: 2018/2/11
Input	:
OutPut	:

*/
#define _CRT_SECURE_NO_WARNINGS
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

#define PI  3.1415926

#define TRUE 1
#define FALSE 0

typedef int BOOL;

struct Line {
	// a1x + a2y >= b
	double a1, a2, b;
	double slope;
	bool symbol;
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
	Line l1, l2;
	Vertex v1;
	bool symbol;
};

typedef struct Line Line;
typedef struct Objfunc Objfunc;
typedef struct Vertex Vertex;

vector<struct Line> originalC;




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
//		printf("Intersection: (%lf, %lf)\n", v1->x, v1->y);
	return true;
}

// Slope line
bool Slope(struct Line *l)
{
	if (fabs(l->a2 - 0.0) < 1e-6)
	{
		if ((l->a1 > 0 && l->a2 < 0) || (l->a1 < 0 && l->a2 > 0))
		{
			l->slope = FLT_MAX;
		}
		else if ((l->a1 < 0 && l->a2 < 0) || (l->a1 > 0 && l->a2 > 0))
		{
			l->slope = -FLT_MAX;
		}
		else
		{
			l->slope = -l->a1 / l->a2;
		}
		return false;
	}
	l->slope = -l->a1 / l->a2;
	return true;
}

// Compare
int compare(const void *a, const void *b)
{
	struct Line *aa = (struct Line *)a;
	struct Line *bb = (struct Line *)b;
	return ((aa->slope > bb->slope) ? 1 : -1);
}

// 
bool Transformation(struct Line lines[], struct Objfunc object, int index, int *numG, int *numH)
{
	double thetaArc = atan(-object.c1 / object.c2);
	double thetaDec = atan(-object.c1 / object.c2) * 180 / PI;

	int i;
	double a1Temp, a2Temp, bTemp;

	(*numG) = 0;
	(*numH) = 0;

	for (i = 0; i < index; i += 1) {
		a1Temp = originalC[i].a1;
		a2Temp = originalC[i].a2;
		bTemp = originalC[i].b;

		lines[i].a1 = cos(thetaArc) * a1Temp + sin(thetaArc) * a2Temp;
		lines[i].a2 = cos(thetaArc) * a2Temp - sin(thetaArc) * a1Temp;
		lines[i].b = bTemp;

		if (lines[i].a2 > 0) {
			(*numG)++;
	//		printf("%d", (*numG));
		}
		else if (lines[i].a2 < 0) {
			(*numH)++;
	//		printf("%d", (*numH));
		}
		else {
			return false;
		}

		Slope(&lines[i]);
		lines[i].symbol = true;

	}

	if ((*numG) + (*numH) != index) {
		//printf("Error in transformation ()!\n");
		exit(-1);
	}

	return true;
}

// 
bool Judge(struct Line l1[], struct Line l2[], struct Line lines[], int numG, int numH)
{
	int index = numG + numH;
	int i, g = 0, h = 0;
	for (i = 0; i < index; i++) {
		if (lines[i].a2 > 0) {
			l1[g].a1 = -lines[i].a1 / lines[i].a2;
			l1[g].a2 = 1;
			l1[g].b = lines[i].b / lines[i].a2;
			Slope(&l1[g]);
			l1[g].slope = -l1[g].slope;
			l1[g].symbol = true;
			g++;
		}
		else if (lines[i].a2 < 0) {
			l2[h].a1 = -lines[i].a1 / lines[i].a2;
			l2[h].a2 = 1;
			l2[h].b = lines[i].b / lines[i].a2;
			Slope(&l2[h]);
			l2[h].slope = -l2[h].slope;
			l2[h].symbol = true;
			h++;
		}
		else {
			return false;
		}
	}
	return true;
}

// Make pairs
bool Pairs(struct Line l1[], struct Line l2[],
	struct Pair pairsG[], struct Pair pairsH[],
	int numG, int numH, int *index,
	double leftBound, double rightBound)
{
	int g, h, gtemp;
	*index = 0;
	for (g = 0; g < numG; g += 1) {
		// drop
		if (l1[g].symbol == false) {
			continue;
		}
		for (gtemp = g + 1; gtemp < numG; gtemp++) {
			if (l1[gtemp].symbol == true) {
				break;
			}
		}
		if (gtemp == numG) break;

		if (abs(l1[g].slope - l1[gtemp].slope) < 1e-6) {
			if (l1[g].b > l1[gtemp].b) {
				l1[gtemp].symbol= false;
			}
			else {
				l1[g].symbol = false;
			}
			g = gtemp - 1;
			continue;
		}
		struct Vertex *p = (struct Vertex *)malloc(sizeof(struct Vertex));
		Intersection(&l1[g], &l1[gtemp], p);
		if (p->x < leftBound || p->x > rightBound) {
			if (abs(l1[g].slope) > abs(l1[gtemp].slope)) {
				l1[gtemp].symbol= false;
			}
			else if (abs(l1[gtemp].slope) < abs(l1[gtemp].slope)) {
				l1[g].symbol = false;
			}
			g = gtemp - 1;
			continue;
		}
		pairsG[(*index)].index = (*index);
		pairsG[(*index)].l1 = l1[g];
		pairsG[(*index)].l2 = l1[gtemp];
		pairsG[(*index)].v1.x = p->x; pairsG[(*index)].v1.y = p->y;
		//printf("Intersection2: (%lf, %lf)\n", p->x, p->y);
	//printf("Value: %lf, %lf\n", l1[g].a1 * p->x + l1[g].b, l1[gtemp].a1 * p->x + l1[gtemp].b);
		(*index)++;
	}

	return true;
}

// sg, Sg, sh, Sh
struct Vertex *TestingLine(struct Pair pairsG[], struct Pair pairsH[],
	struct Line l1[], struct Line l2[],
	int numG, int numH, int numDot)
{
	// Randomly choose a point
	srand((unsigned int)time(NULL));
	int index = rand() % numDot;
	double xPrimeG = pairsG[index].v1.x;   // x' - xPrime
	double yPrimeG = pairsG[index].v1.y;
	double yPrimeH;

	struct Line *sg = (&pairsG[index].l1.a1 < &pairsG[index].l2.a1) ? &pairsG[index].l1 : &pairsG[index].l2;
	struct Line *Sg = (&pairsG[index].l1.a1 >= &pairsG[index].l2.a1) ? &pairsG[index].l1 : &pairsG[index].l2;
	struct Line *sh = NULL;
	struct Line *Sh = NULL;
	
	vector<int> linesG;
	vector<int> linesH;

	for (int i = 0; i < numG; i++) {
		if (l1[i].symbol == true) {
			if (abs(yPrimeG - (l1[i].a1 * xPrimeG + l1[i].b)) > 1e-6 && yPrimeG < (l1[i].a1 * xPrimeG + l1[i].b)) {
				printf("x* y* : %lf %lf \n", xPrimeG, yPrimeG );
				yPrimeG = l1[i].a1 * xPrimeG + l1[i].b;
				sg = &l1[i];
				Sg = &l1[i];
			}
		}
	}
	for (int i = 0; i < numH; i++) {
		if (l2[i].symbol == true) {
			if (sh == NULL || Sh == NULL) {
				sh = &l2[i];
				Sh = &l2[i];
				yPrimeH = l2[i].a1 * xPrimeG + l2[i].b;
			}
			else if (abs(yPrimeH - (l2[i].a1 * xPrimeG + l2[i].b)) > 1e-6 && yPrimeH > (l2[i].a1 * xPrimeG + l2[i].b)) {
				yPrimeH = l2[i].a1 * xPrimeG + l2[i].b;
				sh = &l2[i];
				Sh = &l2[i];
			}
		}
	}
	if (numH == 0) {
		yPrimeH = yPrimeG + 1000.0;
	}

	for (int i = 0; i < numG; i++) {
		double currentLineValueG = l1[i].a1 * xPrimeG + l1[i].b;
		if (l1[i].symbol == false || abs(currentLineValueG - yPrimeG) >= 1e-6) {
			continue;
		}

		if (l1[i].a1 < sg->a1) {
			sg = &l1[i];
		}
		if (l1[i].a1 > Sg->a1) {
			Sg = &l1[i];
		}
	}


	for (int i = 0; i < numH; i++) {
		double currentLineValueH = l2[i].a1 * xPrimeG + l2[i].b;
		if (l2[i].symbol == false || abs(currentLineValueH - yPrimeH) >= 1e-6) {
			continue;
		}

		if (l2[i].a1 < sh->a1) {
			sh = &l2[i];
		}
		if (l2[i].a1 > Sh->a1) {
			Sh = &l2[i];
		}
	}


	if (abs(yPrimeG - yPrimeH) < 1e-6) {
		if (sg->a1 > 0 && sg->a1 >= Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->symbol = false;
			}
			Sg->symbol = false;
			return NULL;
		}
		else if (Sg->a1 < 0 && Sg->a1 <= sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->symbol = false;
			}
			sg->symbol = false;
			return NULL;
		}
		else {
			// x* = x'
			return &(pairsG[index].v1);
		}
	}
	else if (yPrimeG > yPrimeH) {  
		if (sg->a1 > Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->symbol = false;
			}
			Sg->symbol = false;
			return NULL;
		}
		else if (Sg->a1 < sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->symbol = false;
			}
			sg->symbol = false;
			return NULL;
		}
		else if ((sg->a1 - Sh->a1) <= 0 && 0 <= (Sg->a1 - sh->a1)) {
			// no feasible
			printf("No feasible solution!\n");
			exit(0);
			return NULL;
		}
	}
	else if (yPrimeG < yPrimeH) {   
		if (sg->a1 > 0) {
			// x* < x'
			Sg->symbol = false;
			return NULL;
		}
		else if (Sg->a1 < 0) {
			// x* > x'
			sg->symbol = false;
			return NULL;
		}
		else if (sg->a1 <= 0 && 0 <= Sg->a1) {
			// x* = x'
			return &(pairsG[index].v1);
		}
	}
}

void LP(void)
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
		originalC.push_back(lineTemp);
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

	judge = Transformation(lines, object, indexRecord, &numGRecord, &numHRecord);
	if (judge == false) {
	//	printf("Fatal Error at LinearProgramming() - Rotation()!\n");
		exit(-1);
	}

	judge =Judge(I1, I2, lines, numGRecord, numHRecord);
	if (judge == false) {
	//	printf("Fatal Error at LinearProgramming() - Segmentation()!\n");
		exit(-1);
	}

	while (1) {
		judge = Pairs(I1, I2, pairG, pairH, numGRecord, numHRecord, &indexPair, leftBound, rightBound);
		if (judge == false) {
	//		printf("Fatal Error at LinearProgramming() - MakePairs()!\n");
			exit(-1);
		}

		sln = TestingLine(pairG, pairH, I1, I2, numGRecord, numHRecord, indexPair);
		if (sln != NULL) {
			break;
		}
	}

	printf("last answer : %lf %lf", sln->x, sln->y);

	return;

}


int main()
{

	LP();

	return 0;
}

