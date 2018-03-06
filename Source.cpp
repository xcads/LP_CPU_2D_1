//generate the Lines' coefficient 
//# of lines : 999999
//# main line :3x-y>=-2 -x+y<=1 3x+4y<=12
//min: (-0.5,0.5) max:(2.14285714 15/7 , 3.14285714 22/7)
//domain other point(4/15,14/5)

//103553 lines
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits.h>
#include <math.h>
//DBL_EPSILON    2.2204460492503131e-016

struct Line {
	double a, b, c;
};


struct Point {
	double x, y;
};

typedef Line  Line;
typedef Point  Point;

//generate the Lines' coefficient around the feasible domain
//draw a circle and generate the line around the circle
double generateLineCoefficient(Line l1,Line l2, Line l3) {
	
	return true;
}


//test over
bool getCrossPoint(Line &l1, Line &l2,Point &p1)
{
	 p1.x = (l1.b*l2.c - l2.b*l1.c) /( l1.a*l2.b - l2.a*l1.b);
	 p1.y = (l2.c*l1.a - l1.c*l2.a) / (l1.a*l2.b - l2.a*l1.b);
	return true;
	/*double x = (l1.b*l2.c - l2.b*l1.c) / (l1.a*l2.b - l2.a*l1.b);
	//double y = (l1.c*l2.a - l2.c*l1.a) / (l1.a*l2.b - l2.a*l1.b);
	//printf("%lf   %lf", x, y);*/
}

//test over
double circleRadius(Point &p1, Point &p2, Point &p3,Point &center) {
	double r;
	double middleValueX1, middleValueY1, middleValueX2, middleValueY2, Slope1,  C1, Slope2,  C2;
	middleValueX1 = (p1.x + p3.x) / 2;
	middleValueY1 = (p1.y + p3.y) / 2;
	middleValueX2 = (p2.x + p3.x) / 2;
	middleValueY2 = (p2.y + p3.y) / 2;
	Slope1 = (p3.y - p1.y) / (p3.x - p1.x);
	C1 = middleValueY1 - Slope1*middleValueX1;
	Slope2 = (p3.y - p2.y) / (p3.x - p2.x);
	C2 = middleValueY2 - Slope2*middleValueX2;
	center.x = ( C2 -  C1) / (Slope1 - Slope2);
	center.y = Slope1 *center.x +  C1;
	r = sqrt((p1.x - center.x)*(p1.x - center.x) + (p1.y - center.y)*(p1.y - center.y));

	return r;
	//printf("%lf", drawCircle(p1, p2, p3));
}

//(ax0+bx0+c)/(sqrt(a2+b2))
double distance(Point p1, double a , double b ,double c) {
	double distance;
	distance = (p1.x*a + p1.y*b + c) / sqrt((a*a + b*b));

	return distance;
}

int main(void) {
	//# main line :3x-y>=-2 -x+y<=1 3x+4y<=12
	struct Point p1, p2, p3;
	struct Point center;

	struct Line l1 = { 3, -1, 2 };
	struct Line l2 = { -1, 1, -1 };
	struct Line l3 = { 3, 4, -12 };

	getCrossPoint(l1, l2, p1);
	getCrossPoint(l2, l3, p2);
	getCrossPoint(l1, l3, p3);
	
	double r = circleRadius(p1, p2, p3, center);
	//9900 <=
	//step 1 The first quadrant
	/*int a, b, c;
	for (a = 1; a < 10; a++) {
		for (b = 1; b < 10; b++) {
			for (c = -1; c > -100; c--) {
				if (fabs(distance(center, a, b, c)) >= r) {
					FILE *fp;
					fp = fopen("Coefficient.txt", "a");
					if (fp != NULL)
						fprintf(fp, "%d   %d   %d\n", a, b, c);
					printf("#: %d+%d+%d\n", a, b, c);
					fclose(fp);
				}
				
			}
		}
		
	}*/
	//step 2 The 2nd quadrant
	/*int a, b, c;
	for (a = -1; a > -100; a--) {
		for (b = 1; b < 20; b++) {
			for (c = 1; c < 20; c++) {
				if (fabs(distance(center, a, b, c)) >= r) {
					FILE *fp;
					fp = fopen("Coefficient.txt", "a");
					if (fp != NULL)
						fprintf(fp, "%d   %d   %d\n", a, b, c);
					printf("#: %d+%d+%d\n", a, b, c);
					fclose(fp);
				}

			}
		}
	}*/
	//step 3 The 3rd quadrant
	/*int a, b, c;
		for (a = -1; a > -200; a--) {
			for (b = -1; b > -200; b--) {
				for (c = -1; c > -200; c--) {
					if (fabs(distance(center, a, b, c)) >= r) {
						FILE *fp;
						fp = fopen("Coefficient.txt", "a");
						if (fp != NULL)
							fprintf(fp, "%d   %d   %d\n", a, b, c);
						printf("#: %d+%d+%d\n", a, b, c);
					fclose(fp);
				}

			}
		}
	}*/
		//step 4 The 4th quadrant
		int a, b, c;
		for (a = 1; a < 10; a++) {
			for (b = -1; b > -10; b--) {
				for (c = -1; c > -100; c--) {
					if (fabs(distance(center, a, b, c)) >= r) {
						FILE *fp;
						fp = fopen("Coefficient.txt", "a");
						if (fp != NULL)
							fprintf(fp, "%d   %d   %d\n", a, b, c);
						printf("#: %d+%d+%d\n", a, b, c);
						fclose(fp);
					}

				}
			}
		}
	return 0;
}
	

