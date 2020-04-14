#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
using namespace std;
const int SIZE = 8;
int main()
{
	double x[SIZE + 1] = { 0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0 };
	double y[SIZE + 1] = { 1.0,0.979915,0.927295,0.858001,0.785398,0.716844,0.716844,0.716844,0.716844 };
	double step = 0.25;
	//Spline function interpolation
	double a[SIZE] = {};
	double b[SIZE] = {};
	double c[SIZE + 1] = {};
	double d[SIZE] = {};
	double delta[SIZE - 1] = {};
	double lambda[SIZE - 1] = {};
	//sweep method
	for (int i = 0; i < (SIZE - 1); i++) {
		double l = ((y[i + 2] - y[i + 1]) - (y[i + 1] - y[i])) / step;
		if (i == 0) {
			delta[i] = -(step / (4 * step));//-h2/(2(h1+h2))  h1=x1-x0=step
			lambda[i] = 3 * l / (4 * step);
		}
		else {
			delta[i] = -(step / (4 * step + step * delta[i - 1]));
			lambda[i] = (3 * l - step * lambda[i - 1]) / (4 * step + step * delta[i - 1]);
		}
	}
	c[SIZE] = 0;
	for (int i = SIZE - 1; i > 0; i--) {
		c[i] = delta[i - 1] * c[i + 1] + lambda[i - 1];
	}
	for (int i = 0; i < SIZE; i++) {
		double l = (y[i + 1] - y[i]) / step;
		b[i] = l + (2 * c[i + 1] * step + step * c[i]) / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * step);
		a[i] = y[i + 1];
	}
	//Show points

	ofstream fout;
	//fout.open("MyFunc.txt");
	double h = 0.05;

	int k = 1;
	/*
	for (double i = x[0]; i <= x[SIZE]; k++) {
		for (; i < x[k] + h; i = i + h) {
			fout << i << "\t";
			fout << a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);//Spline
			fout << endl;
		}
	}
	*/
	//Trapezoid method
	///*

	k = 1;
	double Integ1 = 0.0, Integ2 = 0.0;
		for (double i = x[0]; i <= x[SIZE]; k++) {
			for (; i < x[k] + h; i = i + h) {
				double y1 = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
				double y2 = a[k - 1] + b[k - 1] * (i + h - x[k]) + c[k] * (i + h - x[k]) * (i + h - x[k]) + d[k - 1] * (i + h - x[k]) * (i + h - x[k]) * (i + h - x[k]);
				Integ1 += h * (y1 + y2) / 2;
			}
		}
		h = h / 2;
		k = 1;
		for (double i = x[0]; i <= x[SIZE]; k++) {
			for (; i < x[k] + h; i = i + h) {
				double y1 = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
				double y2 = a[k - 1] + b[k - 1] * (i + h - x[k]) + c[k] * (i + h - x[k]) * (i + h - x[k]) + d[k - 1] * (i + h - x[k]) * (i + h - x[k]) * (i + h - x[k]);
				Integ2 += h * (y1 + y2) / 2;
			}
		}
	double eps1 = (abs((Integ1 - Integ2) / Integ1));
	cout << "Trapezoid method:\nIntegral = " << Integ2 << endl << "Infelicity = " <<eps1<< endl<<endl;
	fout.open("ans1.dat");
	fout << "Trapezoid method:\nIntegral = " << Integ2 << endl << "Infelicity = " << eps1;
	fout.close();
	//Simpson method
	h = 0.05;
	k = 1;
	Integ1 = 0.0, Integ2 = 0.0;
	k = 1;
	for (double i = x[0]; i <= x[SIZE]; k++) {
		for (; i < x[k] + h; i = i + h) {
			double y1 = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
			double y2 = a[k - 1] + b[k - 1] * (i + h - x[k]) + c[k] * (i + h - x[k]) * (i + h - x[k]) + d[k - 1] * (i + h - x[k]) * (i + h - x[k]) * (i + h - x[k]);
			double y3 = a[k - 1] + b[k - 1] * (i + (h / 2) - x[k]) + c[k] * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]) + d[k - 1] * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]);
			Integ1 += h * (y1 + 4 * y3 + y2) / 6;
		}
	}
	h = h / 2;
	k = 1;
	for (double i = x[0]; i <= x[SIZE]; k++) {
		for (; i < x[k] + h; i = i + h) {
			double y1 = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
			double y2 = a[k - 1] + b[k - 1] * (i + h - x[k]) + c[k] * (i + h - x[k]) * (i + h - x[k]) + d[k - 1] * (i + h - x[k]) * (i + h - x[k]) * (i + h - x[k]);
			double y3 = a[k - 1] + b[k - 1] * (i + (h / 2) - x[k]) + c[k] * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]) + d[k - 1] * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]) * (i + (h / 2) - x[k]);
			Integ2 += h * (y1 + 4 * y3 + y2) / 6;
		}
	}
	double eps2 = (abs((Integ1 - Integ2) / Integ1));
	cout << "Simpson method:\nIntegral = " << Integ2 << endl << "Infelicity = " <<eps2<< endl<<endl;
	fout.open("ans2.dat");
	fout << "Simpson method:\nIntegral = " << Integ2 << endl << "Infelicity = " << eps2;
	fout.close();

	//Rectangle method
	h = 0.05;
	k = 1;
	Integ1 = 0.0, Integ2 = 0.0;
	k = 1;
	for (double i = x[0]; i <= x[SIZE]; k++) {
		for (; i < x[k] + h; i = i + h) {
			double y = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
			Integ1 += y * h;
		}
	}
	h = h / 2;
	k = 1;
	for (double i = x[0]; i <= x[SIZE]; k++) {
		for (; i < x[k] + h; i = i + h) {
			double y = a[k - 1] + b[k - 1] * (i - x[k]) + c[k] * (i - x[k]) * (i - x[k]) + d[k - 1] * (i - x[k]) * (i - x[k]) * (i - x[k]);
			Integ2 += y * h;
		}
	}
	double eps3 = (abs((Integ1 - Integ2) / Integ1));
	cout << "Rectangle method:\nIntegral = " << Integ1 << endl <<"Infelicity = "<<eps3<<endl<< endl;
	fout.open("ans3.dat");
	fout << "Rectangle method:\nIntegral = " << Integ1 << endl << "Infelicity = " << eps3;
	fout.close();
	//system("python3 Lab5.py");
	
	return 0;
}
