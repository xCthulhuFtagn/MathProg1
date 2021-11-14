#include <iostream>
#include <valarray>
#include <cmath>
//#include <functional>

using namespace std;

double f(valarray<double> x){ 
    return -3*pow(x[0],2) - 2*pow(x[1],2) + x[0]*x[1] + 5*x[0] - pow(x[2],2) - x[1]*x[2] + 6*x[2]; 
}
valarray<double> grad_f(valarray<double> x) { 
    valarray<double> ans = {-6*x[0] + x[1] + 5, -4*x[1] + x[0] - x[2], -2*x[2] - x[1] + 6};
    return ans;
}


double norm(valarray<double> v){
    double ans = 0;
    for(auto el : v) ans += pow(el, 2);
    return sqrt(ans);
}

void GradSpusk(valarray<double> x, double epsilon){
    double t = 1, c = 1.0/3;
    while(norm(grad_f(x)) >= epsilon){
        cout << "*************************************************" << endl;
        cout << "x = ";
        for(auto el : x) cout << el << " ";
        cout << endl;
        cout << "f(x) = " << f(x) << endl;
        cout << "norm(grad_f(x)) = " << norm(grad_f(x)) << endl;
        if(f(x) > f(x - t*grad_f(x))){
            x -= t*grad_f(x);
        } else {
            t *= c;
        }
    }
}

void NaiskorSpusk(){

}

void PokoordSpusk(){

}

void Newton(){

}

void SoprPerem(){

}


int main(){
    valarray<double> x = {100, 100, 100};
    GradSpusk(x, 0.1);
}