#include <iostream>
#include <valarray>
#include <cmath>
#include <vector>

using namespace std;

double f(valarray<double> x){ 
    return -3*pow(x[0],2) - 2*pow(x[1],2) + x[0]*x[1] + 5*x[0] - pow(x[2],2) - x[1]*x[2] + 6*x[2]; 
}

void Results(valarray<double> x){
    cout << "*************************************************" << endl;
    cout << "x = { ";
    for(auto coord : x) cout << coord << " ";
    cout << "}" << endl;
    cout << "f(x) = " << f(x) << endl;
}
valarray<double> grad_f(valarray<double> x) { 
    valarray<double> ans = {-6*x[0] + x[1] + 5, -4*x[1] + x[0] - x[2], -2*x[2] - x[1] + 6};
    return ans;
}

// vector<valarray<double>> GuesseIt(valarray <double> x){
//     vector<valarray<double>> Matrix(x.size());
//     Matrix[0] = {-6, 1, 0};
//     Matrix[1] = {1, -4, -1};
//     Matrix[2] = {0, -1, -2};
// }


double norm(valarray<double> v){ return sqrt((v * v).sum()); }

void GradSpusk(valarray<double> x, double epsilon, bool IsMax){
    double t = 1, c = 1.0/3;
    auto cmp = ((!IsMax) ? [](double x, double y)
                    { return x > y; }
                         : [](double x, double y)
                    { return x <= y; });
    while(norm(grad_f(x)) >= epsilon){
        Results(x);
        if(cmp(f(x), f(x + t*grad_f(x)))){
            x += t*grad_f(x);
        } else {
            t *= c;
        }
    }
}

double GetT(double a, double b, valarray<double> arg, double epsilon, bool IsMax){
    auto cmp = ((!IsMax) ? [](double x, double y)
                    { return x > y; }
                         : [](double x, double y)
                    { return x <= y; });
    double mid;
    while(a + epsilon < b){ 
        mid = (a+b)/2;
        if(cmp(f(arg - (mid - epsilon)*grad_f(arg)), f(arg - (mid + epsilon)*grad_f(arg)))){
            b = mid;
        } else{
            a = mid;
        }
    }
    return mid;
}

void NaiskorSpusk(valarray<double> x, double epsilon, bool IsMax){
    double t = 1, c = 1.0/3;
    auto cmp = ((!IsMax) ? [](double x, double y)
                    { return x > y; }
                         : [](double x, double y)
                    { return x <= y; });
    while(norm(grad_f(x)) >= epsilon){
        Results(x);
        t = GetT(0, 100, x, 0.01, IsMax);
        x += t*grad_f(x);
    }
}

void GetOneDimExtr(double diap, valarray<double>& mid, size_t i, double epsilon, bool IsMax){
    auto cmp = ((!IsMax) ? [](double x, double y)
                    { return x > y; }
                         : [](double x, double y)
                    { return x <= y; });

    double lefter = mid[i] - diap, righter = mid[i] + diap;
    valarray<double> l_mid, r_mid;
    while(lefter + epsilon < righter){ // + epsilon to avoid assplosions
        mid[i] = (lefter + righter) / 2;
        l_mid = mid; l_mid[i] -= epsilon;
        r_mid = mid; r_mid[i] += epsilon;
        if(cmp(f(r_mid), f(l_mid))){
            righter = mid[i];
        } else{
            lefter = mid[i];
        }
    }
    cout << lefter << ' ' << righter << endl;
}

void PokoordSpusk(valarray<double> x, double epsilon, bool IsMax){
    valarray<double> dx(0.1, 3);
    while ((dx * dx).sum() >= epsilon) {
        for (auto i = 0; i < x.size(); ++i){
            Results(x);
            dx[i] = x[i];
            GetOneDimExtr(100, x, i, epsilon, IsMax);
            dx[i] = x[i] - dx[i];
            Results(x);
        }
    }
    Results(x);
}

void Newton(){

}

void SoprPerem(){

}


int main(){
    valarray<double> x = {100, 100, 100};
    bool IsMax = true;
    //vector<valarray<double>> Guesse = {{-6, 1, 0}, {1, -4, -1}, {0, -1, -2}};
    NaiskorSpusk(x, 0.01, IsMax);
}