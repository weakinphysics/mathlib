#include <math.h>
#include <assert.h>
#include <bits/stdc++.h>

using namespace std;

double mean(vector <double> &nums){
    double sum = 0;
    int n = nums.size();
    for(int i = 0; i < n; i++) sum += nums[i];
    return sum/n;
}

double variance(vector <double> &nums){
    // standard ensemble calculation for the full population
    int n = nums.size();
    double sum = 0;
    for(int i = 0; i < n; i++) sum += nums[i]*nums[i];
    double mu = mean(nums);
    sum -= n*mu*mu;
    return sum/n;
}

vector <double> cdf(vector <double> &pdf){
    double sum = 0;
    int n = pdf.size();
    vector <double> response;
    for(int i = 0; i < n; i++){
        sum += pdf[i];
        response.push_back(sum);
    }
    return response;
}

double expectation(vector <double> &pdf, vector <double> &nums){
    int n = nums.size();
    assert(pdf.size() == n);
    double mu = 0;
    for(int i = 0; i < n; i++){
        mu += pdf[i]*nums[i];
    }
    return mu;
}

double bin_pow(double number, int power){
    double result = 1.00;
    while(power){
        if(power & 1) result *= number;
        number = number*number;
        power = power >> 1;
    }
    return result;
}

long long calcFactorial()

vector <double> binomial(double p_success, int trial_count){
    // a lot of precision might be lost in this calculation
    for(int i = 0; i < trial_count; i++){

    }
}


class Solution {
private:
    int selectPlace(int number, int power){
        int eval = pow(10, power);
        return ((number/eval)%10);
    }
    
    double mean(vector <int> &nums, int place){
        int n = nums.size();
        for(int i = 0; i < n; i++){
            nums
        }
    }
public:
    long long minimumCost(vector<int>& nums) {
        long long int cost = 0;
        int n = nums.size();
        int magicNumber = 0;
        // for each decimal place calculate the average 
        for(int i = 0; i < n; i++){
            for(int j = 0; j < 10; j++){
                
            }
        }
    }
};

