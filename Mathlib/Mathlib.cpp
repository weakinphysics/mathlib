#include <bits/stdc++.h>
#include <iostream>

using namespace std;


int findLargestSetBit(int n){
    int b = 0;
    int count = 0;
    while(n){
        if(n & 1) b = count;
        count++;
        n = n >> 1;
    }
    return b;
}


int function(vector <int> &nums){
    // find position of ceil 
    int n = nums.size();
    vector <int> positions(n+1);
    for(int i = 0; i < n; i++){
        positions[nums[i]] = i;
    }
    // bhai jugaad hai bohot galat 
    // create a binary search tree and maintain while traversing the array 
    // update the nodes to contain the data you want 
    // return count 
    

}


int main(){
    int t,n;
    cin>>t;
    while(t--){
        cin>>n;
        int lsb = findLargestSetBit(n);
        cout<< ((1 << lsb) - 1) << endl;
    }
    return 0;
}