#include "src/naive.h"
#include "src/formats.h"

#include <iostream>
#include <string>
using namespace std;

int main(int argc, char const *argv[])
{

    string torsoL = "./data/torso1/torso1.mtx";
    string torsob = "./data/torso1/b_for_torso1.mtx";

    auto coo = COO<float>(torsoL.c_str());
    cout<<"COO: "<<coo.m<<","<<coo.n<<" nnz:"<<coo.nnz<<'\n';

    return 0;
}

