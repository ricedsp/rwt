#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

#include "dwt_init.h"


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])

{
  dwtInit(nlhs,plhs,nrhs,prhs,NORMAL_DWT);
}


void thisIsATest() {
	mxArray *plhs[2];
	const mxArray *prhs[4];
	mexFunction(0,plhs,0,prhs);

	ASSERTM("start writing tests", false);
}

void runSuite(){
	cute::suite s;
	//TODO add your test here
	s.push_back(CUTE(thisIsATest));
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "The Suite");
}

int main(){
    runSuite();
    return 0;
}



