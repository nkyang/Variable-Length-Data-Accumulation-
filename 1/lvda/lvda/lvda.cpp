// lvda.cpp : 定义 DLL 应用程序的导出函数。
//

# include "stdafx.h"
# include "mex.h"
# include "matrix.h"
# include "omp.h"

# define N 10000

void cumSignal(short *newSignal, 
	short  *signal, 
	int segmentLen, 
	int totalLen, 
	int flag) 
{
	// flag = -1 多普勒频移为负
	// flag = 0 多普勒频移为0
	// flag = 1 多普勒频移为正
	int newSignalIdx = 0;
	int segmentLenIdx = 0;
	int segmentHalfLen = segmentLen/2;
	int i = 0;
	if (flag == 0) 
	{
		for (i = 0; i < totalLen; i++)
		{
			if (newSignalIdx >= 10000) 
			{
				newSignalIdx = 0;
			}
			newSignal[newSignalIdx] += signal[i];
		}
	}
	else if (flag < 0)
	{
		for (i = 0; i < totalLen; i++)
		{
			if (newSignalIdx >= 10000)
			{
				newSignalIdx = 0;
			}
			if (segmentLenIdx == segmentHalfLen) 
			{
				i++;
				segmentLenIdx++;
			}
			if (segmentLenIdx > segmentLen)
			{
				segmentLenIdx = 0;
			}
			newSignal[newSignalIdx] += signal[i];
			segmentLenIdx++;
			newSignalIdx++;
		}
	}
	else 
	{
		for (i = 0; i < totalLen; i++) 
		{
			if (newSignalIdx >= 10000)
			{
				newSignalIdx = 0;
			}
			if (segmentLenIdx == segmentHalfLen) {
				newSignal[newSignalIdx] += signal[i];
				newSignalIdx++;
			}
			if (segmentLenIdx > segmentLen)
			{
				segmentLenIdx = 0;
			}
			newSignal[newSignalIdx] += signal[i];
			segmentLenIdx++;
			newSignalIdx++;
		}
	}
	//mexPrintf("\n %i", newSignalIdx);
}
void mexFunction(int nlhs,
	mxArray *plhs[],
	int nrhs,
	const mxArray *prhs[])
{
	if (nrhs != 4)
		mexErrMsgIdAndTxt("MATLAB:vlda:invalidNumInputs", "Invalid number of input arguments");
	if (nlhs != 1)
		mexErrMsgIdAndTxt("MATLAB:vlda:invalidNumOutputs", "Invalid number of output arguments");
	short *signal,*newSignal;
	int *totalLen;
	int *segmentLen;
	int *flag;
	int freqNum = mxGetNumberOfElements(prhs[1]);
	const mwSize dims[] = {10000,freqNum};
	signal = mxGetInt16s(prhs[0]);
	segmentLen = mxGetInt32s(prhs[1]);
	totalLen = mxGetInt32s(prhs[2]);
	flag = mxGetInt32s(prhs[3]);
	plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS, mxREAL);
	newSignal = mxGetInt16s(plhs[0]);
#pragma omp parallel for
	for (int i = 0; i < freqNum; i++) {
		cumSignal(newSignal + 10000 * i, signal, segmentLen[i], totalLen[i], flag[i]);
	}
}