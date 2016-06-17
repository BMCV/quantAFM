
#include "opencvmex.hpp"

using namespace cv;


void checkInputs(int nrhs, const mxArray *prhs[])
{
  if (nrhs != 1) 
  {
    mexErrMsgTxt("Incorrect number of inputs. Function expects 1 input.");
  }
  
  if ( mxGetNumberOfDimensions(prhs[0]) > 2 )
  {
    mexErrMsgTxt("Incorrect number of dimensions. Input must be a matrix.");
  }

  if ( !mxIsUint8(prhs[0]) )
  {
    mexErrMsgTxt("Input must be UInt8.");
  }
}



///////////////////////////////////////////////////////////////////////////
// Main entry point to a MEX function
///////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cv::Mat test, test2;
  cv::Mat outputImg;
  // Check inputs to mex function
  checkInputs(nrhs, prhs);
  
  // convert mxArray input (expected to be uint8 array) into CV type
  cv::Ptr< cv::Mat > inputImg = ocvMxArrayToImage_uint8(prhs[0],true);
  test.create(inputImg->size(), CV_8U);
  
  cv::cvtColor( *inputImg, test, CV_GRAY2BGR);
  
  cv::cvtColor(test, test2, CV_BGR2GRAY);
  
  // Allocate output matrix
  int outRows = inputImg->rows;
  int outCols = inputImg->cols;
  //might not be correct OpenCV type for uint8
  //cv::Mat grayImg( (int)outRows, (int)outCols, CV_8UC1 );
  
  /*
   * opencv built-in function for denoising
   *
   * src: 8-bit image
   * target: must be similar to src
   * filter strength: float, higher values do no benefit in our case
   * templateWindowSize: int
   * searchWindowSize: int, too high does not help in our case
   *
   */
  cv::fastNlMeansDenoising(test2, outputImg,  2, 7, 17);

  // Put data back into output MATLAB array
  plhs[0] = ocvMxArrayFromImage_uint8(outputImg);

}
